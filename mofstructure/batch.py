#!/usr/bin/env python3
'''
Run an analysis over many structures so that one structure cannot stop it.

A long batch over a structure database fails in ways that have nothing to do
with the structure being analysed at the time it fails:

    - compiled code (Open Babel, zeo++, numpy) can crash with a segmentation
      fault or abort, which ends the Python process and everything in it;
    - a structure can run for hours, for example a very large cell in the
      deconstruction, the canonical form or the embedding, and a signal based
      timeout cannot interrupt compiled code;
    - memory can grow until the scheduler kills the whole job;
    - results kept in memory and written at the end, or rewritten in full after
      every structure, are lost or corrupted when the job is killed.

`run_batch` therefore runs every structure in a worker process that the parent
can kill. Each worker is started with the spawn method, which is safe on macOS
and on Linux, handles structures one at a time and is replaced after a fixed
number of structures to release memory. The parent enforces a wall-clock
limit per structure, records a crash or a timeout as the result for that
structure and starts a new worker. Every result is appended as one line to a
json-lines checkpoint file and flushed to disk at once, so a job that is
killed loses at most the structures in flight, and repeating the same command
continues where it stopped.

A run can also be split across the tasks of a cluster array job with
`shard`, each task writing its own checkpoint file, and the parts are joined
by `consolidate`.

**checkpoint record:**
    One json object per line::

        {"name": "ABAVIJ", "file": ".../ABAVIJ.cif", "status": "ok",
         "elapsed_s": 1.83, "data": {...}}

    `status` is "ok", "error" (the analysis raised; `detail` holds the
    exception), "timeout" (killed after the wall-clock limit) or "crashed"
    (the worker died; `detail` names the signal or exit code).
'''
from __future__ import annotations

import importlib
import json
import logging
import multiprocessing as mp
import os
import signal
import sys
import time
import traceback
from collections.abc import Callable, Iterable, Sequence
from multiprocessing.connection import wait
from pathlib import Path

logger = logging.getLogger(__name__)

CHECKPOINT_DIR = "_progress"


# ---------------------------------------------------------------------------
# Worker side
# ---------------------------------------------------------------------------

def _limit_memory(memory_gb: float | None) -> None:
    '''
    Cap the address space of the worker so that a runaway structure fails in
    its own process instead of drawing the scheduler's out-of-memory killer
    onto the whole job. Linux enforces the limit; macOS does not, and the call
    is skipped there.
    '''
    if not memory_gb:
        return
    try:
        import resource

        limit = int(memory_gb * 1024 ** 3)
        resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    except (ImportError, ValueError, OSError):
        logger.debug("memory limit not supported on this platform")


def _resolve(task: str | Callable) -> Callable:
    '''Import a task given as "module:function", or return a callable.'''
    if callable(task):
        return task
    module, _, name = task.partition(":")
    return getattr(importlib.import_module(module), name)


def _worker(conn, task, kwargs, memory_gb, quiet) -> None:
    '''
    Worker loop: receive (name, path), run the task, send the result back.

    A None message ends the worker. Anything the task raises is returned as
    an error record, so only a crash in compiled code can end the process
    early, and that is detected by the parent.
    '''
    # The parent decides when to stop, so the worker leaves an interrupt to
    # the parent.
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    # A task can start children of its own, as porosity does for zeo++.
    # Killing only the worker on a timeout would leave such a child running
    # with nothing left to stop it, so the worker leads a process group the
    # parent kills as a whole. This also keeps a Ctrl-C at the terminal away
    # from the workers; the parent stops them on its way out.
    if hasattr(os, "setpgrp"):
        os.setpgrp()
    _limit_memory(memory_gb)
    if quiet:
        import warnings

        # zeo++, voro++ and Open Babel print from compiled code, which no
        # Python setting can silence, and over a large folder that output
        # buries the progress lines of a cluster log. Pointing the worker's
        # file descriptors at /dev/null silences them and the child
        # interpreters they start; errors still reach the checkpoint record.
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, 1)
        os.dup2(devnull, 2)
        warnings.filterwarnings("ignore")
        logging.disable(logging.WARNING)
        try:
            from rdkit import RDLogger

            RDLogger.DisableLog("rdApp.*")
        except ImportError:
            pass
    function = _resolve(task)
    while True:
        try:
            message = conn.recv()
        except EOFError:
            return
        if message is None:
            return
        name, path = message
        start = time.perf_counter()
        try:
            data = function(path, **kwargs)
            reply = {"status": "ok", "data": data}
        except MemoryError:
            reply = {"status": "error", "detail": "MemoryError"}
        except Exception as exc:  # noqa: BLE001 - reported, never swallowed
            reply = {
                "status": "error",
                "detail": f"{type(exc).__name__}: {exc}"[:500],
                "traceback": traceback.format_exc(limit=8)[-2000:],
            }
        reply["elapsed_s"] = round(time.perf_counter() - start, 3)
        try:
            conn.send((name, reply))
        except (TypeError, ValueError, AttributeError) as exc:
            # The task returned something that cannot be pickled.
            conn.send((name, {"status": "error",
                              "detail": f"unpicklable result: {exc}",
                              "elapsed_s": reply["elapsed_s"]}))


# ---------------------------------------------------------------------------
# Checkpoint files
# ---------------------------------------------------------------------------

def shard_suffix(shard: tuple[int, int] | None) -> str:
    '''File suffix for one shard of a split run, empty for a whole run.'''
    if not shard:
        return ""
    index, count = shard
    return f".part-{index:04d}-of-{count:04d}"


def parse_shard(text: str | None) -> tuple[int, int] | None:
    '''
    Parse "i/N" into (i, N) with 0 <= i < N. Also accepts the value of
    SLURM_ARRAY_TASK_ID when given as "slurm/N".
    '''
    if not text:
        return None
    index, _, count = text.partition("/")
    if index.lower() in ("slurm", "auto"):
        index = os.environ.get("SLURM_ARRAY_TASK_ID",
                               os.environ.get("PBS_ARRAY_INDEX", ""))
    index, count = int(index), int(count)
    if not 0 <= index < count:
        raise ValueError(f"shard index {index} is outside 0..{count - 1}")
    return index, count


def select_shard(items: Sequence[tuple[str, str]],
                 shard: tuple[int, int] | None) -> list[tuple[str, str]]:
    '''
    The items one shard is responsible for. Items are sorted by name and dealt
    out in turn, so the split does not depend on the order the files were
    listed in and every shard receives a similar mix of sizes.
    '''
    ordered = sorted(items)
    if not shard:
        return ordered
    index, count = shard
    return ordered[index::count]


def read_checkpoints(paths: Iterable[Path]) -> dict[str, dict]:
    '''
    Read checkpoint files into name -> record, later lines winning. A final
    line cut short by a killed job is ignored.
    '''
    records: dict[str, dict] = {}
    for path in paths:
        if not path.exists():
            continue
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                try:
                    record = json.loads(line)
                except json.JSONDecodeError:
                    continue
                records[record["name"]] = record
    return records


class Checkpoint:
    '''Append-only json-lines writer that flushes every record to disk.'''

    def __init__(self, path: Path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.handle = open(self.path, "a", encoding="utf-8")

    def write(self, record: dict) -> None:
        self.handle.write(json.dumps(record, default=str) + "\n")
        self.handle.flush()
        os.fsync(self.handle.fileno())

    def close(self) -> None:
        self.handle.close()


# ---------------------------------------------------------------------------
# Parent side
# ---------------------------------------------------------------------------

def _describe_exit(code: int | None) -> str:
    if code is None:
        return "worker vanished"
    if code < 0:
        try:
            return f"killed by {signal.Signals(-code).name}"
        except ValueError:
            return f"killed by signal {-code}"
    return f"exit code {code}"


class _Slot:
    '''One worker process and the structure it is working on.'''

    def __init__(self, context, task, kwargs, memory_gb, quiet):
        self.args = (task, kwargs, memory_gb, quiet)
        self.context = context
        self.start_worker()

    def start_worker(self):
        parent, child = self.context.Pipe()
        self.process = self.context.Process(
            target=_worker, args=(child, *self.args), daemon=True)
        self.process.start()
        child.close()
        self.conn = parent
        self.current = None
        self.started = 0.0
        self.done = 0

    def assign(self, item):
        self.current = item
        self.started = time.monotonic()
        self.conn.send(item)

    def kill(self):
        '''Kill the worker together with any child its task started.'''
        try:
            os.killpg(self.process.pid, signal.SIGKILL)
        except (AttributeError, ProcessLookupError, PermissionError):
            # no process groups (Windows), or the group is already gone
            self.process.kill()

    def stop(self, force=False):
        try:
            if force:
                self.kill()
            else:
                self.conn.send(None)
        except (BrokenPipeError, OSError):
            pass
        self.process.join(timeout=10)
        if self.process.is_alive():
            self.kill()
            self.process.join(timeout=10)
        self.conn.close()


def run_batch(
    items: Sequence[tuple[str, str]],
    task: str | Callable,
    checkpoint: Path,
    *,
    kwargs: dict | None = None,
    workers: int = 1,
    timeout: float | None = 3600,
    memory_gb: float | None = None,
    recycle: int = 50,
    skip: Iterable[str] = (),
    quiet: bool = True,
    on_record: Callable[[dict, int, int], None] | None = None,
) -> dict[str, int]:
    '''
    Run `task(path, **kwargs)` on every (name, path) item in isolated workers.

    **parameters:**
        - items: sequence of (name, path)

        - task: "module:function" or a picklable top-level function
            Called as task(path, **kwargs); must return json-serialisable
            data.

        - checkpoint: pathlib.Path
            json-lines file results are appended to.

        - kwargs: python dictionary
            Keyword arguments for the task.

        - workers: int
            Number of worker processes. 0 runs everything in the calling
            process without isolation or time limit, for debugging.

        - timeout: float or None
            Wall-clock seconds per structure before its worker is killed.

        - memory_gb: float or None
            Address-space limit per worker (Linux).

        - recycle: int
            Structures a worker handles before it is replaced.

        - skip: iterable of str
            Names already done, for resuming.

        - quiet: bool
            Silence warnings and library logging in the workers.

        - on_record: callable(record, done, total), optional
            Called after every structure, for progress output.

    **returns:**
        python dictionary
            Count of records by status.
    '''
    kwargs = kwargs or {}
    skip = set(skip)
    queue = [item for item in items if item[0] not in skip]
    total = len(queue)
    tally: dict[str, int] = {}
    if not queue:
        return tally

    if workers == 0:
        return _run_in_process(queue, task, kwargs, checkpoint, tally,
                               on_record)

    context = mp.get_context("spawn")
    writer = Checkpoint(checkpoint)
    slots = [_Slot(context, task, kwargs, memory_gb, quiet)
             for _ in range(max(1, min(workers, total)))]
    files = dict(queue)
    done = 0

    def finish(slot, reply):
        nonlocal done
        name = slot.current[0]
        record = {"name": name, "file": files[name], **reply}
        writer.write(record)
        tally[record["status"]] = tally.get(record["status"], 0) + 1
        done += 1
        slot.current = None
        slot.done += 1
        if on_record:
            on_record(record, done, total)

    try:
        pending = list(reversed(queue))
        while pending or any(slot.current for slot in slots):
            for slot in slots:
                if slot.current is None and pending:
                    if slot.done >= recycle or not slot.process.is_alive():
                        slot.stop()
                        slot.start_worker()
                    item = pending.pop()
                    try:
                        slot.assign(item)
                    except (BrokenPipeError, OSError):
                        slot.stop(force=True)
                        slot.start_worker()
                        slot.assign(item)

            busy = [slot for slot in slots if slot.current]
            ready = wait([slot.conn for slot in busy] +
                         [slot.process.sentinel for slot in busy], timeout=1.0)
            now = time.monotonic()
            for slot in busy:
                if slot.conn in ready or slot.conn.poll():
                    try:
                        _name, reply = slot.conn.recv()
                    except (EOFError, OSError):
                        slot.process.join(timeout=5)
                        reply = {"status": "crashed",
                                 "detail": _describe_exit(
                                     slot.process.exitcode),
                                 "elapsed_s": round(now - slot.started, 3)}
                        finish(slot, reply)
                        slot.stop(force=True)
                        slot.start_worker()
                        continue
                    finish(slot, reply)
                elif not slot.process.is_alive():
                    reply = {"status": "crashed",
                             "detail": _describe_exit(slot.process.exitcode),
                             "elapsed_s": round(now - slot.started, 3)}
                    finish(slot, reply)
                    slot.stop(force=True)
                    slot.start_worker()
                elif timeout and now - slot.started > timeout:
                    reply = {"status": "timeout",
                             "detail": f"exceeded {timeout:g} s",
                             "elapsed_s": round(now - slot.started, 3)}
                    finish(slot, reply)
                    slot.stop(force=True)
                    slot.start_worker()
    finally:
        for slot in slots:
            slot.stop(force=bool(slot.current))
        writer.close()
    return tally


def _run_in_process(queue, task, kwargs, checkpoint, tally, on_record):
    '''
    Run the items in the calling process, for debugging and tests. Exceptions
    are recorded as for a worker, but a crash in compiled code ends the run
    and no time limit is applied.
    '''
    function = _resolve(task)
    writer = Checkpoint(checkpoint)
    try:
        for done, (name, path) in enumerate(queue, start=1):
            start = time.perf_counter()
            try:
                reply = {"status": "ok", "data": function(path, **kwargs)}
            except Exception as exc:  # noqa: BLE001 - recorded
                reply = {"status": "error",
                         "detail": f"{type(exc).__name__}: {exc}"[:500]}
            reply["elapsed_s"] = round(time.perf_counter() - start, 3)
            record = {"name": name, "file": path, **reply}
            writer.write(record)
            tally[record["status"]] = tally.get(record["status"], 0) + 1
            if on_record:
                on_record(record, done, len(queue))
    finally:
        writer.close()
    return tally


def progress_printer(label: str = "") -> Callable[[dict, int, int], None]:
    '''
    Progress callback printing one flushed line per structure with an
    estimate of the time left. Output is flushed so that a log file on a
    cluster shows progress as it happens rather than in buffered blocks.
    '''
    start = time.monotonic()

    def show(record: dict, done: int, total: int) -> None:
        elapsed = time.monotonic() - start
        left = elapsed / done * (total - done) if done else 0
        detail = record.get("detail", "")
        print(f"{label}{done}/{total} {record['name']:<30} "
              f"{record['status']:<8} {record.get('elapsed_s', 0):8.1f} s  "
              f"eta {left / 3600:6.2f} h  {detail[:60]}",
              flush=True)

    return show


def structure_items(folder_or_files: Sequence[str],
                    suffixes: Sequence[str] = (".cif",),
                    naming: str = "stem") -> list[tuple[str, str]]:
    '''
    List (name, path) for the structures in folders or file lists.

    Hidden files, which on macOS include the "._" companions an archive
    leaves behind, are skipped.

    **parameters:**
        - folder_or_files: paths to folders or files

        - suffixes: file suffixes read from a folder

        - naming: "stem" uses the file name without its last suffix;
          "first_dot" keeps the text before the first dot, the convention
          the database commands have always used.
    '''
    paths: list[Path] = []
    for entry in folder_or_files:
        entry = Path(entry)
        if entry.is_dir():
            paths.extend(
                p for p in entry.iterdir()
                if p.is_file() and not p.name.startswith(".")
                and p.suffix.lower() in suffixes
            )
        elif entry.exists():
            paths.append(entry)
    items = []
    for path in paths:
        name = path.stem if naming == "stem" else path.name.split(".")[0]
        items.append((name, str(path)))
    return items


if __name__ == "__main__":  # pragma: no cover
    sys.exit("mofstructure.batch is a library; use the mofstructure commands")
