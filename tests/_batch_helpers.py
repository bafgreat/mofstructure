"""Tasks that misbehave on purpose, for tests/test_batch.py."""
import ctypes
import os
import subprocess
import sys
import time


def behave(path, **kwargs):
    mode = os.path.basename(path)
    if mode.startswith("ok"):
        return {"mode": mode, **kwargs}
    if mode.startswith("raise"):
        raise ValueError("deliberate failure")
    if mode.startswith("abort"):
        os.abort()
    if mode.startswith("segfault"):
        ctypes.string_at(0)
    if mode.startswith("hang"):
        time.sleep(600)
    if mode.startswith("child"):
        # a long child of the task, as zeo++ is for porosity; it records its
        # pid so the test can check it did not outlive the worker
        child = subprocess.Popen([sys.executable, "-c",
                                  "import time; time.sleep(600)"])
        with open(kwargs["pid_file"], "w") as handle:
            handle.write(str(child.pid))
        child.wait()
    return None
