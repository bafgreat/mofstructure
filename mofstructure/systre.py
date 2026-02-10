#!/usr/bin/env python3
"""
systre.py
=========

Run Systre (Gavrog) on CGD (PERIODIC_GRAPH) files to identify RCSR topology.

Supports:
- CGD file path
- CGD text (string containing "PERIODIC_GRAPH")
- Structure file path (CIF, etc.) -> generate CGD first via TopologyExtractor
- Folder input (batch)
- ASE Atoms object
- pymatgen Structure object (optional dependency)

Also checks for Java and provides OS-specific install guidance.

Author: Dr. Dinga Wonanke
Status: production
"""

from __future__ import annotations

import os
import sys
import shutil
import tempfile
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple, Union

# importlib.resources for packaged jar/arc access
try:
    from importlib.resources import files as ir_files  # py>=3.9
except Exception:  # pragma: no cover
    ir_files = None  # type: ignore

from ase.atoms import Atoms
from ase.io import read as ase_read

from mofstructure.topology import TopologyExtractor


CGDLike = Union[str, Path]
InputLike = Union[str, Path, Atoms, "PymatgenStructure"]  # pymatgen is optional


def _is_cgd_text(s: str) -> bool:
    return ("PERIODIC_GRAPH" in s) and ("EDGES" in s)


def _java_install_hints() -> str:
    plat = sys.platform.lower()
    # Keep this simple + robust (no web calls); Adoptium is the safest general recommendation.
    if plat.startswith("linux"):
        return (
            "Java not found.\n\n"
            "Install a JRE/JDK (recommended: Temurin/OpenJDK):\n"
            "  - Adoptium (Temurin): https://adoptium.net/\n"
            "Or via your package manager, e.g.:\n"
            "  - Ubuntu/Debian: sudo apt-get install default-jre\n"
            "  - Fedora/RHEL:   sudo dnf install java-17-openjdk\n"
            "  - Arch:          sudo pacman -S jre-openjdk\n"
        )
    if plat == "darwin":
        return (
            "Java not found.\n\n"
            "Install a JRE/JDK (recommended: Temurin/OpenJDK):\n"
            "  - Adoptium (Temurin): https://adoptium.net/\n"
            "Or via Homebrew:\n"
            "  - brew install --cask temurin\n"
        )
    if plat.startswith("win"):
        return (
            "Java not found.\n\n"
            "Install a JRE/JDK (recommended: Temurin/OpenJDK):\n"
            "  - Adoptium (Temurin): https://adoptium.net/\n"
            "After install, ensure 'java' is on PATH (restart terminal)."
        )
    return (
        "Java not found.\n\n"
        "Install a JRE/JDK:\n"
        "  - Adoptium (Temurin): https://adoptium.net/\n"
    )


def find_java(java: str = "java") -> str:
    """
    Return a usable java executable path/name. Raises RuntimeError if not found.
    """
    exe = shutil.which(java) if os.path.sep not in java else (java if os.path.exists(java) else None)
    if not exe:
        raise RuntimeError(_java_install_hints())
    return exe


def _resource_path(rel: str) -> str:
    """
    Return absolute filesystem path to a packaged resource (jar/arc).
    Works for normal installs and editable installs. For zipped installs,
    importlib.resources gives a virtual path; we materialize to a temp file if needed.
    """
    # Common case: source tree / editable install -> direct file path exists
    direct = Path(__file__).resolve().parent / rel
    if direct.exists():
        return str(direct)

    if ir_files is None:
        raise RuntimeError(f"Cannot locate resource {rel}. importlib.resources not available.")

    pkg_root = ir_files("mofstructure")
    candidate = pkg_root.joinpath(rel)

    # If candidate is a real file on disk, return it
    try:
        cand_path = Path(str(candidate))
        if cand_path.exists():
            return str(cand_path)
    except Exception:
        pass

    # Otherwise, materialize to temp (rare, but safe)
    data = candidate.read_bytes()
    suffix = Path(rel).suffix
    fd, tmp = tempfile.mkstemp(prefix="mofstructure_resource_", suffix=suffix)
    os.close(fd)
    with open(tmp, "wb") as f:
        f.write(data)
    return tmp


def systre_command(
    *,
    java: str = "java",
    xmx: str = "1024m",
    jar_relpath: str = "bin/Systre-19.6.0.jar",
    rcsr_relpath: str = "db/RCSRnets-2019-06-01.arc",
) -> List[str]:
    """
    Build the base Systre command list.
    """
    java_exe = find_java(java)
    jar = _resource_path(jar_relpath)
    rcsr = _resource_path(rcsr_relpath)
    return [
        java_exe,
        f"-Xmx{xmx}",
        "-cp",
        jar,
        "org.gavrog.apps.systre.SystreCmdline",
        rcsr,
    ]


def run_systre_on_cgd_path(
    cgd_path: Union[str, Path],
    *,
    java: str = "java",
    timeout_s: int = 30,
    xmx: str = "1024m",
) -> subprocess.CompletedProcess:
    """
    Run Systre on a CGD file path; returns CompletedProcess.
    """
    cmd = systre_command(java=java, xmx=xmx) + [str(cgd_path)]
    return subprocess.run(
        cmd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=timeout_s,
    )


def parse_systre_topology(stdout: str) -> str:
    """
    Parse Systre stdout and return:
      - RCSR net name (e.g. "pcu", "dia", ...)
      - or one of: "UNKNOWN", "ERROR", "TIMEOUT", "MISMATCH"
    This follows the logic you showed from mofid but keeps it self-contained.
    """
    topologies: List[str] = []
    current_component = 0
    topology_line = False
    repeat_line = False

    for raw_line in stdout.splitlines():
        line = raw_line.strip()

        if topology_line:
            topology_line = False
            parts = line.split()
            # Expected: "Name: pcu"
            if len(parts) >= 2 and parts[0] == "Name:":
                topologies.append(parts[1])
            else:
                return "ERROR"
            continue

        if repeat_line:
            repeat_line = False
            parts = line.split()
            if not parts or parts[0] != "Name:":
                return "ERROR"
            # Expected: "Name: refcode_clean_component_2"
            comps = line.split("_")
            if len(comps) >= 2 and comps[-2] == "component":
                try:
                    idx = int(comps[-1]) - 1  # Systre is 1-indexed
                    topologies.append(topologies[idx])
                except Exception:
                    return "ERROR"
            else:
                return "ERROR"
            continue

        if "ERROR" in line:
            return "ERROR"
        if "Structure was found in archive" in line:
            topology_line = True
            continue
        if line == "Structure is new for this run.":
            topologies.append("UNKNOWN")
            continue
        if line == "Structure already seen in this run.":
            repeat_line = True
            continue
        if "Processing component " in line:
            current_component += 1
            # Basic sanity: should extract one topology per component
            if len(topologies) != (current_component - 1):
                # Don't hard-fail; some Systre versions differ slightly
                pass
            continue

    if not topologies:
        return "ERROR"
    first = topologies[0]
    for t in topologies[1:]:
        if t != first:
            return "MISMATCH"
    return first


@dataclass
class SystreResult:
    topology: str
    stdout: str
    stderr: str
    cgd_path: Optional[str] = None


def _atoms_from_pymatgen(obj) -> Atoms:
    """
    Convert pymatgen Structure to ASE Atoms.
    """
    try:
        from pymatgen.io.ase import AseAtomsAdaptor  # type: ignore
    except Exception as exc:  # pragma: no cover
        raise ImportError(
            "pymatgen is not installed, but a pymatgen Structure was provided. "
            "Install pymatgen or pass an ASE Atoms object instead."
        ) from exc
    return AseAtomsAdaptor.get_atoms(obj)


def _make_cgd_from_input(
    x: InputLike,
    *,
    method: str = "sbus",
    name: str = "net",
    top_k_regions: int = 1,
    connect_mode: str = "clique",
    dedup_parallel_edges: bool = True,
    dedup_mode: str = "shift",
) -> str:
    """
    Accept Atoms, pymatgen Structure, or a structure file path; produce CGD text.
    """
    if isinstance(x, Atoms):
        topo = TopologyExtractor(ase_atoms=x)
        return topo.build_cgd(
            method=method,
            name=name,
            top_k_regions=top_k_regions,
            connect_mode=connect_mode,
            dedup_parallel_edges=dedup_parallel_edges,
            dedup_mode=dedup_mode,
        )

    # If it's a path-like
    if isinstance(x, (str, Path)):
        p = Path(x)
        if p.exists() and p.is_file():
            topo = TopologyExtractor(filename=str(p))
            return topo.build_cgd(
                method=method,
                name=name,
                top_k_regions=top_k_regions,
                connect_mode=connect_mode,
                dedup_parallel_edges=dedup_parallel_edges,
                dedup_mode=dedup_mode,
            )
        raise FileNotFoundError(f"Input path not found: {p}")

    # Assume pymatgen Structure-like
    return _make_cgd_from_input(
        _atoms_from_pymatgen(x),
        method=method,
        name=name,
        top_k_regions=top_k_regions,
        connect_mode=connect_mode,
        dedup_parallel_edges=dedup_parallel_edges,
        dedup_mode=dedup_mode,
    )


def identify_topology(
    x: Union[CGDLike, Atoms, "PymatgenStructure"],
    *,
    input_is_cgd: Optional[bool] = None,
    method: str = "all_node",
    name: str = "net",
    top_k_regions: int = 1,
    connect_mode: str = "clique",
    dedup_parallel_edges: bool = True,
    dedup_mode: str = "shift",
    # systre options
    java: str = "java",
    timeout_s: int = 30,
    xmx: str = "1024m",
    keep_tmp: bool = False,
) -> SystreResult:
    """
    High-level API:
      - If x is a CGD file path -> run directly
      - If x is CGD text -> write temp CGD and run
      - If x is Atoms / pymatgen Structure / structure file -> generate CGD first then run

    Returns SystreResult with parsed topology plus stdout/stderr.
    """
    tmp_path: Optional[str] = None

    try:
        # Decide how to interpret x
        if isinstance(x, (str, Path)):
            s = str(x)
            p = Path(s)
            if input_is_cgd is True:
                # treat as cgd path or cgd text
                if p.exists() and p.is_file():
                    proc = run_systre_on_cgd_path(p, java=java, timeout_s=timeout_s, xmx=xmx)
                    topo = parse_systre_topology(proc.stdout)
                    return SystreResult(topo, proc.stdout, proc.stderr, cgd_path=str(p))
                if _is_cgd_text(s):
                    fd, tmp_path = tempfile.mkstemp(prefix="mofstructure_", suffix=".cgd")
                    os.close(fd)
                    Path(tmp_path).write_text(s, encoding="utf-8")
                    proc = run_systre_on_cgd_path(tmp_path, java=java, timeout_s=timeout_s, xmx=xmx)
                    topo = parse_systre_topology(proc.stdout)
                    return SystreResult(topo, proc.stdout, proc.stderr, cgd_path=tmp_path)
                raise ValueError("input_is_cgd=True but input is neither a CGD file nor CGD text.")

            # auto-detect:
            if p.exists() and p.is_file() and p.suffix.lower() == ".cgd":
                proc = run_systre_on_cgd_path(p, java=java, timeout_s=timeout_s, xmx=xmx)
                topo = parse_systre_topology(proc.stdout)
                return SystreResult(topo, proc.stdout, proc.stderr, cgd_path=str(p))

            if _is_cgd_text(s):
                fd, tmp_path = tempfile.mkstemp(prefix="mofstructure_", suffix=".cgd")
                os.close(fd)
                Path(tmp_path).write_text(s, encoding="utf-8")
                proc = run_systre_on_cgd_path(tmp_path, java=java, timeout_s=timeout_s, xmx=xmx)
                topo = parse_systre_topology(proc.stdout)
                return SystreResult(topo, proc.stdout, proc.stderr, cgd_path=tmp_path)

            # otherwise treat as structure file
            if p.exists() and p.is_file():
                cgd_text = _make_cgd_from_input(
                    p,
                    method=method,
                    name=name,
                    top_k_regions=top_k_regions,
                    connect_mode=connect_mode,
                    dedup_parallel_edges=dedup_parallel_edges,
                    dedup_mode=dedup_mode,
                )
                fd, tmp_path = tempfile.mkstemp(prefix="mofstructure_", suffix=".cgd")
                os.close(fd)
                Path(tmp_path).write_text(cgd_text, encoding="utf-8")
                proc = run_systre_on_cgd_path(tmp_path, java=java, timeout_s=timeout_s, xmx=xmx)
                topo = parse_systre_topology(proc.stdout)
                return SystreResult(topo, proc.stdout, proc.stderr, cgd_path=tmp_path)

            raise FileNotFoundError(f"Path not found: {p}")

        # ASE Atoms / pymatgen Structure
        cgd_text = _make_cgd_from_input(
            x,
            method=method,
            name=name,
            top_k_regions=top_k_regions,
            connect_mode=connect_mode,
            dedup_parallel_edges=dedup_parallel_edges,
            dedup_mode=dedup_mode,
        )
        fd, tmp_path = tempfile.mkstemp(prefix="mofstructure_", suffix=".cgd")
        os.close(fd)
        Path(tmp_path).write_text(cgd_text, encoding="utf-8")
        proc = run_systre_on_cgd_path(tmp_path, java=java, timeout_s=timeout_s, xmx=xmx)
        topo = parse_systre_topology(proc.stdout)
        return SystreResult(topo, proc.stdout, proc.stderr, cgd_path=tmp_path)

    except subprocess.TimeoutExpired as exc:
        return SystreResult("TIMEOUT", stdout=getattr(exc, "stdout", "") or "", stderr=getattr(exc, "stderr", "") or "", cgd_path=tmp_path)

    finally:
        if tmp_path and (not keep_tmp):
            try:
                os.remove(tmp_path)
            except Exception:
                pass


def identify_topology_batch(
    folder: Union[str, Path],
    *,
    patterns: Tuple[str, ...] = (".cgd", ".cif", ".vasp", ".poscar"),
    recursive: bool = True,
    **kwargs,
) -> Dict[str, SystreResult]:
    """
    Batch topology identification for a folder of CGD or structure files.
    Returns mapping: filepath -> SystreResult.
    """
    folder = Path(folder)
    if not folder.exists() or not folder.is_dir():
        raise NotADirectoryError(f"Not a folder: {folder}")

    results: Dict[str, SystreResult] = {}
    it = folder.rglob("*") if recursive else folder.glob("*")
    for p in it:
        if not p.is_file():
            continue
        if p.suffix.lower() not in {x.lower() for x in patterns}:
            continue
        try:
            res = identify_topology(p, **kwargs)
        except Exception as exc:
            results[str(p)] = SystreResult("ERROR", stdout="", stderr=str(exc), cgd_path=str(p) if p.suffix.lower() == ".cgd" else None)
        else:
            results[str(p)] = res
    return results
