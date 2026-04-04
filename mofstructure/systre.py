#!/usr/bin/env python3
"""
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
from typing import Dict, List, Optional, Tuple, Union

try:
    from importlib.resources import files as ir_files  # py>=3.9
except Exception:  # pragma: no cover
    ir_files = None  # type: ignore

from ase.atoms import Atoms

from mofstructure.generate_cgd import TopologyExtractor

CGDLike = Union[str, Path]
InputLike = Union[str, Path, Atoms, "PymatgenStructure"]  # pymatgen is optional


def _is_cgd_text(s: str) -> bool:
    '''
    Check whether a string appears to contain CGD text.

    The test is intentionally simple and only checks whether the string
    contains the standard PERIODIC_GRAPH and EDGES blocks expected in
    a Systre-readable CGD file.

    **parameters:**
        s: str
            Input string to test.

    **returns:**
        bool
            True if the string appears to be CGD text, otherwise False.
    '''
    return ("PERIODIC_GRAPH" in s) and ("EDGES" in s)


def _java_install_hints() -> str:
    '''
    Return platform-specific installation hints when Java is not found.

    Systre requires a working Java executable. This helper generates a
    user-friendly installation message depending on the current operating
    system.

    **returns:**
        str
            A formatted help message describing how to install Java.
    '''
    plat = sys.platform.lower()
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
    '''
    Find a usable Java executable.

    The function first checks whether the supplied executable name or path
    exists and is accessible. If Java is not found, a RuntimeError is raised
    with platform-specific installation guidance.

    **parameters:**
        java: str
            Name of the Java executable or an explicit path to it.
            Default is "java".

    **returns:**
        str
            Path or executable name of a working Java command.

    **raises:**
        RuntimeError:
            If no Java executable can be found.
    '''
    exe = shutil.which(java) if os.path.sep not in java else (java if os.path.exists(java) else None)
    if not exe:
        raise RuntimeError(_java_install_hints())
    return exe


def _resource_path(rel: str) -> str:
    '''
    Return the absolute filesystem path of a packaged resource.

    This helper is used to locate packaged resources such as the Systre JAR
    file or the RCSR archive. It works for standard installs, editable installs,
    and environments where importlib.resources returns a virtual path. In the
    latter case, the resource is materialized into a temporary file.

    **parameters:**
        rel: str
            Relative path to the resource inside the package.

    **returns:**
        str
            Absolute filesystem path to the resource.

    **raises:**
        RuntimeError:
            If the resource cannot be located and importlib.resources is not
            available.
    '''
    direct = Path(__file__).resolve().parent / rel
    if direct.exists():
        return str(direct)

    if ir_files is None:
        raise RuntimeError(f"Cannot locate resource {rel}. importlib.resources not available.")

    pkg_root = ir_files("mofstructure")
    candidate = pkg_root.joinpath(rel)

    try:
        cand_path = Path(str(candidate))
        if cand_path.exists():
            return str(cand_path)
    except Exception:
        pass

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
    '''
    Build the command-line call required to run Systre.

    This function resolves the Java executable, locates the packaged Systre
    JAR file and the RCSR archive, and returns the base command list that can
    be passed to `subprocess.run`.

    **parameters:**
        java: str
            Java executable name or explicit path.
            Default is "java".

        xmx: str
            Maximum Java heap size.
            Example: "1024m".

        jar_relpath: str
            Relative path to the packaged Systre JAR file.

        rcsr_relpath: str
            Relative path to the packaged RCSR archive file.

    **returns:**
        list
            Command list suitable for subprocess execution.
    '''
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
    '''
    Run Systre on a CGD file located on disk.

    **parameters:**
        cgd_path: str or Path
            Path to a CGD file.

        java: str
            Java executable name or explicit path.

        timeout_s: int
            Timeout in seconds for the Systre run.

        xmx: str
            Maximum Java heap size.

    **returns:**
        subprocess.CompletedProcess
            The completed subprocess result containing stdout, stderr,
            and the process return code.
    '''
    cmd = systre_command(java=java, xmx=xmx) + [str(cgd_path)]
    return subprocess.run(
        cmd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=timeout_s,
    )


def parse_systre_topology(stdout: str) -> str:
    '''
    Parse the standard output of Systre and extract the identified topology.

    The function looks for net names reported by Systre and returns a single
    consensus topology if all processed components match. If multiple components
    produce different topologies, the function returns "MISMATCH". Other sentinel
    values such as "UNKNOWN" or "ERROR" are also returned where appropriate.

    **parameters:**
        stdout: str
            Standard output produced by a Systre run.

    **returns:**
        str
            One of the following:
                - an RCSR topology name, e.g. "pcu", "dia"
                - "UNKNOWN"
                - "ERROR"
                - "TIMEOUT"
                - "MISMATCH"
    '''
    topologies: List[str] = []
    current_component = 0
    topology_line = False
    repeat_line = False

    for raw_line in stdout.splitlines():
        line = raw_line.strip()

        if topology_line:
            topology_line = False
            parts = line.split()
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
            comps = line.split("_")
            if len(comps) >= 2 and comps[-2] == "component":
                try:
                    idx = int(comps[-1]) - 1
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
    '''
    Container class for the result of a Systre topology identification.

    **attributes:**
        topology: str
            Identified topology name or sentinel value such as "ERROR" or "UNKNOWN".

        stdout: str
            Standard output from the Systre process.

        stderr: str
            Standard error from the Systre process.

        cgd_path: str or None
            Path to the CGD file used for the Systre run, if available.
    '''
    topology: str
    stdout: str
    stderr: str
    cgd_path: Optional[str] = None


def _atoms_from_pymatgen(obj) -> Atoms:
    '''
    Convert a pymatgen Structure object to an ASE Atoms object.

    **parameters:**
        obj:
            A pymatgen Structure object.

    **returns:**
        ASE atoms object

    **raises:**
        ImportError:
            If pymatgen is not installed but a pymatgen Structure was provided.
    '''
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
    method: str = "all_node",
    name: str = "net",
) -> str:
    '''
    Generate CGD text from a supported structure-like input.

    This helper accepts:
        - ASE Atoms objects
        - pymatgen Structure objects
        - structure file paths readable by ASE

    It then calls the new topology module to generate a CGD representation.

    **parameters:**
        x:
            Structure-like input. This can be an ASE atoms object, a pymatgen
            Structure object, or a structure file path.

        method: str
            Topology extraction mode passed to TopologyExtractor.build_cgd().

        name: str
            CGD graph ID.

    **returns:**
        str
            CGD text.
    '''
    if isinstance(x, Atoms):
        topo = TopologyExtractor(ase_atoms=x)
        return topo.build_cgd(
            method=method,
            name=name,
        )

    if isinstance(x, (str, Path)):
        p = Path(x)
        if p.exists() and p.is_file():
            topo = TopologyExtractor(filename=str(p))
            return topo.build_cgd(
                method=method,
                name=name,
            )
        raise FileNotFoundError(f"Input path not found: {p}")

    return _make_cgd_from_input(
        _atoms_from_pymatgen(x),
        method=method,
        name=name,
    )


def identify_topology(
    x: Union[CGDLike, Atoms, "PymatgenStructure"],
    *,
    input_is_cgd: Optional[bool] = None,
    method: str = "all_node",
    name: str = "net",
    dedup_parallel_edges: bool = True,  # ignored
    dedup_mode: str = "shift",          # ignored
    java: str = "java",
    timeout_s: int = 30,
    xmx: str = "1024m",
    keep_tmp: bool = False,
) -> SystreResult:
    '''
    Identify the topology of a structure or CGD input using Systre.

    This is the main high-level API of the module. It accepts:
        1) a CGD file path,
        2) CGD text,
        3) an ASE Atoms object,
        4) a pymatgen Structure object,
        5) a structure file path such as CIF.

    If the input is not already a CGD file or CGD text, a CGD representation
    is generated first using `TopologyExtractor`.

    **parameters:**
        x:
            Input object. This can be a CGD path, CGD text, ASE atoms object,
            pymatgen Structure object, or a structure file path.

        input_is_cgd: bool or None
            If True, force the input to be treated as CGD.
            If False or None, the input is auto-detected.

        method: str
            Topology extraction mode used when a CGD must first be generated.

        name: str
            CGD graph ID used when generating a CGD representation.

        dedup_parallel_edges: bool
            Deprecated and ignored. Retained only for API compatibility.

        dedup_mode: str
            Deprecated and ignored. Retained only for API compatibility.

        java: str
            Java executable name or explicit path.

        timeout_s: int
            Timeout in seconds for the Systre run.

        xmx: str
            Maximum Java heap size.

        keep_tmp: bool
            If True, temporary CGD files are kept on disk.

    **returns:**
        SystreResult
            Object containing the identified topology, stdout, stderr,
            and optionally the CGD file path.
    '''
    _ = (dedup_parallel_edges, dedup_mode)
    tmp_path: Optional[str] = None

    try:
        if isinstance(x, (str, Path)):
            s = str(x)
            p = Path(s)

            if input_is_cgd is True:
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

            if p.exists() and p.is_file():
                cgd_text = _make_cgd_from_input(
                    p,
                    method=method,
                    name=name,
                )
                fd, tmp_path = tempfile.mkstemp(prefix="mofstructure_", suffix=".cgd")
                os.close(fd)
                Path(tmp_path).write_text(cgd_text, encoding="utf-8")
                proc = run_systre_on_cgd_path(tmp_path, java=java, timeout_s=timeout_s, xmx=xmx)
                topo = parse_systre_topology(proc.stdout)
                return SystreResult(topo, proc.stdout, proc.stderr, cgd_path=tmp_path)

            raise FileNotFoundError(f"Path not found: {p}")

        cgd_text = _make_cgd_from_input(
            x,
            method=method,
            name=name,
        )
        fd, tmp_path = tempfile.mkstemp(prefix="mofstructure_", suffix=".cgd")
        os.close(fd)
        Path(tmp_path).write_text(cgd_text, encoding="utf-8")
        proc = run_systre_on_cgd_path(tmp_path, java=java, timeout_s=timeout_s, xmx=xmx)
        topo = parse_systre_topology(proc.stdout)
        return SystreResult(topo, proc.stdout, proc.stderr, cgd_path=tmp_path)

    except subprocess.TimeoutExpired as exc:
        return SystreResult(
            "TIMEOUT",
            stdout=getattr(exc, "stdout", "") or "",
            stderr=getattr(exc, "stderr", "") or "",
            cgd_path=tmp_path,
        )

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
    '''
    Identify topologies for all supported files in a folder.

    The function scans a directory for CGD files or structure files, runs
    `identify_topology()` on each of them, and returns a dictionary mapping
    file paths to `SystreResult` objects.

    **parameters:**
        folder: str or Path
            Path to the folder containing input files.

        patterns: tuple
            Allowed filename suffixes to process.

        recursive: bool
            If True, recurse into subfolders.

        **kwargs:
            Additional keyword arguments passed directly to `identify_topology()`.

    **returns:**
        python dictionary
            Mapping:
                filepath -> SystreResult

    **raises:**
        NotADirectoryError:
            If the provided folder path is not a valid directory.
    '''
    folder = Path(folder)
    if not folder.exists() or not folder.is_dir():
        raise NotADirectoryError(f"Not a folder: {folder}")

    results: Dict[str, SystreResult] = {}
    allowed = {x.lower() for x in patterns}
    it = folder.rglob("*") if recursive else folder.glob("*")

    for p in it:
        if not p.is_file():
            continue
        if p.suffix.lower() not in allowed:
            continue
        try:
            res = identify_topology(p, **kwargs)
        except Exception as exc:
            results[str(p)] = SystreResult(
                "ERROR",
                stdout="",
                stderr=str(exc),
                cgd_path=str(p) if p.suffix.lower() == ".cgd" else None,
            )
        else:
            results[str(p)] = res

    return results