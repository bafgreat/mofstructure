#!/usr/bin/env python3
from __future__ import annotations
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"

import argparse
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

from ase.atoms import Atoms
from mofstructure.generate_cgd import TopologyExtractor

try:
    from importlib.resources import files as ir_files
except Exception:
    ir_files = None


logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )


InputLike = Union[str, Path, Atoms, "PymatgenStructure"]
Coord3D = Tuple[float, float, float]
Cell6 = Tuple[float, float, float, float, float, float]
Edge3D = Tuple[Coord3D, Coord3D]


@dataclass(frozen=True)
class RelaxedComponent:
    '''
    Container for one relaxed topology component parsed from Systre output.

    **parameters:**
        component_index: int
            Index of the component in the Systre output.

        dimension: int
            Periodicity of the component.

        cell: tuple
            Relaxed cell parameters in the form
            (a, b, c, alpha, beta, gamma).

        nodes: python dictionary
            Mapping:
                node_label -> (x, y, z)

        edges: list
            List of relaxed edges as coordinate pairs.

        edge_centers: list
            List of relaxed edge centres.

        coordination_number: python dictionary
            Mapping:
                node_label -> coordination number

        given_space_group: str, optional
            Given space group reported by Systre.

        ideal_space_group: str, optional
            Ideal space group reported by Systre.

        rcsr_name: str, optional
            RCSR name reported by Systre.

        td10: int, optional
            TD10 value reported by Systre.
    '''
    component_index: int
    dimension: int
    cell: Cell6
    nodes: Dict[str, Coord3D]
    edges: List[Edge3D]
    edge_centers: List[Coord3D]
    coordination_number: Dict[str, int]
    given_space_group: Optional[str] = None
    ideal_space_group: Optional[str] = None
    rcsr_name: Optional[str] = None
    td10: Optional[int] = None

    @property
    def is_named_net(self) -> bool:
        '''
        Check whether the relaxed component has an RCSR name.

        **returns:**
            bool
        '''
        return bool(self.rcsr_name)

    @property
    def n_nodes(self) -> int:
        '''
        Return the number of relaxed nodes.

        **returns:**
            int
        '''
        return len(self.nodes)

    @property
    def n_edges(self) -> int:
        '''
        Return the number of relaxed edges.

        **returns:**
            int
        '''
        return len(self.edges)

    def normalized_payload(self, decimals: int = 8) -> Dict[str, object]:
        '''
        Return a normalized payload that can be serialized and hashed.

        The payload contains only topology-relevant relaxed information.
        Floating-point values are rounded to make the hash more stable.

        **parameters:**
            decimals: int
                Number of decimal places used when rounding floats.

        **returns:**
            python dictionary
                JSON-serializable normalized topology payload.
        '''
        nodes = [
            {
                "label": str(label),
                "xyz": _round_coord(xyz, decimals),
                "coordination_number": int(self.coordination_number.get(label, 0)),
            }
            for label, xyz in sorted(self.nodes.items(), key=lambda kv: _sort_key(kv[0]))
        ]

        edges = [
            {
                "start": _round_coord(p, decimals),
                "end": _round_coord(q, decimals),
            }
            for p, q in _canonicalize_edges(self.edges, decimals=decimals)
        ]

        edge_centers = [
            _round_coord(center, decimals)
            for center in sorted(
                (_round_coord(center, decimals) for center in self.edge_centers),
                key=lambda xyz: (xyz[0], xyz[1], xyz[2]),
            )
        ]

        return {
            "dimension": int(self.dimension),
            "td10": None if self.td10 is None else int(self.td10),
            "ideal_space_group": self.ideal_space_group,
            "cell": _round_cell(self.cell, decimals),
            "nodes": nodes,
            "edges": edges,
            "edge_centers": edge_centers,
        }

    def topology_key(self, decimals: int = 8) -> str:
        '''
        Return a deterministic JSON string for hashing or indexing.

        **parameters:**
            decimals: int
                Number of decimal places used when rounding floats.

        **returns:**
            str
                Canonical JSON string.
        '''
        return _json_key(self.normalized_payload(decimals=decimals))

    def topology_hash(self, algorithm: str = "sha256", decimals: int = 8) -> str:
        '''
        Compute a stable topology hash from the normalized payload.

        **parameters:**
            algorithm: str
                Hash algorithm supported by hashlib.

            decimals: int
                Number of decimal places used when rounding floats.

        **returns:**
            str
                Hex digest.
        '''
        return _hash_text(self.topology_key(decimals=decimals), algorithm=algorithm)

    def crystal2_text(self, name: Optional[str] = None, include_edge_centers: bool = True) -> str:
        '''
        Return the CRYSTAL2 CGD representation as text.

        **parameters:**
            name: str, optional
                Name written to the CRYSTAL2 block.

            include_edge_centers: bool
                If True, write edge centres as commented metadata lines.

        **returns:**
            str
                CRYSTAL2 CGD text.
        '''
        return crystal2_text_from_component(
            self,
            name=name or self.rcsr_name or f"component_{self.component_index}",
            include_edge_centers=include_edge_centers,
        )

    def crystal2_hash(
        self,
        algorithm: str = "sha256",
        name: Optional[str] = None,
        include_edge_centers: bool = True,
    ) -> str:
        '''
        Compute a hash directly from the CRYSTAL2 text.

        **parameters:**
            algorithm: str
                Hash algorithm supported by hashlib.

            name: str, optional
                Name written to the CRYSTAL2 block before hashing.

            include_edge_centers: bool
                If True, include commented edge-centre metadata.

        **returns:**
            str
                Hex digest.
        '''
        return _hash_text(
            self.crystal2_text(name=name, include_edge_centers=include_edge_centers),
            algorithm=algorithm,
        )


@dataclass(frozen=True)
class SystreResult:
    '''
    Container class for a Systre topology-identification result.

    **parameters:**
        topology: str
            Identified topology name or sentinel value such as
            "UNKNOWN", "ERROR", or "TIMEOUT".

        stdout: str
            Standard output from Systre.

        stderr: str
            Standard error from Systre.

        cgd_path: str, optional
            Path to the CGD file passed to Systre.
    '''
    topology: str
    stdout: str
    stderr: str
    cgd_path: Optional[str] = None


class SystreTopology:
    '''
    High-level class for running Systre, parsing relaxed output, and
    generating stable topology hashes.

    **parameters:**
        input_object:
            Structure-like input. Supported values are:
                - ASE atoms object
                - structure filename
                - CGD filepath
                - CGD text
                - pymatgen Structure object

        input_is_cgd: bool, optional
            If True, force the input to be treated as CGD.

        method: str
            Topology extraction mode passed to `TopologyExtractor.build_cgd()`
            when a CGD must first be generated.

        name: str
            CGD graph ID used when generating the initial CGD.

        java: str, optional
            Java executable name or explicit path. If None, jdk4py is tried first,
            then the system PATH.

        timeout_s: int
            Timeout for the Systre call in seconds.

        xmx: str
            Maximum Java heap size.

        keep_tmp: bool
            If True, temporary CGD files are kept on disk.
    '''
    def __init__(
        self,
        input_object: InputLike,
        *,
        input_is_cgd: Optional[bool] = None,
        method: str = "all_node",
        name: str = "net",
        java: Optional[str] = None,
        timeout_s: int = 30,
        xmx: str = "1024m",
        keep_tmp: bool = False,
    ):
        self.input_object = input_object
        self.input_is_cgd = input_is_cgd
        self.method = method
        self.name = name
        self.java = java
        self.timeout_s = timeout_s
        self.xmx = xmx
        self.keep_tmp = keep_tmp
        self._cached_result: Optional[SystreResult] = None
        self._cached_components: Optional[List[RelaxedComponent]] = None

    def identify(self, force: bool = False) -> SystreResult:
        '''
        Run Systre and return the topology-identification result.

        Results are cached after the first call unless `force=True`.

        **parameters:**
            force: bool
                If True, rerun Systre even if a cached result exists.

        **returns:**
            SystreResult
        '''
        if self._cached_result is None or force:
            self._cached_result = identify_topology(
                self.input_object,
                input_is_cgd=self.input_is_cgd,
                method=self.method,
                name=self.name,
                java=self.java,
                timeout_s=self.timeout_s,
                xmx=self.xmx,
                keep_tmp=self.keep_tmp,
            )
            self._cached_components = None
        return self._cached_result

    def extract_relaxed_components(
        self,
        stdout: Optional[str] = None,
        force: bool = False,
    ) -> List[RelaxedComponent]:
        '''
        Parse relaxed components from Systre stdout.

        **parameters:**
            stdout: str, optional
                Systre stdout. If not provided, Systre is executed first.

            force: bool
                If True, rerun Systre when `stdout` is not provided.

        **returns:**
            list
                List of parsed relaxed components.
        '''
        if stdout is not None:
            return extract_relaxed_components(stdout)

        if self._cached_components is None or force:
            self._cached_components = extract_relaxed_components(self.identify(force=force).stdout)
        return self._cached_components

    def best_component(self, stdout: Optional[str] = None) -> Optional[RelaxedComponent]:
        '''
        Return the most informative relaxed component.

        Current strategy:
            rank by number of edges, then by number of nodes.

        **parameters:**
            stdout: str, optional
                Systre stdout. If not provided, Systre is executed first.

        **returns:**
            RelaxedComponent or None
        '''
        return _best_component(self.extract_relaxed_components(stdout=stdout))

    def named_components(self, stdout: Optional[str] = None) -> List[RelaxedComponent]:
        '''
        Return all relaxed components that have an RCSR name.

        **parameters:**
            stdout: str, optional
                Systre stdout. If not provided, Systre is executed first.

        **returns:**
            list
                List of named relaxed components.
        '''
        return [comp for comp in self.extract_relaxed_components(stdout=stdout) if comp.rcsr_name]

    def unnamed_components(self, stdout: Optional[str] = None) -> List[RelaxedComponent]:
        '''
        Return all relaxed components without an RCSR name.

        **parameters:**
            stdout: str, optional
                Systre stdout. If not provided, Systre is executed first.

        **returns:**
            list
                List of unnamed relaxed components.
        '''
        return [comp for comp in self.extract_relaxed_components(stdout=stdout) if not comp.rcsr_name]

    def topology_summary(self, decimals: int = 8, component_index: Optional[int] = None) -> Optional[Dict[str, object]]:
        '''
        Return a normalized topology summary for hashing or storage.

        **parameters:**
            decimals: int
                Number of decimal places used when rounding floats.

            component_index: int, optional
                Component index to use. If not provided, the best component is used.

        **returns:**
            python dictionary or None
        '''
        comp = self._select_component(component_index=component_index)
        return None if comp is None else comp.normalized_payload(decimals=decimals)

    def topology_hash(
        self,
        decimals: int = 8,
        component_index: Optional[int] = None,
        algorithm: str = "sha256",
    ) -> Optional[str]:
        '''
        Return the relaxed-topology hash for one parsed component.

        **parameters:**
            decimals: int
                Number of decimal places used when rounding floats.

            component_index: int, optional
                Component index to use. If not provided, the best component is used.

            algorithm: str
                Hash algorithm supported by hashlib.

        **returns:**
            str or None
                Hex digest, or None if no relaxed component was parsed.
        '''
        comp = self._select_component(component_index=component_index)
        return None if comp is None else comp.topology_hash(algorithm=algorithm, decimals=decimals)

    def hashes_for_unnamed_components(
        self,
        decimals: int = 8,
        algorithm: str = "sha256",
    ) -> Dict[int, str]:
        '''
        Return hashes for components that do not have an RCSR name.

        **parameters:**
            decimals: int
                Number of decimal places used when rounding floats.

            algorithm: str
                Hash algorithm supported by hashlib.

        **returns:**
            python dictionary
                Mapping:
                    component_index -> topology hash
        '''
        return {
            comp.component_index: comp.topology_hash(algorithm=algorithm, decimals=decimals)
            for comp in self.unnamed_components()
        }

    def crystal2_text(
        self,
        *,
        component_index: Optional[int] = None,
        include_edge_centers: bool = True,
        fallback_to_input_cgd: bool = False,
    ) -> Optional[str]:
        '''
        Return CRYSTAL2 CGD text for one parsed relaxed component.

        **parameters:**
            component_index: int, optional
                Component index to export. If not provided, the best component is used.

            include_edge_centers: bool
                If True, include commented edge-centre metadata.

            fallback_to_input_cgd: bool
                If True, convert the original PERIODIC_GRAPH CGD to a CRYSTAL2 wrapper
                when relaxed data is unavailable.

        **returns:**
            str or None
        '''
        comp = self._select_component(component_index=component_index)
        if comp is not None:
            return comp.crystal2_text(name=comp.rcsr_name or self.name, include_edge_centers=include_edge_centers)

        if not fallback_to_input_cgd:
            return None

        cgd_text = _input_to_cgd_text(
            self.input_object,
            input_is_cgd=self.input_is_cgd,
            method=self.method,
            name=self.name,
        )
        if cgd_text is None:
            return None
        return crystal2_fallback_from_periodic_graph(cgd_text, name=self.name)

    def crystal2_hash(
        self,
        *,
        component_index: Optional[int] = None,
        include_edge_centers: bool = True,
        fallback_to_input_cgd: bool = False,
        algorithm: str = "sha256",
    ) -> Optional[str]:
        '''
        Return a hash computed from the CRYSTAL2 CGD text.

        **parameters:**
            component_index: int, optional
                Component index to export. If not provided, the best component is used.

            include_edge_centers: bool
                If True, include commented edge-centre metadata.

            fallback_to_input_cgd: bool
                If True, convert the original PERIODIC_GRAPH CGD to a CRYSTAL2 wrapper
                when relaxed data is unavailable.

            algorithm: str
                Hash algorithm supported by hashlib.

        **returns:**
            str or None
        '''
        text = self.crystal2_text(
            component_index=component_index,
            include_edge_centers=include_edge_centers,
            fallback_to_input_cgd=fallback_to_input_cgd,
        )
        return None if text is None else _hash_text(text, algorithm=algorithm)

    def write_crystal2(
        self,
        out_path: Union[str, Path],
        *,
        component_index: Optional[int] = None,
        include_edge_centers: bool = True,
        fallback_to_input_cgd: bool = False,
    ) -> Path:
        '''
        Write CRYSTAL2 CGD text to disk.

        **parameters:**
            out_path: str or Path
                Output file path.

            component_index: int, optional
                Component index to export. If not provided, the best component is used.

            include_edge_centers: bool
                If True, include commented edge-centre metadata.

            fallback_to_input_cgd: bool
                If True, convert the original PERIODIC_GRAPH CGD to a CRYSTAL2 wrapper
                when relaxed data is unavailable.

        **returns:**
            pathlib.Path
                Written file path.

        **raises:**
            ValueError:
                If no CRYSTAL2 text can be generated.
        '''
        text = self.crystal2_text(
            component_index=component_index,
            include_edge_centers=include_edge_centers,
            fallback_to_input_cgd=fallback_to_input_cgd,
        )
        if text is None:
            raise ValueError("No CRYSTAL2 text could be generated from this input.")

        out_path = Path(out_path)
        out_path.write_text(text, encoding="utf-8")
        return out_path

    def _select_component(self, component_index: Optional[int] = None) -> Optional[RelaxedComponent]:
        '''
        Select one relaxed component.

        **parameters:**
            component_index: int, optional
                Explicit component index.

        **returns:**
            RelaxedComponent or None
        '''
        components = self.extract_relaxed_components()
        if component_index is None:
            return _best_component(components)
        for comp in components:
            if comp.component_index == component_index:
                return comp
        return None


def _round_coord(xyz: Coord3D, decimals: int) -> Coord3D:
    '''
    Round a 3D coordinate tuple.

    **parameters:**
        xyz: tuple
            Coordinate tuple.

        decimals: int
            Number of decimal places.

    **returns:**
        tuple
            Rounded coordinate tuple.
    '''
    return tuple(round(float(value), decimals) for value in xyz)


def _round_cell(cell: Cell6, decimals: int) -> Cell6:
    '''
    Round a six-parameter cell tuple.

    **parameters:**
        cell: tuple
            Cell tuple.

        decimals: int
            Number of decimal places.

    **returns:**
        tuple
            Rounded cell tuple.
    '''
    return tuple(round(float(value), decimals) for value in cell)


def _json_key(payload: Dict[str, object]) -> str:
    '''
    Serialize a payload into canonical JSON text.

    **parameters:**
        payload: python dictionary
            JSON-serializable payload.

    **returns:**
        str
            Canonical JSON string.
    '''
    return json.dumps(payload, sort_keys=True, separators=(",", ":"))


def _hash_text(text: str, algorithm: str = "sha256") -> str:
    '''
    Hash a text string.

    **parameters:**
        text: str
            Input text.

        algorithm: str
            Hash algorithm supported by hashlib.

    **returns:**
        str
            Hex digest.
    '''
    h = hashlib.new(algorithm)
    h.update(text.encode("utf-8"))
    return h.hexdigest()


def _best_component(components: Sequence[RelaxedComponent]) -> Optional[RelaxedComponent]:
    '''
    Select the most informative relaxed component.

    Ranking:
        1) number of edges
        2) number of nodes

    **parameters:**
        components: sequence
            Parsed relaxed components.

    **returns:**
        RelaxedComponent or None
    '''
    if not components:
        return None
    return sorted(components, key=lambda comp: (comp.n_edges, comp.n_nodes), reverse=True)[0]


def _is_cgd_text(text: str) -> bool:
    '''
    Check whether a string looks like CGD text.

    **parameters:**
        text: str
            Input string.

        **returns:**
            bool
    '''
    up = text.upper()
    return ("PERIODIC_GRAPH" in up or "CRYSTAL" in up) and ("EDGES" in up or "EDGE" in up)


def _java_install_hints() -> str:
    '''
    Return platform-specific Java installation hints.

    **returns:**
        str
    '''
    return (
        "Java runtime not found.\n\n"
        "Install one of the following:\n"
        "  - pip install jdk4py\n"
        "  - or install a system JRE/JDK and ensure 'java' is on PATH\n"
        "  - Adoptium (Temurin): https://adoptium.net/\n"
    )


def _resolve_jdk4py_java() -> Optional[str]:
    '''
    Try to resolve the Java executable provided by jdk4py.

    This helper supports several access patterns because package APIs may vary.

    **returns:**
        str or None
            Java executable path if found.
    '''
    try:
        import jdk4py  # type: ignore
    except Exception:
        return None

    candidates = []

    for attr in ("JAVA", "java", "JAVA_EXECUTABLE", "java_executable"):
        value = getattr(jdk4py, attr, None)
        if value:
            candidates.append(value)

    for attr in ("JAVA_HOME", "java_home"):
        value = getattr(jdk4py, attr, None)
        if value:
            home = Path(str(value))
            candidates.append(home / "bin" / ("java.exe" if sys.platform.startswith("win") else "java"))

    for value in candidates:
        path = Path(str(value))
        if path.exists():
            return str(path)

    return None


def find_java(java: Optional[str] = None) -> str:
    '''
    Find a working Java executable.

    Resolution order:
        1) explicit `java` argument
        2) jdk4py bundled Java
        3) system PATH lookup

    **parameters:**
        java: str, optional
            Java executable name or explicit path.

    **returns:**
        str
            Working Java command.

    **raises:**
        RuntimeError:
            If Java cannot be located.
    '''
    if java:
        exe = shutil.which(java) if os.path.sep not in java else (java if os.path.exists(java) else None)
        if exe:
            return exe

    jdk4py_java = _resolve_jdk4py_java()
    if jdk4py_java:
        return jdk4py_java

    exe = shutil.which("java")
    if exe:
        return exe

    raise RuntimeError(_java_install_hints())


def _resource_path(rel: str) -> str:
    '''
    Return the absolute filesystem path of a packaged resource.

    **parameters:**
        rel: str
            Relative resource path inside the package.

    **returns:**
        str
            Absolute path to the resource.
    '''
    direct = Path(__file__).resolve().parent / rel
    if direct.exists():
        return str(direct)

    if ir_files is None:
        raise RuntimeError(f"Cannot locate resource {rel}. importlib.resources not available.")

    pkg_root = ir_files("mofstructure")
    candidate = pkg_root.joinpath(rel)

    try:
        candidate_path = Path(str(candidate))
        if candidate_path.exists():
            return str(candidate_path)
    except Exception:
        pass

    data = candidate.read_bytes()
    suffix = Path(rel).suffix
    fd, tmp = tempfile.mkstemp(prefix="mofstructure_resource_", suffix=suffix)
    os.close(fd)
    with open(tmp, "wb") as handle:
        handle.write(data)
    return tmp


def systre_command(
    *,
    java: Optional[str] = None,
    xmx: str = "1024m",
    jar_relpath: str = "bin/Systre-19.6.0.jar",
    rcsr_relpath: str = "db/RCSRnets-2019-06-01.arc",
) -> List[str]:
    '''
    Build the Systre command-line call.

    **parameters:**
        java: str, optional
            Java executable name or explicit path.

        xmx: str
            Maximum Java heap size.

        jar_relpath: str
            Relative path to the packaged Systre JAR file.

        rcsr_relpath: str
            Relative path to the packaged RCSR archive.

    **returns:**
        list
            Command list suitable for `subprocess.run`.
    '''
    return [
        find_java(java),
        f"-Xmx{xmx}",
        "-cp",
        _resource_path(jar_relpath),
        "org.gavrog.apps.systre.SystreCmdline",
        _resource_path(rcsr_relpath),
    ]


def run_systre_on_cgd_path(
    cgd_path: Union[str, Path],
    *,
    java: Optional[str] = None,
    timeout_s: int = 30,
    xmx: str = "1024m",
) -> subprocess.CompletedProcess:
    '''
    Run Systre on a CGD file stored on disk.

    **parameters:**
        cgd_path: str or Path
            Path to a CGD file.

        java: str, optional
            Java executable name or explicit path.

        timeout_s: int
            Timeout in seconds.

        xmx: str
            Maximum Java heap size.

    **returns:**
        subprocess.CompletedProcess
    '''
    cmd = systre_command(java=java, xmx=xmx) + [str(cgd_path)]
    return subprocess.run(
        cmd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=timeout_s,
    )


def _completed_process_to_result(
    proc: subprocess.CompletedProcess,
    *,
    cgd_path: Optional[str] = None,
) -> SystreResult:
    '''
    Convert a completed Systre subprocess into a SystreResult.

    **parameters:**
        proc: subprocess.CompletedProcess
            Completed subprocess result.

        cgd_path: str, optional
            Path to the CGD file used.

    **returns:**
        SystreResult
    '''
    return SystreResult(
        topology=parse_systre_topology(proc.stdout),
        stdout=proc.stdout,
        stderr=proc.stderr,
        cgd_path=cgd_path,
    )


def _run_systre_with_cgd_text(
    cgd_text: str,
    *,
    java: Optional[str],
    timeout_s: int,
    xmx: str,
    keep_tmp: bool,
) -> SystreResult:
    '''
    Run Systre from CGD text.

    **parameters:**
        cgd_text: str
            CGD text.

        java: str, optional
            Java executable name or explicit path.

        timeout_s: int
            Timeout in seconds.

        xmx: str
            Maximum Java heap size.

        keep_tmp: bool
            If True, keep the temporary CGD file.

    **returns:**
        SystreResult
    '''
    tmp_path = None
    try:
        fd, tmp_path = tempfile.mkstemp(prefix="mofstructure_", suffix=".cgd")
        os.close(fd)
        Path(tmp_path).write_text(cgd_text, encoding="utf-8")
        proc = run_systre_on_cgd_path(tmp_path, java=java, timeout_s=timeout_s, xmx=xmx)
        return _completed_process_to_result(proc, cgd_path=tmp_path)
    except subprocess.TimeoutExpired as exc:
        return SystreResult(
            "TIMEOUT",
            stdout=getattr(exc, "stdout", "") or "",
            stderr=getattr(exc, "stderr", "") or "",
            cgd_path=tmp_path,
        )
    finally:
        if tmp_path and not keep_tmp:
            try:
                os.remove(tmp_path)
            except Exception:
                pass


def _input_to_cgd_text(
    x: InputLike,
    *,
    input_is_cgd: Optional[bool] = None,
    method: str = "all_node",
    name: str = "net",
) -> Optional[str]:
    '''
    Convert supported input into CGD text.

    **parameters:**
        x:
            Supported input object.

        input_is_cgd: bool, optional
            If True, force CGD interpretation.

        method: str
            TopologyExtractor method.

        name: str
            Graph name.

        **returns:**
            str or None
    '''
    if isinstance(x, (str, Path)):
        s = str(x)
        p = Path(s)

        if input_is_cgd is True:
            if p.exists() and p.is_file():
                return p.read_text(encoding="utf-8")
            if _is_cgd_text(s):
                return s
            raise ValueError("input_is_cgd=True but input is neither a CGD file nor CGD text.")

        if p.exists() and p.is_file() and p.suffix.lower() == ".cgd":
            return p.read_text(encoding="utf-8")

        if _is_cgd_text(s):
            return s

    try:
        return _make_cgd_from_input(x, method=method, name=name)
    except Exception:
        return None


def _split_component_blocks(stdout: str) -> List[Tuple[int, str]]:
    '''
    Split Systre stdout into per-component blocks.

    **parameters:**
        stdout: str
            Systre standard output.

    **returns:**
        list
            List of tuples:
                (component_index, block_text)
    '''
    if "Processing component " not in stdout:
        return [(1, stdout)]

    parts = re.split(r"\n\s*Processing component\s+(\d+):\s*\n", stdout)
    if len(parts) < 3:
        return [(1, stdout)]

    pre = parts[0]
    out = []
    for i in range(1, len(parts) - 1, 2):
        idx = int(parts[i])
        block = parts[i + 1]
        out.append((idx, pre + "\nProcessing component %d:\n" % idx + block))
    return out


def _parse_td10(block: str) -> Optional[int]:
    '''
    Parse the TD10 value from a Systre block.

    **parameters:**
        block: str
            One Systre component block.

    **returns:**
        int or None
    '''
    match = re.search(r"(?m)^\s*TD10\s*=\s*(\d+)\s*$", block)
    return int(match.group(1)) if match else None


def _parse_rcsr_name(block: str) -> Optional[str]:
    '''
    Parse the RCSR name from a Systre block.

    **parameters:**
        block: str
            One Systre component block.

    **returns:**
        str or None
    '''
    match = re.search(
        r"Structure\s+was\s+identified\s+with\s+RCSR\s+symbol:\s*(?:\r?\n)+\s*Name:\s*([A-Za-z0-9_\-]+)",
        block,
        flags=re.IGNORECASE,
    )
    if match:
        return match.group(1)

    match = re.search(
        r"Structure\s+was\s+found\s+in\s+archive.*?(?:\r?\n)+\s*Name:\s*([A-Za-z0-9_\-]+)",
        block,
        flags=re.IGNORECASE | re.DOTALL,
    )
    return match.group(1) if match else None


def _parse_dimension(block: str) -> Optional[int]:
    '''
    Parse the periodic dimension from a Systre block.

    **parameters:**
        block: str
            One Systre component block.

    **returns:**
        int or None
    '''
    for pattern in (
        r"described\s+as\s+(\d+)\s*-\s*periodic",
        r"described\s+as\s+(\d+)\s+periodic",
        r"dimension\s*=\s*(\d+)",
    ):
        match = re.search(pattern, block, flags=re.IGNORECASE)
        if match:
            return int(match.group(1))
    return None


def _parse_relaxed_cell(block: str) -> Optional[Cell6]:
    '''
    Parse relaxed cell parameters from a Systre block.

    For 2D nets, the returned cell is normalized to:
        (a, b, 1.0, 90.0, 90.0, gamma)

    **parameters:**
        block: str
            One Systre component block.

    **returns:**
        tuple or None
            (a, b, c, alpha, beta, gamma)
    '''
    match = re.search(
        r"Relaxed cell parameters:\s*(?:\r?\n)\s*"
        r"a\s*=\s*([0-9.\-Ee+]+)\s*,\s*b\s*=\s*([0-9.\-Ee+]+)\s*,\s*c\s*=\s*([0-9.\-Ee+]+)\s*(?:\r?\n)\s*"
        r"alpha\s*=\s*([0-9.\-Ee+]+)\s*,\s*beta\s*=\s*([0-9.\-Ee+]+)\s*,\s*gamma\s*=\s*([0-9.\-Ee+]+)",
        block,
        flags=re.IGNORECASE,
    )
    if match:
        return tuple(float(match.group(i)) for i in range(1, 7))

    match = re.search(
        r"Relaxed cell parameters:\s*(?:\r?\n)\s*"
        r"a\s*=\s*([0-9.\-Ee+]+)\s*,\s*b\s*=\s*([0-9.\-Ee+]+)\s*,\s*gamma\s*=\s*([0-9.\-Ee+]+)",
        block,
        flags=re.IGNORECASE,
    )
    if match:
        return (
            float(match.group(1)),
            float(match.group(2)),
            1.0,
            90.0,
            90.0,
            float(match.group(3)),
        )

    match = re.search(
        r"Relaxed cell parameters:\s*(?:\r?\n|\s+)"
        r"a\s*=\s*([0-9.\-Ee+]+)\s*,\s*b\s*=\s*([0-9.\-Ee+]+)\s*,\s*c\s*=\s*([0-9.\-Ee+]+)\s*"
        r"(?:,|\s+)\s*alpha\s*=\s*([0-9.\-Ee+]+)\s*,\s*beta\s*=\s*([0-9.\-Ee+]+)\s*,\s*gamma\s*=\s*([0-9.\-Ee+]+)",
        block,
        flags=re.IGNORECASE,
    )
    if match:
        return tuple(float(match.group(i)) for i in range(1, 7))

    return None


def _parse_space_groups(block: str) -> Tuple[Optional[str], Optional[str]]:
    '''
    Parse the given and ideal space groups from a Systre block.

    **parameters:**
        block: str
            One Systre component block.

    **returns:**
        tuple
            (given_space_group, ideal_space_group)
    '''
    match_given = re.search(
        r"Given\s+space\s+group\s+is\s+([A-Za-z0-9\-_/]+)\s*\.",
        block,
        flags=re.IGNORECASE,
    )
    match_ideal = re.search(
        r"Ideal\s+space\s+group\s+is\s+([A-Za-z0-9\-_/]+)\s*\.",
        block,
        flags=re.IGNORECASE,
    )
    return (
        match_given.group(1) if match_given else None,
        match_ideal.group(1) if match_ideal else None,
    )


def _parse_edge_centers(block: str) -> List[Coord3D]:
    '''
    Parse edge centres from a Systre block.

    **parameters:**
        block: str
            One Systre component block.

    **returns:**
        list
            List of edge-centre coordinates.
    '''
    match = re.search(
        r"Edge\s+centers:\s*(.*?)(?:\r?\n)\s*(?:Edge statistics:|Angle statistics:|Shortest non-bonded distance|Degrees of freedom:|Finished structure|\Z)",
        block,
        flags=re.DOTALL | re.IGNORECASE,
    )
    if not match:
        return []

    centers = []
    for line in match.group(1).strip().splitlines():
        line = line.strip()
        if not line:
            continue

        match3 = re.match(r"([0-9.\-Ee+]+)\s+([0-9.\-Ee+]+)\s+([0-9.\-Ee+]+)", line)
        if match3:
            centers.append((float(match3.group(1)), float(match3.group(2)), float(match3.group(3))))
            continue

        match2 = re.match(r"([0-9.\-Ee+]+)\s+([0-9.\-Ee+]+)", line)
        if match2:
            centers.append((float(match2.group(1)), float(match2.group(2)), 0.0))

    return centers


def _parse_relaxed_positions(block: str) -> Optional[Dict[str, Coord3D]]:
    '''
    Parse relaxed node positions from a Systre block.

    **parameters:**
        block: str
            One Systre component block.

    **returns:**
        python dictionary or None
            Mapping:
                node_label -> (x, y, z)
    '''
    match = re.search(r"Relaxed positions:\s*(.*?)(?:\r?\n)\s*Edges:", block, flags=re.DOTALL | re.IGNORECASE)
    if not match:
        return None

    nodes = {}
    for line in match.group(1).strip().splitlines():
        line = line.strip()

        match3 = re.match(r"Node\s+(\S+):\s+([0-9.\-Ee+]+)\s+([0-9.\-Ee+]+)\s+([0-9.\-Ee+]+)", line)
        if match3:
            nodes[match3.group(1)] = (
                float(match3.group(2)),
                float(match3.group(3)),
                float(match3.group(4)),
            )
            continue

        match2 = re.match(r"Node\s+(\S+):\s+([0-9.\-Ee+]+)\s+([0-9.\-Ee+]+)", line)
        if match2:
            nodes[match2.group(1)] = (
                float(match2.group(2)),
                float(match2.group(3)),
                0.0,
            )

    return nodes


def _parse_coordination_numbers(block: str) -> Dict[str, int]:
    '''
    Parse node coordination numbers from the coordination-sequence block.

    **parameters:**
        block: str
            One Systre component block.

    **returns:**
        python dictionary
            Mapping:
                node_label -> coordination number
    '''
    match = re.search(
        r"Coordination\s+sequences:\s*(.*?)(?:\r?\n\s*\r?\n|\r?\n\s*TD10\s*=|\Z)",
        block,
        flags=re.DOTALL | re.IGNORECASE,
    )
    if not match:
        return {}

    out = {}
    for line in match.group(1).splitlines():
        line = line.strip()
        row = re.match(r"Node\s+(\S+):\s+(\d+)\b", line)
        if row:
            out[row.group(1)] = int(row.group(2))
    return out


def _parse_edges(block: str) -> Optional[List[Edge3D]]:
    '''
    Parse relaxed edges from a Systre block.

    **parameters:**
        block: str
            One Systre component block.

    **returns:**
        list or None
            List of relaxed edges.
    '''
    match = re.search(
        r"Edges:\s*(.*?)(?:\r?\n)\s*(?:Edge centers:|Edge statistics:|Angle statistics:|Shortest non-bonded distance|Degrees of freedom:|Finished structure|\Z)",
        block,
        flags=re.DOTALL | re.IGNORECASE,
    )
    if not match:
        return None

    edges = []
    for line in match.group(1).strip().splitlines():
        line = line.strip()

        match3 = re.match(
            r"([0-9.\-Ee+]+)\s+([0-9.\-Ee+]+)\s+([0-9.\-Ee+]+)\s+<->\s+([0-9.\-Ee+]+)\s+([0-9.\-Ee+]+)\s+([0-9.\-Ee+]+)",
            line,
        )
        if match3:
            edges.append(
                (
                    (float(match3.group(1)), float(match3.group(2)), float(match3.group(3))),
                    (float(match3.group(4)), float(match3.group(5)), float(match3.group(6))),
                )
            )
            continue

        match2 = re.match(
            r"([0-9.\-Ee+]+)\s+([0-9.\-Ee+]+)\s+<->\s+([0-9.\-Ee+]+)\s+([0-9.\-Ee+]+)",
            line,
        )
        if match2:
            edges.append(
                (
                    (float(match2.group(1)), float(match2.group(2)), 0.0),
                    (float(match2.group(3)), float(match2.group(4)), 0.0),
                )
            )

    return edges


def parse_relaxed_component(block: str, component_index: int) -> Optional[RelaxedComponent]:
    '''
    Parse one relaxed component from a Systre block.

    **parameters:**
        block: str
            One Systre component block.

        component_index: int
            Component index.

    **returns:**
        RelaxedComponent or None
    '''
    cell = _parse_relaxed_cell(block)
    nodes = _parse_relaxed_positions(block)
    edges = _parse_edges(block)
    if not cell or not nodes or edges is None:
        return None

    given_sg, ideal_sg = _parse_space_groups(block)
    return RelaxedComponent(
        component_index=component_index,
        dimension=_parse_dimension(block) or 0,
        cell=cell,
        nodes=nodes,
        edges=edges,
        edge_centers=_parse_edge_centers(block),
        coordination_number=_parse_coordination_numbers(block),
        given_space_group=given_sg,
        ideal_space_group=ideal_sg,
        rcsr_name=_parse_rcsr_name(block),
        td10=_parse_td10(block),
    )


def extract_relaxed_components(stdout: str) -> List[RelaxedComponent]:
    '''
    Parse all relaxed components from Systre stdout.

    **parameters:**
        stdout: str
            Systre standard output.

    **returns:**
        list
            List of parsed relaxed components.
    '''
    components = []
    for idx, block in _split_component_blocks(stdout):
        comp = parse_relaxed_component(block, component_index=idx)
        if comp is not None:
            components.append(comp)
    return components


def parse_systre_topology(stdout: str) -> str:
    '''
    Parse the identified topology name from Systre stdout.

    If several components are present and they do not all map to the same
    topology, the function returns "MISMATCH".

    **parameters:**
        stdout: str
            Systre standard output.

    **returns:**
        str
            Identified topology name or one of the sentinel values:
                - "UNKNOWN"
                - "ERROR"
                - "MISMATCH"
    '''
    topologies = []
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
                    topologies.append(topologies[int(comps[-1]) - 1])
                except Exception:
                    return "ERROR"
            else:
                return "ERROR"
            continue

        if "ERROR" in line:
            return "ERROR"
        if "Structure was found in archive" in line:
            topology_line = True
        elif line == "Structure is new for this run.":
            topologies.append("UNKNOWN")
        elif line == "Structure already seen in this run.":
            repeat_line = True

    if not topologies:
        return "ERROR"

    first = topologies[0]
    return first if all(topo == first for topo in topologies[1:]) else "MISMATCH"


def _atoms_from_pymatgen(obj) -> Atoms:
    '''
    Convert a pymatgen Structure object to an ASE atoms object.

    **parameters:**
        obj:
            pymatgen Structure object.

    **returns:**
        ASE atoms object

    **raises:**
        ImportError:
            If pymatgen is not installed.
    '''
    try:
        from pymatgen.io.ase import AseAtomsAdaptor
    except Exception as exc:
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

    **parameters:**
        x:
            ASE atoms object, pymatgen Structure object, or structure filepath.

        method: str
            Topology mode passed to `TopologyExtractor.build_cgd()`.

        name: str
            CGD graph ID.

    **returns:**
        str
            CGD text.
    '''
    if isinstance(x, Atoms):
        return TopologyExtractor(ase_atoms=x).build_cgd(method=method, name=name)

    if isinstance(x, (str, Path)):
        path = Path(x)
        if path.exists() and path.is_file():
            return TopologyExtractor(filename=str(path)).build_cgd(method=method, name=name)
        raise FileNotFoundError(f"Input path not found: {path}")

    return _make_cgd_from_input(_atoms_from_pymatgen(x), method=method, name=name)


def identify_topology(
    x: InputLike,
    *,
    input_is_cgd: Optional[bool] = None,
    method: str = "all_node",
    name: str = "net",
    java: Optional[str] = None,
    timeout_s: int = 30,
    xmx: str = "1024m",
    keep_tmp: bool = False,
) -> SystreResult:
    '''
    Identify the topology of a structure or CGD input using Systre.

    **parameters:**
        x:
            Input object. Supported values are:
                - CGD filepath
                - CGD text
                - ASE atoms object
                - pymatgen Structure object
                - structure filepath readable by ASE

        input_is_cgd: bool, optional
            If True, force the input to be treated as CGD.

        method: str
            Topology mode used when a CGD must first be generated.

        name: str
            CGD graph ID.

        java: str, optional
            Java executable name or explicit path.

        timeout_s: int
            Timeout in seconds.

        xmx: str
            Maximum Java heap size.

        keep_tmp: bool
            If True, temporary CGD files are kept.

    **returns:**
        SystreResult
    '''
    if isinstance(x, (str, Path)):
        s = str(x)
        p = Path(s)

        if input_is_cgd is not True and p.exists() and p.is_file() and p.suffix.lower() == ".cgd":
            try:
                proc = run_systre_on_cgd_path(p, java=java, timeout_s=timeout_s, xmx=xmx)
                return _completed_process_to_result(proc, cgd_path=str(p))
            except subprocess.TimeoutExpired as exc:
                return SystreResult(
                    "TIMEOUT",
                    stdout=getattr(exc, "stdout", "") or "",
                    stderr=getattr(exc, "stderr", "") or "",
                    cgd_path=str(p),
                )

        if input_is_cgd is True or _is_cgd_text(s):
            cgd_text = _input_to_cgd_text(
                x,
                input_is_cgd=input_is_cgd,
                method=method,
                name=name,
            )
            if cgd_text is None:
                raise ValueError("Could not interpret the provided input as CGD text.")
            return _run_systre_with_cgd_text(
                cgd_text,
                java=java,
                timeout_s=timeout_s,
                xmx=xmx,
                keep_tmp=keep_tmp,
            )

    cgd_text = _make_cgd_from_input(x, method=method, name=name)
    return _run_systre_with_cgd_text(
        cgd_text,
        java=java,
        timeout_s=timeout_s,
        xmx=xmx,
        keep_tmp=keep_tmp,
    )


def identify_topology_batch(
    folder: Union[str, Path],
    *,
    patterns: Tuple[str, ...] = (".cgd", ".cif", ".vasp", ".poscar"),
    recursive: bool = True,
    **kwargs,
) -> Dict[str, SystreResult]:
    '''
    Identify topologies for all supported files in a folder.

    **parameters:**
        folder: str or Path
            Input folder.

        patterns: tuple
            Allowed filename suffixes.

        recursive: bool
            If True, recurse into subfolders.

        **kwargs:
            Additional keyword arguments passed to `identify_topology()`.

    **returns:**
        python dictionary
            Mapping:
                filepath -> SystreResult
    '''
    folder = Path(folder)
    if not folder.exists() or not folder.is_dir():
        raise NotADirectoryError(f"Not a folder: {folder}")

    results = {}
    allowed = {suffix.lower() for suffix in patterns}
    iterator = folder.rglob("*") if recursive else folder.glob("*")

    for path in iterator:
        if not path.is_file() or path.suffix.lower() not in allowed:
            continue
        try:
            results[str(path)] = identify_topology(path, **kwargs)
        except Exception as exc:
            results[str(path)] = SystreResult(
                "ERROR",
                stdout="",
                stderr=str(exc),
                cgd_path=str(path) if path.suffix.lower() == ".cgd" else None,
            )

    return results


def _sort_key(value: str):
    '''
    Create a stable sort key for node labels.

    Numeric labels are sorted numerically before non-numeric labels.

    **parameters:**
        value: str
            Node label.

    **returns:**
        tuple
    '''
    try:
        return (0, int(value))
    except ValueError:
        return (1, value)


def _canonicalize_edges(edges: Sequence[Edge3D], decimals: int = 8) -> List[Edge3D]:
    '''
    Canonicalize relaxed edges for stable hashing.

    The function rounds coordinates, normalizes endpoint order within each edge,
    and removes duplicates.

    **parameters:**
        edges: list
            List of relaxed edges.

        decimals: int
            Number of decimal places used when rounding floats.

    **returns:**
        list
            Sorted, deduplicated canonical edges.
    '''
    canon = set()
    for p, q in edges:
        rp = _round_coord(p, decimals)
        rq = _round_coord(q, decimals)
        canon.add((rp, rq) if rp <= rq else (rq, rp))
    return sorted(canon, key=lambda item: (item[0], item[1]))


def _node_degrees_from_edges(comp: RelaxedComponent) -> Dict[Coord3D, int]:
    '''
    Build fallback node degrees from relaxed edge endpoints.

    **parameters:**
        comp: RelaxedComponent
            Parsed relaxed component.

    **returns:**
        python dictionary
            Mapping:
                coordinate -> degree
    '''
    degrees: Dict[Coord3D, int] = {}
    for p, q in comp.edges:
        degrees[p] = degrees.get(p, 0) + 1
        degrees[q] = degrees.get(q, 0) + 1
    return degrees


def crystal2_text_from_component(
    comp: RelaxedComponent,
    *,
    name: str,
    include_edge_centers: bool = True,
) -> str:
    '''
    Return CRYSTAL2 CGD text from a parsed relaxed component.

    **parameters:**
        comp: RelaxedComponent
            Parsed relaxed component.

        name: str
            Name written to the CRYSTAL2 block.

        include_edge_centers: bool
            If True, write edge centres as commented metadata lines.

    **returns:**
        str
            CRYSTAL2 CGD text.
    '''
    node_items = sorted(comp.nodes.items(), key=lambda kv: _sort_key(kv[0]))
    fallback_degrees = _node_degrees_from_edges(comp)
    a, b, c, alpha, beta, gamma = comp.cell

    lines = [
        "# Generated from Systre stdout (relaxed) -> CRYSTAL2",
        f"# Component: {comp.component_index}",
        f"# Systre-dimension: {comp.dimension}",
        f"# TD10: {comp.td10 if comp.td10 is not None else 'N/A'}",
    ]
    if comp.rcsr_name:
        lines.append(f"# RCSR: {comp.rcsr_name}")

    lines.extend([
        "CRYSTAL",
        f"  NAME {name}",
        f"  GROUP {comp.ideal_space_group or 'P1'}",
        f"  CELL {a:.6f} {b:.6f} {c:.6f} {alpha:.6f} {beta:.6f} {gamma:.6f}",
    ])

    for idx, (label, (x, y, z)) in enumerate(node_items, start=1):
        degree = comp.coordination_number.get(label)
        if degree is None:
            degree = fallback_degrees.get((x, y, z), 0)
        lines.append(f"  NODE {idx} {degree}  {x:.6f} {y:.6f} {z:.6f}")

    for (x1, y1, z1), (x2, y2, z2) in _canonicalize_edges(comp.edges, decimals=8):
        lines.append(f"  EDGE  {x1:.6f} {y1:.6f} {z1:.6f}   {x2:.6f} {y2:.6f} {z2:.6f}")

    if include_edge_centers:
        for x, y, z in sorted(comp.edge_centers, key=lambda xyz: (xyz[0], xyz[1], xyz[2])):
            lines.append(f"  # EDGE_CENTER  {x:.6f} {y:.6f} {z:.6f}")

    lines.append("END")
    return "\n".join(lines) + "\n"


def crystal2_fallback_from_periodic_graph(periodic_graph_cgd: str, *, name: str = "net") -> str:
    '''
    Convert a PERIODIC_GRAPH CGD into a CRYSTAL2 wrapper.

    This function is intended as a fallback when Systre did not emit relaxed
    coordinates, but the pipeline still needs a CRYSTAL2-like text block.

    **parameters:**
        periodic_graph_cgd: str
            Input PERIODIC_GRAPH CGD text.

        name: str
            Name written to the CRYSTAL2 block.

    **returns:**
        str
            CRYSTAL2 CGD text.
    '''
    out = [
        "# Fallback: CRYSTAL2 written from input PERIODIC_GRAPH (no relaxed output detected)",
        "CRYSTAL",
        f"  NAME {name}",
        "  GROUP P1",
    ]

    for line in periodic_graph_cgd.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if stripped.upper().startswith(("PERIODIC_GRAPH", "ID ", "NAME ", "GROUP ")):
            continue
        out.append("  " + stripped)

    out.append("END")
    return "\n".join(out) + "\n"


def main(argv: Optional[Sequence[str]] = None) -> int:
    '''
    Command-line entry point for topology identification and relaxed export.

    **parameters:**
        argv: list, optional
            Command-line arguments.

    **returns:**
        int
            Exit status code.
    '''
    parser = argparse.ArgumentParser(
        description="Identify topology with Systre and optionally export relaxed CRYSTAL2 CGD text."
    )
    parser.add_argument("input", help="Structure path or CGD path/text.")
    parser.add_argument("--input-is-cgd", action="store_true", help="Treat input as CGD path or CGD text.")
    parser.add_argument("--method", default="all_node", help="TopologyExtractor method used when generating CGD.")
    parser.add_argument("--name", default="net", help="Graph name used when generating CGD.")
    parser.add_argument("--java", default=None, help="Java executable or path. If omitted, jdk4py is tried first.")
    parser.add_argument("--timeout", type=int, default=30, help="Timeout in seconds.")
    parser.add_argument("--xmx", default="1024m", help="Java heap size, for example 1024m.")
    parser.add_argument("--component-index", type=int, default=None, help="Relaxed component index to export or hash.")
    parser.add_argument("--print-summary", action="store_true", help="Print normalized topology summary as JSON.")
    parser.add_argument("--print-crystal2", action="store_true", help="Print CRYSTAL2 CGD text.")
    parser.add_argument("--write-crystal2", default=None, help="Write CRYSTAL2 CGD text to this file.")
    parser.add_argument("--hash-topology", action="store_true", help="Print relaxed-topology hash.")
    parser.add_argument("--hash-crystal2", action="store_true", help="Print CRYSTAL2 text hash.")
    parser.add_argument("--fallback-to-input-cgd", action="store_true", help="Fallback to original CGD when relaxed output is unavailable.")
    parser.add_argument("--include-edge-centers", action="store_true", help="Include commented edge-centre metadata in CRYSTAL2 output.")

    args = parser.parse_args(argv)

    runner = SystreTopology(
        args.input,
        input_is_cgd=args.input_is_cgd,
        method=args.method,
        name=args.name,
        java=args.java,
        timeout_s=args.timeout,
        xmx=args.xmx,
    )

    result = runner.identify()
    print(result.topology)

    if args.print_summary:
        print(json.dumps(runner.topology_summary(component_index=args.component_index), indent=2, sort_keys=True))

    if args.hash_topology:
        print(runner.topology_hash(component_index=args.component_index))

    if args.print_crystal2:
        text = runner.crystal2_text(
            component_index=args.component_index,
            include_edge_centers=args.include_edge_centers,
            fallback_to_input_cgd=args.fallback_to_input_cgd,
        )
        if text is None:
            raise SystemExit("No CRYSTAL2 CGD text could be generated.")
        print(text, end="")

    if args.write_crystal2:
        runner.write_crystal2(
            args.write_crystal2,
            component_index=args.component_index,
            include_edge_centers=args.include_edge_centers,
            fallback_to_input_cgd=args.fallback_to_input_cgd,
        )

    if args.hash_crystal2:
        print(
            runner.crystal2_hash(
                component_index=args.component_index,
                include_edge_centers=args.include_edge_centers,
                fallback_to_input_cgd=args.fallback_to_input_cgd,
            )
        )

    return 0


__all__ = [
    "RelaxedComponent",
    "SystreResult",
    "SystreTopology",
    "identify_topology",
    "identify_topology_batch",
    "extract_relaxed_components",
    "parse_relaxed_component",
    "parse_systre_topology",
    "crystal2_text_from_component",
    "crystal2_fallback_from_periodic_graph",
    "find_java",
    "run_systre_on_cgd_path",
    "systre_command",
    "main",
]


if __name__ == "__main__":
    raise SystemExit(main())
