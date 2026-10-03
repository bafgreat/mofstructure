#!/usr/bin/env python3
'''
One call from a structure to its topology.

Everything needed to identify the net of a framework already exists in this
package - `generate_cgd` deconstructs a structure into a quotient graph and
`graph_net` identifies one - but nothing joined them, so a caller had to
deconstruct, write a CGD, parse it back, build a `PeriodicGraph` and call
`identify` by hand, and had to know which deconstruction a given material
wanted before starting. This module is that missing layer.

**one answer shape for every material.**
A MOF, a COF and a zeolite are deconstructed along completely different
lines - at the metal cluster, at the linkage bond, at the bridging oxygen -
but what comes back is the same record in all three cases, because the thing
being reported is a net and a net does not remember what it was made of. That
is what makes the results comparable: a MOF and a COF built on the same net
answer with the same `key_hash`, which they would not if each material class
had its own reply.

**what is always there, and what is not.**
`key` is the identification. It exists for any locally stable net that
finishes computing, and it is unique: nets with the same key are the same
net, and 17409 independent nets were checked to give 17409 distinct keys.
A *name* is a different matter, and depends on whether an archive happens to
carry the net. Three situations produce no key at all, and each is reported
as itself rather than as a failure: a net that is unstable, for which no
canonical form exists at all; a fragment with no periodicity, which is not a
net; and a net too large to finish inside the budget.

**references:**
    Delgado-Friedrichs, O. & O'Keeffe, M. (2003). Acta Cryst. A59, 351-360.
        doi:10.1107/S0108767303012017
    O'Keeffe, M., Peskov, M. A., Ramsden, S. J. & Yaghi, O. M. (2008).
        Acc. Chem. Res. 41, 1782-1789.  (RCSR)
    Ramsden, S. J., Robins, V. & Hyde, S. T. (2009). Acta Cryst. A65, 81-108.
        (EPINET)
'''
from __future__ import annotations

import logging
import signal
from pathlib import Path
from collections.abc import Sequence

from ase.atoms import Atoms
from ase.io import read

from mofstructure import mofdeconstructor
from mofstructure.generate_cgd import (
    ZEOLITE_T_ELEMENTS,
    TopologyExtractor,
    tetrahedral_t_elements,
    zeolite_t_edges,
)
from mofstructure.graph_net import identify
from mofstructure.graph_net.archive import KEY_VERSION
from mofstructure.graph_net.periodic_graph import PeriodicGraph

logger = logging.getLogger(__name__)

InputLike = str | Atoms

#: The deconstruction each kind of framework is identified with by default.
#: A COF is cut at its linkage bond and a zeolite at its bridging oxygen, and
#: in neither case is there a second sensible place to cut.
DEFAULT_METHOD = {
    "zeolite": "zeol",
    "mof": "all_node",
    "cof": "cof",
}

#: The alternative nets a MOF has, in preference order. No one net is more
#: true than another; they answer different questions. `sbus` contracts each
#: metal cluster to a vertex, `all_node` keeps the atoms of a rod apart, so
#: MIL-53 reads as rna rather than pcu, and `ligand_cluster` gives an
#: incidence net sensitive to defects.
MOF_METHODS = ("sbus", "all_node", "single_node", "ligand_cluster")

#: Archives the headline `topology` is drawn from, by material, most
#: preferred first. MOFs and COFs are named from RCSR alone, so the column can
#: be grouped without asking which naming system a row used. A zeolite leads
#: with its IZA framework-type code and falls back to RCSR, since a dozen or
#: so framework nets carry an RCSR symbol with no IZA type assigned. `names`
#: carries every known name in every case.
NAME_PREFERENCE_BY_MATERIAL = {
    "mof": ("rcsr",),
    "cof": ("rcsr",),
    "zeolite": ("iza", "rcsr"),
}

#: Reported as the topology of a net that was identified but that no archive
#: names. A value rather than None: the key is present, so the net is
#: classified and two structures carrying it are comparable.
UNNAMED_TOPOLOGY = "unknown"

#: Reported as the topology of a structure that produced no net at all. With
#: `UNNAMED_TOPOLOGY` the field has three readings and no empty one: a name,
#: "unknown" for a net identified but unnamed, and this. How it failed is
#: `status`.
FAILED_TOPOLOGY = "error"

#: Elements that make a framework a MOF rather than a COF. A zeolite is
#: recognised before this is consulted, so the tetrahedral metals it shares
#: with this list - Zn, Co, Fe, Ti - do not cause a misreading.
_METALS = frozenset(
    {
        "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
        "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Rb", "Sr", "Y", "Zr", "Nb", "Mo",
        "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Cs", "Ba", "La",
        "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
        "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au",
        "Hg", "Tl", "Pb", "Bi", "Th", "U",
    }
)


class _Timeout(BaseException):
    '''
    Raised when identification exceeds the budget it was given.

    Derived from BaseException so that the per-component error handling in
    `graph_net.identify`, which catches Exception, does not record a timeout
    as an ordinary failure and report it as "unidentified".
    '''


def _alarm(signum, frame):  # noqa: ARG001
    '''
    SIGALRM handler bounding the cost of one structure.

    **parameters:**
        - signum: int
            Signal number.

        - frame: frame
            Interrupted stack frame.
    '''
    raise _Timeout()


def as_atoms(structure: InputLike) -> Atoms:
    '''
    Accept either a structure file or an ASE atoms object.

    **parameters:**
        - structure: str or ASE atoms object
            Path to anything ASE can read, or the atoms themselves.

    **returns:**
        ASE atoms object
    '''
    if isinstance(structure, Atoms):
        return structure
    return read(structure)


def _has_zeolite_framework(atoms: Atoms) -> bool:
    '''
    Check whether T-O-T contraction contains a three-periodic component.

    This detects frameworks with exchange cations or carbon-bearing guests that
    fail the element-only zeolite check.

    **parameters:**
        - atoms: ASE atoms object

    **returns:**
        bool
    '''
    t_elements = tetrahedral_t_elements(atoms)
    if not t_elements:
        return False
    try:
        n_vertices, edges, _ = zeolite_t_edges(atoms, t_elements=t_elements)
    except Exception:  # noqa: BLE001 - contraction failure leaves the type unresolved
        return False
    if not edges:
        return False
    net = PeriodicGraph.build(
        3, n_vertices, [(u, v, (sx, sy, sz)) for u, v, sx, sy, sz in edges]
    )
    return any(component.periodicity() == 3 for component in net.components())


def _metals_are_organic_bound(atoms: Atoms, graph: dict) -> bool:
    '''
    Check whether all metals match an organic building-unit environment.

    A metal qualifies if the porphyrin detector identifies it or all its
    neighbours are carbon. Every metal must qualify, so a porphyrin linker alone
    does not cause a MOF with separate metal nodes to be classified as a COF.
    Return False when no metals are present.

    **parameters:**
        - atoms: ASE atoms object

        - graph: atom index -> list of bonded atom indices

    **returns:**
        bool
    '''
    symbols = atoms.get_chemical_symbols()
    metals = [i for i, symbol in enumerate(symbols) if symbol in _METALS]
    if not metals:
        return False
    in_porphyrin = set(mofdeconstructor.metal_in_porphyrin2(atoms, graph))
    for index in metals:
        if index in in_porphyrin:
            continue
        # Treat carbon-only metal environments as part of an organic unit.
        neighbours = [symbols[j] for j in graph[index]]
        if neighbours and set(neighbours) == {"C"}:
            continue
        return False
    return True


def classify(structure: InputLike) -> str:
    '''
    Classify a structure as a zeolite, MOF or COF.

    Check zeolites first, using composition followed by T-O-T periodicity when
    needed. This allows exchange cations and carbon-bearing guests. Remaining
    structures with metals are MOFs unless all metals match the organic-bound
    criterion. Structures without metals are COFs.

    **parameters:**
        - structure: str or ASE atoms object
            Structure to classify.

    **returns:**
        str
            One of "zeolite", "mof" or "cof".
    '''
    atoms = as_atoms(structure)
    symbols = set(atoms.get_chemical_symbols())
    framework = symbols - {"H", "O"}
    if "C" not in symbols and framework and framework <= set(ZEOLITE_T_ELEMENTS):
        return "zeolite"

    # Perceive connectivity only when composition is insufficient.
    if "O" in symbols and symbols & set(ZEOLITE_T_ELEMENTS):
        if _has_zeolite_framework(atoms):
            return "zeolite"

    if symbols & _METALS:
        graph, _ = mofdeconstructor.compute_ase_neighbour(atoms)
        if _metals_are_organic_bound(atoms, graph):
            return "cof"
        return "mof"
    return "cof"


def quotient_graph(cgd_text: str, name: str = "net") -> PeriodicGraph:
    '''
    Read a CGD PERIODIC_GRAPH block back into a quotient graph.

    The CGD that `generate_cgd` writes is an explicit list of
    `(tail, head, shift)` edges rather than a geometric description, so this
    conversion loses nothing and no coordinates are involved.

    **parameters:**
        - cgd_text: str
            CGD content holding a PERIODIC_GRAPH block.

        - name: str
            Name to give the resulting graph.

    **returns:**
        PeriodicGraph

    **raises:**
        ValueError
            If the text holds no EDGES block, which is what an empty
            deconstruction produces.
    '''
    edges: list = []
    in_edges = False
    highest = 0
    for raw in cgd_text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        upper = line.upper()
        if upper.startswith("EDGES"):
            in_edges = True
            continue
        if upper.startswith("END"):
            in_edges = False
            continue
        if not in_edges:
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        tail, head = int(parts[0]), int(parts[1])
        highest = max(highest, tail, head)
        edges.append((tail - 1, head - 1, tuple(int(x) for x in parts[2:])))
    if not edges:
        raise ValueError("the deconstruction produced no edges")
    return PeriodicGraph.build(
        dim=len(edges[0][2]), n_vertices=highest, edges=edges, name=name
    )


def is_cgd(structure: InputLike) -> bool:
    '''
    Whether the input is a CGD periodic graph rather than a structure.

    Accepts a path ending in `.cgd` and CGD text passed directly, so a net
    written by this package can be read straight back without going through a
    structure file it never came from.

    **parameters:**
        - structure: str or ASE atoms object
            Candidate input.

    **returns:**
        bool
    '''
    if isinstance(structure, Atoms):
        return False
    text = str(structure)
    if "PERIODIC_GRAPH" in text or "CRYSTAL" in text:
        return True
    return text.lower().endswith(".cgd")


def _identify_into(
    result: dict[str, object],
    graph: PeriodicGraph,
    timeout: int | None,
    descriptors: bool,
    symmetry: bool,
) -> dict[str, object]:
    '''
    Identify a quotient graph and fold the answer into a result record.

    Shared by the structure path and the CGD path so that both report the
    same fields, the same statuses and the same summary rules.

    **parameters:**
        - result: python dictionary
            Record to fill in.

        - graph: PeriodicGraph
            Net to identify.

        - timeout: int or None
            Seconds allowed, None for no limit.

        - descriptors: bool
            Compute coordination sequences and symbols.

        - symmetry: bool
            Derive the maximal space group.

    **returns:**
        python dictionary
            The same record, completed. `topology` is never empty: a net
            that was keyed but that no archive names carries
            `UNNAMED_TOPOLOGY`, and one that produced no key carries
            `FAILED_TOPOLOGY`, so an identified net is never reported as
            though it had failed and neither is mistaken for the other.
    '''
    if timeout:
        signal.signal(signal.SIGALRM, _alarm)
        signal.setitimer(signal.ITIMER_REAL, timeout)
    try:
        result["components"] = identify(
            graph,
            descriptors=descriptors,
            symmetry=symmetry,
            # The material is already on the record when a structure was read,
            # and absent when a bare net was, which is the right default: a
            # CGD net that was never a structure has no material to restrict
            # its naming by.
            prefer=NAME_PREFERENCE_BY_MATERIAL.get(result.get("material")),
        )
    except _Timeout:
        result["status"] = "timeout"
        result["detail"] = f"identification exceeded {timeout}s"
        return result
    except Exception as exc:  # noqa: BLE001 - reported, never swallowed
        result["status"] = "error"
        result["detail"] = f"{type(exc).__name__}: {exc}"
        return result
    finally:
        if timeout:
            signal.setitimer(signal.ITIMER_REAL, 0)

    def agreed(field: str):
        '''The value all components share, or None if they differ.'''
        values = {
            component.get(field)
            for component in result["components"]
            if component.get(field)
        }
        return values.pop() if len(values) == 1 else None

    result["topology"] = agreed("topology")
    result["key"] = agreed("key")
    result["key_hash"] = agreed("key_hash")
    result["n_components"] = len(result["components"])

    # A component that produced no key has not been identified, whatever the
    # deconstruction managed. The commonest reason is an unstable net, which
    # has no canonical form at all, so saying "ok" here would report a
    # failure as a success.
    if not any(component.get("key") for component in result["components"]):
        errors = [component.get("error")
                  for component in result["components"]
                  if component.get("error")]
        result["status"] = "unidentified"
        result["detail"] = errors[0] if errors else (
            "no component produced a canonical key")

    # A net that was keyed but that no archive names has still been
    # identified. The test is whether any component was keyed, not whether the
    # record carries one: an interpenetrated structure whose components are
    # different nets has no single key, yet each component was identified.
    if any(component.get("key") for component in result["components"]):
        result["topology"] = result.get("topology") or UNNAMED_TOPOLOGY
    else:
        result["topology"] = FAILED_TOPOLOGY

    return result


def analyse(
    structure: InputLike,
    *,
    method: str = "auto",
    timeout: int | None = 300,
    descriptors: bool = False,
    symmetry: bool = False,
    name: str = "net",
) -> dict[str, object]:
    '''
    Identify the topology of a framework.

    **parameters:**
        - structure: str or ASE atoms object
            Path to anything ASE can read, or the atoms themselves.

        - method: str
            Deconstruction to use. "auto" classifies the structure first;
            otherwise any method `TopologyExtractor.build_cgd` accepts.

        - timeout: int or None
            Seconds allowed for the identification, None for no limit. The
            cost of the canonical form is driven by the number of vertices
            and by degenerate coordination rather than by the size of the
            structure, and a large net can take hours, so a bounded answer
            that says it ran out of time is more useful than an unbounded
            one that never returns.

        - descriptors: bool
            Compute coordination sequences and point and vertex symbols.
            These are the useful output for a net with no name, and the most
            expensive part of the call.

        - symmetry: bool
            Derive the maximal space group and intrinsic chirality.

        - name: str
            Name to carry through to the quotient graph.

    **returns:**
        python dictionary
            Keys `status`, `material`, `method`, `key_version`, `components`
            and, when every component agrees, the convenience fields
            `topology`, `key` and `key_hash`. `status` is "ok" when a canonical
            key was produced, and otherwise names what stopped it:
            "deconstruction_failed" when the framework could not be reduced to
            building units, "no_net" when the deconstruction left no periodic
            graph, "unidentified" when a graph was built but has no canonical
            form, which happens for an unstable net, and "timeout" or "error".
            A named net is not required for "ok". `topology` reads as one
            of three things and never as nothing: a name, `UNNAMED_TOPOLOGY`
            for a net that was identified but that no archive names, and
            `FAILED_TOPOLOGY` for a structure that produced no net, whatever
            `status` says stopped it. Which archive supplies a name depends
            on the material, following `NAME_PREFERENCE_BY_MATERIAL`. The
            per-component records are those `graph_net.identify` returns.
    '''
    result: dict[str, object] = {
        "status": "ok",
        "material": None,
        "method": method,
        "key_version": KEY_VERSION,
        # Every path out of this routine carries a topology, including the
        # ones that return before a net exists, so the field is seeded with
        # the failing value and overwritten once there is something to say.
        "topology": FAILED_TOPOLOGY,
        "key": None,
        "key_hash": None,
        "n_components": 0,
        "components": [],
    }
    # A CGD is already a quotient graph, so it needs no structure and no
    # deconstruction: identifying one is the natural way to ask about a net
    # that came from somewhere else, or to re-read a net this package wrote.
    if is_cgd(structure):
        result["method"] = "cgd"
        try:
            text = (
                str(structure)
                if "\n" in str(structure)
                else Path(structure).read_text()
            )
            graph = quotient_graph(text, name=name)
        except (ValueError, OSError) as exc:
            result["status"] = "no_net"
            result["detail"] = str(exc)
            return result
        return _identify_into(result, graph, timeout, descriptors, symmetry)

    try:
        atoms = as_atoms(structure)
    except Exception as exc:  # noqa: BLE001 - reported, never swallowed
        result["status"] = "error"
        result["detail"] = f"{type(exc).__name__}: {exc}"
        return result

    # The material is what the structure is, not what the caller asked for,
    # so it is determined the same way whether or not a method was named.
    material = classify(atoms)
    if method == "auto":
        method = DEFAULT_METHOD[material]
    result["material"] = material
    result["method"] = method

    try:
        cgd = TopologyExtractor(ase_atoms=atoms).build_cgd(method=method, name=name)
    except Exception as exc:  # noqa: BLE001
        result["status"] = "deconstruction_failed"
        result["detail"] = f"{type(exc).__name__}: {exc}"
        return result

    try:
        graph = quotient_graph(cgd, name=name)
    except ValueError as exc:
        result["status"] = "no_net"
        result["detail"] = str(exc)
        return result

    return _identify_into(result, graph, timeout, descriptors, symmetry)


def analyse_methods(
    structure: InputLike,
    methods: Sequence[str] | None = None,
    **options,
) -> dict[str, dict[str, object]]:
    '''
    Identify one structure under every deconstruction that suits it.

    This is a MOF question. A MOF has more than one defensible net and they
    answer different questions, so where a single answer is wanted `analyse`
    picks the conventional one and this reports them all side by side. A COF
    and a zeolite have one cut each and answer with that alone. The structure
    is read once and reused, since reading dominates the cost for small nets.

    **parameters:**
        - structure: str or ASE atoms object
            Structure to identify.

        - methods: sequence of str, optional
            Deconstructions to run. Defaults to `MOF_METHODS` for a MOF and
            to the single sensible method for a COF or a zeolite.

        - options:
            Passed through to `analyse`.

    **returns:**
        python dictionary
            Mapping method name -> the record `analyse` returns for it. A
            method that does not apply is present with its failing status
            rather than omitted, so the absence of an answer is visible.
    '''
    atoms = as_atoms(structure)
    if methods is None:
        material = classify(atoms)
        # Only a MOF has alternatives. A COF and a zeolite each have one cut,
        # so "every method" is that one method rather than a sweep that would
        # hand a zeolite to the COF linkage finder and get nothing back.
        methods = MOF_METHODS if material == "mof" else (DEFAULT_METHOD[material],)
    return {
        method: analyse(atoms, method=method, **options) for method in methods
    }


def analyse_many(
    structures: Sequence[InputLike],
    **options,
) -> list[dict[str, object]]:
    '''
    Identify a sequence of structures, one record each.

    A structure that fails is recorded and the run continues, since the
    common use is a directory of thousands where a handful will not
    deconstruct.

    **parameters:**
        - structures: sequence
            Paths or atoms objects.

        - options:
            Passed through to `analyse`.

    **returns:**
        list
            One record per structure, in input order, each carrying `source`.
    '''
    out = []
    for structure in structures:
        record = analyse(structure, **options)
        record["source"] = (
            structure if isinstance(structure, str) else getattr(
                structure, "info", {}
            ).get("name", "atoms")
        )
        out.append(record)
    return out


__all__ = [
    "DEFAULT_METHOD",
    "MOF_METHODS",
    "is_cgd",
    "analyse",
    "analyse_methods",
    "analyse_many",
    "as_atoms",
    "classify",
    "quotient_graph",
]
