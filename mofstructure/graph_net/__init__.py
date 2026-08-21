#!/usr/bin/env python3
'''
Pure-Python identification and analysis of periodic nets.

`graph_net` replaces the parts of the Systre pipeline that required a Java
runtime. Given the quotient graph of a net it computes a canonical key, looks
that key up against the RCSR archive, derives the maximal symmetry the net can
achieve, and writes a geometric realisation as CGD. Nothing in the path uses
the unit cell, the space group or any relaxed geometry of the structure the
net came from: identification is combinatorial, and symmetry and geometry are
*derived from* the topology rather than supplied to it.

**the pieces.**
    `periodic_graph`  the quotient graph, lattice reduction, components
    `primitive`       reduction to the primitive cell
    `barycentric`     exact equilibrium placement and stability
    `canonical`       the canonical key
    `archive`         RCSR lookup
    `invariants`      coordination sequences, TD10, point and vertex symbols
    `symmetry`        maximal space group and intrinsic chirality
    `embedding`       geometric realisation and CGD output
    `lattice`         exact integer lattice arithmetic underneath all of it

**on keys.**
The canonical form is *not* byte-compatible with Systre's. The published
algorithm leaves several tie-breaks open and the conventions in the shipped
archive differ from the printed ones, so rather than reverse-engineer an
unspecified implementation this package uses its own canonical form and
re-keys the RCSR archive with it. What is guaranteed, and tested, is the
property that matters: one key for every representation of a net, and
different keys for different nets. A key produced here must not be compared
against a key produced by Systre.

**references:**
    Delgado-Friedrichs, O. & O'Keeffe, M. (2003). Acta Cryst. A59, 351-360.
        doi:10.1107/S0108767303012017
    Delgado-Friedrichs, O. (2005). Discrete Comput. Geom. 33, 67-81.
        doi:10.1007/s00454-004-1147-x
    Delgado-Friedrichs, O. & O'Keeffe, M. (2005). J. Solid State Chem. 178,
        2480-2485. doi:10.1016/j.jssc.2005.06.011
    O'Keeffe, M., Peskov, M. A., Ramsden, S. J. & Yaghi, O. M. (2008).
        Acc. Chem. Res. 41, 1782-1789. doi:10.1021/ar800124u
'''
from __future__ import annotations

from collections.abc import Sequence

from mofstructure.graph_net.archive import KEY_VERSION, key_hash, lookup
from mofstructure.graph_net.archive import describe as describe_key
from mofstructure.graph_net.barycentric import UnstableNetError
from mofstructure.graph_net.canonical import canonical_form, canonical_key
from mofstructure.graph_net.embedding import ideal_embedding, to_cgd
from mofstructure.graph_net.invariants import describe
from mofstructure.graph_net.periodic_graph import PeriodicGraph
from mofstructure.graph_net.primitive import primitive
from mofstructure.graph_net.symmetry import space_group

__version__ = "0.1.0"


def identify(
    graph: PeriodicGraph,
    descriptors: bool = True,
    symmetry: bool = True,
    prefer: Sequence[str] | None = None,
) -> list[dict[str, object]]:
    '''
    Identify every component of a net and describe it.

    A framework may deconstruct into several disjoint nets, as an
    interpenetrated structure does, so each connected component is handled
    separately and the result is a list. A component that cannot be keyed is
    reported with its reason rather than dropped, because silently omitting a
    component is how an interpenetrated structure gets mistaken for a simple
    one.

    **parameters:**
        - graph: PeriodicGraph
            Net to identify; it may be disconnected.

        - descriptors: bool
            Compute coordination sequences and point and vertex symbols. These
            are the useful output for a net with no RCSR name, and the most
            expensive part of the call.

        - symmetry: bool
            Derive the maximal space group and intrinsic chirality.

        - prefer: sequence of str, optional
            Archives the headline `topology` may come from, most preferred
            first; passed straight through to `archive.describe`.

    **returns:**
        list
            One dictionary per component, each carrying `periodicity`,
            `n_vertices`, `n_edges`, `interpenetration`, the canonical `key`
            with its `key_hash` and `key_version`, the headline name
            `topology` with the `topology_source` it came from and every
            known name in `names`, all of which are empty for a net no
            archive has named. When requested, `descriptors` and `symmetry`
            follow. A component that could not be keyed carries `error`
            instead of a key, which happens when it is unstable or has no
            periodicity at all.
    '''
    results: list[dict[str, object]] = []
    for component in graph.components():
        entry: dict[str, object] = {
            "periodicity": component.periodicity(),
            "n_vertices": component.n_vertices,
            "n_edges": component.n_edges,
            "interpenetration": component.covering_multiplicity(),
        }
        try:
            key = canonical_key(component)
        except UnstableNetError as exc:
            entry["error"] = f"unstable net: {exc}"
            entry["collisions"] = exc.collisions
            results.append(entry)
            continue
        except Exception as exc:  # noqa: BLE001 - reported, never swallowed
            entry["error"] = f"{type(exc).__name__}: {exc}"
            results.append(entry)
            continue

        # The key is the identification. Whether anyone has named the net it
        # denotes is a separate question, and a miss is an ordinary answer:
        # two structures with the same key have the same topology whether or
        # not that topology is in an archive, which is what makes them
        # comparable. So the key and its digest are always reported, and the
        # names are reported beside them with the source they came from.
        entry["key"] = key
        try:
            entry.update(describe_key(key, prefer=prefer))
        except FileNotFoundError:
            entry["topology"] = None
            entry["topology_source"] = None
            entry["names"] = {}
            entry["key_hash"] = key_hash(key)
            entry["key_version"] = KEY_VERSION
            entry["warning"] = "lookup table missing; run tools/build_graph_net_archive.py"

        if descriptors:
            entry["descriptors"] = describe(component)
        if symmetry:
            try:
                entry["symmetry"] = space_group(component)
            except Exception as exc:  # noqa: BLE001
                entry["symmetry"] = {"error": f"{type(exc).__name__}: {exc}"}
        results.append(entry)
    return results


__all__ = [
    "PeriodicGraph",
    "UnstableNetError",
    "canonical_form",
    "canonical_key",
    "describe",
    "ideal_embedding",
    "KEY_VERSION",
    "describe_key",
    "identify",
    "key_hash",
    "lookup",
    "primitive",
    "space_group",
    "to_cgd",
]
