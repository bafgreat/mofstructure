#!/usr/bin/env python3
'''
Lookup of a canonical key against the RCSR archive.

The Reticular Chemistry Structure Resource collects the named nets of
reticular chemistry, and the archive shipped with Systre stores one canonical
key per named net. Identification of a framework therefore reduces to a
dictionary lookup once the key has been computed.

`graph_net` uses its own canonical form rather than reproducing Gavrog's byte
sequence, so the archive is re-keyed once by `tools/build_graph_net_archive.py`
and the result cached in `db/graph_net_rcsr.json`. The original
`db/RCSRnets-2019-06-01.arc` remains the source of truth and is still shipped,
both because it carries the RCSR identifiers and because the re-keyed table
must be reproducible from it.

**references:**
    O'Keeffe, M., Peskov, M. A., Ramsden, S. J. & Yaghi, O. M. (2008).
        Acc. Chem. Res. 41, 1782-1789. doi:10.1021/ar800124u
    Delgado-Friedrichs, O. & O'Keeffe, M. (2003). Acta Cryst. A59, 351-360.
        doi:10.1107/S0108767303012017
'''
from __future__ import annotations

import hashlib
import json
from functools import lru_cache
from collections.abc import Sequence
from pathlib import Path

from mofstructure.filetyper import load_dict_msgpack

_DB_DIR = Path(__file__).resolve().parent.parent / "db"
SYSTRE_ARCHIVE = _DB_DIR / "RCSRnets-2019-06-01.arc"
GRAPH_NET_TABLE = _DB_DIR / "graph_net_rcsr.json"
MERGED_TABLE = _DB_DIR / "graph_net_archive.json"
PACKED_TABLE = _DB_DIR / "graph_net_archive.msgpack"

#: Version of the canonical form the shipped tables were built with. The
#: canonical form is this package's own convention, not Systre's, and several
#: of its tie-breaks are free choices; changing any of them silently
#: invalidates every stored key and every fingerprint derived from one.
#: Recording the version next to the key makes such a change detectable
#: instead of leaving it to be discovered by a mismatch much later.
KEY_VERSION = "graph_net/1"

#: Order in which archives are preferred when a net is named by more than
#: one. A name reticular chemistry has chosen is more informative than a
#: framework-type code, which is more informative than an index into a
#: systematic enumeration.
NAME_PREFERENCE = ("rcsr", "iza", "epinet")


def parse_systre_archive(path: Path | None = None) -> list[tuple[str, str]]:
    '''
    Read (rcsr_id, systre_key) pairs from a Systre `.arc` file.

    **parameters:**
        - path: pathlib.Path or None
            Archive file; the bundled RCSR archive when omitted.

    **returns:**
        list
            List of (rcsr_id, systre_key) pairs in file order.
    '''
    path = path or SYSTRE_ARCHIVE
    entries: list[tuple[str, str]] = []
    current: str | None = None
    for line in path.read_text().splitlines():
        if line.startswith("key"):
            current = line.split(None, 1)[1].strip()
        elif line.startswith("id") and current is not None:
            entries.append((line.split(None, 1)[1].strip(), current))
            current = None
    return entries


@lru_cache(maxsize=1)
def _table() -> dict[str, dict[str, str]]:
    '''
    Load and cache the lookup table and the provenance of each entry.

    The packed table is the one the package ships and is preferred. It holds
    the RCSR nets and, alongside them, the systematically enumerated EPINET
    nets that CrystalNets distributes, which between them name six times as
    many nets as the RCSR archive alone. The JSON tables the build tools
    write are intermediates, not shipped, so a checkout that has only those
    still works but is not the normal case.

    **returns:**
        python dictionary
            Mapping canonical key -> {archive: name}.
    '''
    if PACKED_TABLE.exists():
        payload = load_dict_msgpack(str(PACKED_TABLE))
        stored = payload.get("key_version")
        if stored and stored != KEY_VERSION:
            raise ValueError(
                f"{PACKED_TABLE.name} was packed under canonical form "
                f"{stored!r} but this build produces {KEY_VERSION!r}; the "
                "keys are not comparable, so repack the table"
            )
        return payload.get("entries", {})
    if MERGED_TABLE.exists():
        payload = json.loads(MERGED_TABLE.read_text())
        if payload.get("names"):
            return payload["names"]
        provenance = payload.get("provenance", {})
        return {
            key: {provenance.get(key, "rcsr"): name}
            for key, name in payload.get("keys", {}).items()
        }
    if not GRAPH_NET_TABLE.exists():
        raise FileNotFoundError(
            f"none of {PACKED_TABLE}, {MERGED_TABLE} or {GRAPH_NET_TABLE} "
            "is present; regenerate with "
            "python tools/build_graph_net_archive.py and pack it with "
            "python tools/pack_graph_net_archive.py"
        )
    payload = json.loads(GRAPH_NET_TABLE.read_text())
    return {key: {"rcsr": name} for key, name in payload.get("keys", {}).items()}


def key_hash(key: str, algorithm: str = "sha256") -> str:
    '''
    Stable fixed-width handle for a canonical key.

    The key itself is the identity of a net and should be stored, since only
    it can be looked up again if the archive later grows. A hash of it is
    convenient as a database column or join field, and it is a hash *of the
    key* rather than of a CGD file for a concrete reason: a CGD describes one
    representation, so the same net written with a different vertex order,
    origin or supercell hashes differently, while its canonical key - and
    therefore this digest - does not.

    **parameters:**
        - key: str
            Canonical key as produced by `canonical.canonical_key`.

        - algorithm: str
            Any name accepted by `hashlib.new`.

    **returns:**
        str
            `<version>:<algorithm>:<hexdigest>`, the version included so a
            digest made under a different canonical form is never mistaken
            for one made under this one.
    '''
    digest = hashlib.new(algorithm, key.encode("utf-8")).hexdigest()
    return f"{KEY_VERSION}:{algorithm}:{digest}"


def describe(
    key: str,
    prefer: Sequence[str] | None = None,
) -> dict[str, str | None]:
    '''
    Everything the archive can say about a canonical key.

    The name and its source are reported separately because they are
    different claims. An RCSR symbol says the net is one reticular chemistry
    has named; an EPINET symbol says it matches an entry in a systematic
    enumeration of nets obtained by reticulating hyperbolic tilings. Both are
    useful, and conflating them would overstate the second.

    A net with no name at all is an ordinary outcome, not a failure: the key
    still identifies it uniquely, which is what makes two structures
    comparable to each other whether or not anyone has named their topology.

    One net may be named by more than one archive: a zeolite framework
    carries an RCSR symbol and an IZA framework-type code at once, so ABW is
    `sra` to the former and `ABW` to the latter. `topology` therefore reports
    a single headline name and `topology_source` says which archive it came
    from, with `names` carrying every name known for the net. That shape
    keeps one field to read for "what is this", states the provenance
    explicitly, and takes a further archive without changing its structure -
    where one field per archive would need a new field and every caller
    would have to learn about it.

    Preference runs RCSR, then IZA, then EPINET: a name reticular chemistry
    has chosen says more than a framework code, which in turn says more than
    a position in a systematic enumeration.

    **parameters:**
        - key: str
            Canonical key as produced by `canonical.canonical_key`.

        - prefer: sequence of str, optional
            Archives the headline name may be drawn from, most preferred
            first. Defaults to `NAME_PREFERENCE`. A caller that knows what
            the material is can narrow it, so that a MOF is never headlined
            with a zeolite framework code; `names` still carries every name,
            so narrowing the headline discards nothing.

    **returns:**
        python dictionary
            Keys `key`, `key_hash`, `key_version`, `topology`,
            `topology_source` and `names`. The last three are None, None and
            an empty dictionary for a net no archive has named, which is an
            ordinary outcome rather than a failure.
    '''
    names = _table().get(key) or {}
    topology = None
    source = None
    for candidate in (NAME_PREFERENCE if prefer is None else tuple(prefer)):
        if candidate in names:
            topology = names[candidate]
            source = candidate
            break
    return {
        "key": key,
        "key_hash": key_hash(key),
        "key_version": KEY_VERSION,
        "topology": topology,
        "topology_source": source,
        "names": dict(names),
    }


def lookup(key: str) -> str | None:
    '''
    RCSR identifier for a canonical key, or None when the net is not named.

    A miss is a perfectly ordinary outcome and means the net is new to the
    archive rather than that anything failed; Systre reports the same
    situation as "structure is new for this run".

    **parameters:**
        - key: str
            Canonical key as produced by `canonical.canonical_key`.

    **returns:**
        str or None
    '''
    names = _table().get(key) or {}
    for candidate in NAME_PREFERENCE:
        if candidate in names:
            return names[candidate]
    return None


def table_size() -> int:
    '''
    Number of distinct canonical keys currently in the lookup table.

    **returns:**
        int
    '''
    return len(_table())


__all__ = [
    "KEY_VERSION",
    "NAME_PREFERENCE",
    "MERGED_TABLE",
    "PACKED_TABLE",
    "SYSTRE_ARCHIVE",
    "GRAPH_NET_TABLE",
    "describe",
    "key_hash",
    "parse_systre_archive",
    "lookup",
    "table_size",
]
