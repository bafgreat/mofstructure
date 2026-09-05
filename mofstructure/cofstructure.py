#!/usr/bin/python
'''
COF linkage detection, building-unit extraction and topology construction.

Named finders identify supported linkage patterns. The bridge fallback cuts
unlike-element bonds in short acyclic chains between ring cores. Bond
perception removes excess contacts before ring detection and retains lattice
offsets for topology construction.

Cuts lie outside chemical rings, including attachments to boron and triazine
rings. Ring-forming linkages such as imide, dioxin, phenazine, benzoxazole,
benzimidazole, benzothiazole, thiazole and quinoline are unsupported. The
resulting fragments retain open valences.

Units with at least three connections become vertices; ditopic units are
contracted into edges.
'''
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"

import logging
from collections import Counter
from dataclasses import dataclass, field
from collections.abc import Sequence

import numpy as np
from ase.atoms import Atoms
from ase.data import covalent_radii

from mofstructure import generate_cgd, mofdeconstructor


logger = logging.getLogger(__name__)

Bond = tuple[int, int]

# A C=C is around 1.34 A and a C-C joining two sp2 carbons around 1.48 A, so
# anything below this is taken as the double bond of an olefin linkage.
OLEFIN_MAX_BOND_LENGTH = 1.42

# A framework graph closes pore rings as well as chemical ones: COF-1 closes a
# 14-membered B4C8O2 cycle alongside its 6-membered boroxine. The guards below
# mean chemical rings only, so larger ones are not looked for at all. Real COF
# monomers close rings of at most 6; the smallest pore ring is far above 8.
MAX_CHEMICAL_RING_SIZE = 8

# ASE's default skin adds 0.6 A to pair cutoffs, admitting contacts that can
# create false rings. This empirical heavy-atom cutoff allows C-C to 1.75 A;
# hydrogen contacts and isolated atoms are handled separately below.
MAX_BOND_RADII_FACTOR = 1.15


@dataclass
class COFPerception:
    '''
    Connectivity, ring membership and hydrogen counts of one structure.

    The linkage finders all need the same few lookups, and ring perception is
    the expensive part, so it is done once and passed around. The connectivity
    comes from ``compute_ase_neighbour_with_offsets`` because that is also what
    ``generate_cgd.kept_bond_graph`` uses downstream; perceiving bonds one way
    here and another way there would let a bond be cut that the topology code
    never saw.

    **attributes:**
        atoms: the ASE atoms object

        graph: atom index -> list of bonded atom indices

        bond_matrix: adjacency matrix matching ``graph``

        bond_offsets: (i, j) -> list of lattice offsets of that bond

        symbols: chemical symbol per atom index

        rings: list of rings, each a list of atom indices

        atom_in_ring: atom index -> True when the atom is in any ring

        atom_to_rings: atom index -> indices of the rings containing it

        heavy: atom index -> its non-hydrogen neighbours

        n_hydrogen: atom index -> how many hydrogens it carries
    '''
    atoms: Atoms
    graph: dict[int, list[int]]
    bond_matrix: np.ndarray
    bond_offsets: dict[Bond, list[tuple[int, int, int]]]
    symbols: list[str]
    rings: list[list[int]]
    atom_in_ring: dict[int, bool]
    atom_to_rings: dict[int, list[int]]
    heavy: dict[int, list[int]] = field(default_factory=dict)
    n_hydrogen: dict[int, int] = field(default_factory=dict)

    def in_ring(self, index: int) -> bool:
        '''True when the atom belongs to at least one ring.'''
        return mofdeconstructor.is_atom_in_ring(int(index), self.atom_in_ring)

    def share_ring(self, first: int, second: int) -> bool:
        '''True when both atoms lie in one and the same ring.'''
        return bool(
            set(self.atom_to_rings.get(int(first), []))
            & set(self.atom_to_rings.get(int(second), []))
        )

    def degree(self, index: int) -> int:
        '''Number of bonded neighbours, hydrogens included.'''
        return len(self.graph.get(int(index), []))

    def is_terminal_oxygen(self, index: int) -> bool:
        '''True for a carbonyl oxygen: one heavy neighbour and no hydrogen.'''
        return (
            self.symbols[int(index)] == 'O'
            and len(self.heavy[int(index)]) == 1
            and self.n_hydrogen[int(index)] == 0
        )

    def distance(self, first: int, second: int) -> float:
        '''Bond length across the periodic boundary.'''
        return float(self.atoms.get_distance(int(first), int(second), mic=True))


def _shortest_cycle_through_bond(graph, start, goal, max_size):
    '''
    Return the atoms of the smallest ring closed by the bond ``start``-``goal``.

    A breadth-first search from ``start`` to ``goal`` that is forbidden to use
    the bond itself finds the shortest alternative route between them, and that
    route plus the bond is the smallest ring the bond lies on.

    **parameters:**
        - graph: atom index -> list of bonded atom indices

        - start, goal: the two atoms of the bond

        - max_size: largest ring to look for

    **returns:**
        list of atom indices, or None when the bond lies on no ring this small
    '''
    parents = {start: None}
    frontier = [start]
    depth = 0

    while frontier and depth < max_size - 1:
        depth += 1
        next_frontier = []
        for node in frontier:
            for neighbour in graph[node]:
                if node == start and neighbour == goal:
                    # the bond being tested, and any second periodic image of it
                    continue
                if neighbour in parents:
                    continue
                parents[neighbour] = node
                if neighbour == goal:
                    cycle = []
                    current = goal
                    while current is not None:
                        cycle.append(current)
                        current = parents[current]
                    return cycle
                next_frontier.append(neighbour)
        frontier = next_frontier

    return None


def small_rings(graph, max_size: int = MAX_CHEMICAL_RING_SIZE) -> list[list[int]]:
    '''
    Find the chemical rings of a framework graph.

    ``networkx.minimum_cycle_basis`` returns the right rings but costs about a
    minute on a seven-hundred-atom cell, and ``cycle_basis`` is instant but
    returns fundamental cycles, which on a fused aromatic system are not the
    small rings at all. Taking the smallest ring through each bond instead is
    the ring set chemistry means, and it is a bounded breadth-first search per
    bond rather than a cycle basis of the whole graph.

    A ring here is a ring of atom indices. A bond that crosses a periodic
    boundary is folded into the same index space, so in a cell small enough for
    the framework to wrap around within ``max_size`` bonds this would report a
    ring that is not one. No COF cell is anywhere near that small.

    **parameters:**
        - graph: atom index -> list of bonded atom indices

        - max_size: largest ring to keep

    **returns:**
        list of rings, each a list of atom indices
    '''
    rings = {}
    for start in sorted(graph):
        for goal in graph[start]:
            if goal <= start:
                continue
            cycle = _shortest_cycle_through_bond(graph, start, goal, max_size)
            if cycle is not None:
                rings.setdefault(frozenset(cycle), cycle)
    return list(rings.values())


def prune_perceived_bonds(
    ase_atom: Atoms,
    graph,
    bond_matrix,
    bond_offsets,
    max_radii_factor: float = MAX_BOND_RADII_FACTOR,
):
    '''
    Remove excess contacts from ASE neighbour-list connectivity.

    Hydrogen retains its nearest heavy neighbour, or its nearest contact if no
    heavy neighbour exists. Heavy-heavy bonds are retained up to
    ``max_radii_factor * (r_i + r_j)``. The shortest attachment is restored for
    atoms left isolated, including heavy atoms left without a heavy neighbour.

    The attachment rule preserves long substituent bonds, such as the 1.67 A
    methoxy C-O bond in TpOMe-PaNO2.

    **parameters:**
        - ase_atom: ASE atoms object

        - graph: atom index -> list of bonded atom indices

        - bond_matrix: adjacency matrix matching ``graph``

        - bond_offsets: (i, j) -> list of lattice offsets of that bond

        - max_radii_factor: longest bond kept, as a multiple of r_i + r_j

    **returns:**
        the same three structures with the spurious bonds removed
    '''
    positions = ase_atom.get_positions()
    cell = np.asarray(ase_atom.get_cell())
    numbers = ase_atom.get_atomic_numbers()
    radii = covalent_radii[numbers]

    # (i, j, offset) -> length, for i <= j so each bond is measured once
    lengths = {}
    for (first, second), offsets in bond_offsets.items():
        if second < first:
            continue
        for offset in offsets:
            shift = np.asarray(offset, dtype=float) @ cell
            lengths[(first, second, tuple(int(x) for x in offset))] = float(
                np.linalg.norm(positions[second] + shift - positions[first])
            )

    shortest = {}
    nearest_heavy = {}
    for key, distance in lengths.items():
        first, second, _ = key
        for index, partner in ((first, second), (second, first)):
            if index not in shortest or distance < shortest[index][1]:
                shortest[index] = (key, distance)
            if numbers[partner] == 1:
                continue
            if index not in nearest_heavy or distance < nearest_heavy[index][1]:
                nearest_heavy[index] = (key, distance)

    def is_bond(key):
        first, second, _ = key
        if numbers[first] == 1 or numbers[second] == 1:
            # Prefer a heavy neighbour to avoid pairing nearby C-H hydrogens.
            return all(
                numbers[index] != 1
                or nearest_heavy.get(index, shortest[index])[0] == key
                for index in (first, second)
            )
        return (
            lengths[key]
            <= max_radii_factor * (radii[first] + radii[second])
        )

    kept = {key for key in lengths if is_bond(key)}
    # Restore attachments lost to the distance cutoff.
    attached = {index for key in kept for index in key[:2]}
    heavy_attached = {
        index for key in kept for index in key[:2]
        if numbers[key[0]] != 1 and numbers[key[1]] != 1
    }
    for index, (key, _) in shortest.items():
        if index not in attached:
            kept.add(key)
    for index, (key, _) in nearest_heavy.items():
        if numbers[index] != 1 and index not in heavy_attached:
            kept.add(key)

    new_offsets = {}
    new_graph = {index: [] for index in graph}
    new_matrix = np.zeros_like(bond_matrix)
    for first, second, offset in sorted(kept):
        back = tuple(-x for x in offset)
        new_offsets.setdefault((first, second), []).append(offset)
        new_graph[first].append(second)
        new_matrix[first, second] = 1
        new_matrix[second, first] = 1
        if (second, first, back) in kept and second != first:
            # already counted from the other end of the same bond
            continue
        new_offsets.setdefault((second, first), []).append(back)
        if second != first:
            new_graph[second].append(first)

    return new_graph, new_matrix, new_offsets


def perceive(ase_atom: Atoms) -> COFPerception:
    '''
    Build the shared perception of one structure.

    **parameters:**
        - ase_atom: ASE atoms object

    **returns:**
        COFPerception
    '''
    graph, bond_matrix, bond_offsets = \
        mofdeconstructor.compute_ase_neighbour_with_offsets(ase_atom)
    graph, bond_matrix, bond_offsets = prune_perceived_bonds(
        ase_atom, graph, bond_matrix, bond_offsets
    )
    graph = {int(k): [int(v) for v in vals] for k, vals in graph.items()}
    rings = small_rings(graph, MAX_CHEMICAL_RING_SIZE)

    atom_in_ring = {index: False for index in graph}
    atom_to_rings = {index: [] for index in graph}
    for ring_index, ring in enumerate(rings):
        for atom in ring:
            atom_in_ring[int(atom)] = True
            atom_to_rings[int(atom)].append(ring_index)

    symbols = list(ase_atom.get_chemical_symbols())
    heavy = {}
    n_hydrogen = {}
    for index, neighbours in graph.items():
        heavy[index] = [j for j in neighbours if symbols[j] != 'H']
        n_hydrogen[index] = len(neighbours) - len(heavy[index])

    return COFPerception(
        atoms=ase_atom,
        graph=graph,
        bond_matrix=bond_matrix,
        bond_offsets=bond_offsets,
        symbols=symbols,
        rings=rings,
        atom_in_ring=atom_in_ring,
        atom_to_rings=atom_to_rings,
        heavy=heavy,
        n_hydrogen=n_hydrogen,
    )


def find_imine_bonds(perception: COFPerception, **_options) -> list[Bond]:
    '''
    Find C=N bonds in imine, ketimine, hydrazone and azine linkages.

    Candidate nitrogen must be acyclic, hydrogen-free and bonded to two heavy
    atoms. Candidate carbon must be acyclic, have two or three heavy neighbours
    and at most one hydrogen. This permits CIFs with missing methine hydrogens;
    the nitrogen coordination requirement excludes terminal nitriles.

    Cuts retain nitrogen on the amine-derived fragment. In a hydrazone the N-N
    bond remains intact; in an azine both C=N bonds are cut. Ring atoms are
    excluded to preserve pyridine, bipyridine and porphyrin cores.

    **parameters:**
        - perception: COFPerception

    **returns:**
        list of (nitrogen, carbon) bonds to cut
    '''
    bonds = []
    for index, symbol in enumerate(perception.symbols):
        if symbol != 'N' or perception.in_ring(index):
            continue
        if perception.n_hydrogen[index] != 0:
            continue
        neighbours = perception.heavy[index]
        if len(neighbours) != 2:
            continue

        candidates = [
            j for j in neighbours
            if perception.symbols[j] == 'C'
            and not perception.in_ring(j)
            and len(perception.heavy[j]) in (2, 3)
            and perception.n_hydrogen[j] <= 1
        ]
        if not candidates:
            continue
        if len(candidates) > 1:
            # Both neighbours look like they could carry the double bond, so
            # fall back on the bond length: the C=N is the shorter of the two.
            candidates.sort(key=lambda j: perception.distance(index, j))
            first, second = candidates[0], candidates[1]
            if abs(perception.distance(index, first)
                   - perception.distance(index, second)) < 0.05:
                logger.info(
                    "cofstructure: nitrogen %d has two equally plausible C=N "
                    "partners, leaving it uncut.", index
                )
                continue
        bonds.append((index, candidates[0]))
    return bonds


def find_ketoenamine_bonds(perception: COFPerception, **_options) -> list[Bond]:
    '''
    Find the C-N of a beta-ketoenamine linkage.

    An imine made from 1,3,5-triformylphloroglucinol tautomerises: the phenol
    proton moves onto the nitrogen, the ring loses its aromaticity and keeps
    three exocyclic ketones, and what was a C=N becomes a C-N single bond with
    the double bond shifted into the ring::

           O
           ||
        ---C---C(ring)=CH-NH-Ar
              /
             C(=O)

    The nitrogen therefore carries a hydrogen and no double bond, so
    ``find_imine_bonds`` cannot see it. TpPa, TpBD and the rest of the Tp family
    are among the most reported COFs, so missing them would leave a large part
    of the literature undeconstructed.

    The signature required is a hydrogen-bearing nitrogen bonded to an exocyclic
    sp2 carbon, which is in turn bonded to a ring carbon whose ring neighbour
    carries a carbonyl oxygen. A plain diarylamine has the same nitrogen but not
    the beta carbonyl, so it is left alone.

    **parameters:**
        - perception: COFPerception

    **returns:**
        list of (nitrogen, carbon) bonds to cut
    '''
    bonds = []
    for index, symbol in enumerate(perception.symbols):
        if symbol != 'N' or perception.in_ring(index):
            continue
        if perception.n_hydrogen[index] != 1:
            continue
        neighbours = perception.heavy[index]
        if len(neighbours) != 2:
            continue
        if any(perception.symbols[j] != 'C' for j in neighbours):
            continue

        for carbon in neighbours:
            if perception.in_ring(carbon) or perception.degree(carbon) != 3:
                continue
            ring_carbons = [
                j for j in perception.heavy[carbon]
                if j != index
                and perception.symbols[j] == 'C'
                and perception.in_ring(j)
            ]
            if len(ring_carbons) != 1:
                continue
            if _has_neighbouring_ring_carbonyl(perception, ring_carbons[0]):
                bonds.append((index, carbon))
                break
    return bonds


def _has_neighbouring_ring_carbonyl(perception: COFPerception, carbon: int) -> bool:
    '''
    True when a ring neighbour of ``carbon`` carries a carbonyl oxygen.

    This is the beta-keto half of the beta-ketoenamine signature, and it is what
    separates the tautomerised Tp core from an ordinary secondary aromatic
    amine.
    '''
    for ring_index in perception.atom_to_rings.get(int(carbon), []):
        ring = set(perception.rings[ring_index])
        for neighbour in perception.heavy[int(carbon)]:
            if neighbour not in ring or perception.symbols[neighbour] != 'C':
                continue
            for candidate in perception.heavy[neighbour]:
                if candidate in ring:
                    continue
                if perception.is_terminal_oxygen(candidate):
                    return True
    return False


def find_azo_bonds(perception: COFPerception, **_options) -> list[Bond]:
    '''
    Find the N-N bond in aryl azo, azoxy and azodioxy linkages.

    Each nitrogen must be acyclic, hydrogen-free and bonded to one nitrogen and
    one ring carbon. Terminal oxygen substituents are allowed and remain on the
    nitrogen-containing fragments. This includes the azodioxy linkage in NPN
    frameworks. Azines have exocyclic carbon neighbours and are handled by
    ``find_imine_bonds``.

    **parameters:**
        - perception: COFPerception

    **returns:**
        list of (nitrogen, nitrogen) bonds to cut
    '''
    bonds = []
    for index, symbol in enumerate(perception.symbols):
        if symbol != 'N' or perception.in_ring(index):
            continue
        for partner in perception.heavy[index]:
            if partner <= index or perception.symbols[partner] != 'N':
                continue
            if perception.in_ring(partner):
                continue
            if not _is_aryl_azo_nitrogen(perception, index):
                continue
            if not _is_aryl_azo_nitrogen(perception, partner):
                continue
            bonds.append((index, partner))
    return bonds


def _is_aryl_azo_nitrogen(perception: COFPerception, index: int) -> bool:
    '''
    Check for hydrogen-free nitrogen bonded to one N and one ring C.

    Terminal oxygens are ignored when counting neighbours. Nitro groups fail
    the nitrogen-neighbour requirement.
    '''
    if perception.n_hydrogen[index] != 0:
        return False
    neighbours = [
        j for j in perception.heavy[index]
        if not perception.is_terminal_oxygen(j)
    ]
    if len(neighbours) != 2:
        return False
    carbons = [
        j for j in neighbours
        if perception.symbols[j] == 'C' and perception.in_ring(j)
    ]
    nitrogens = [j for j in neighbours if perception.symbols[j] == 'N']
    return len(carbons) == 1 and len(nitrogens) == 1


def find_amide_bonds(perception: COFPerception, **_options) -> list[Bond]:
    '''
    Find the C(=O)-N of an amide linkage.

    The carbonyl carbon must be acyclic and carry one terminal oxygen, one
    nitrogen and one carbon. Two guards keep this from over-cutting: a ring
    nitrogen means a lactam or a cyclic imide, where one cut would not separate
    anything, and a nitrogen bonded to another nitrogen means a hydrazide, whose
    C(=O)-N belongs inside the monomer of a hydrazone COF rather than being the
    linkage.

    **parameters:**
        - perception: COFPerception

    **returns:**
        list of (carbon, nitrogen) bonds to cut
    '''
    bonds = []
    for index, symbol in enumerate(perception.symbols):
        if symbol != 'C' or perception.in_ring(index):
            continue
        if perception.degree(index) != 3:
            continue
        neighbours = perception.heavy[index]
        if len(neighbours) != 3:
            continue

        oxygens = [j for j in neighbours if perception.is_terminal_oxygen(j)]
        nitrogens = [j for j in neighbours if perception.symbols[j] == 'N']
        carbons = [j for j in neighbours if perception.symbols[j] == 'C']
        if len(oxygens) != 1 or len(nitrogens) != 1 or len(carbons) != 1:
            continue

        nitrogen = nitrogens[0]
        if perception.in_ring(nitrogen):
            continue
        if any(perception.symbols[j] == 'N' for j in perception.heavy[nitrogen]):
            continue
        bonds.append((index, nitrogen))
    return bonds


def find_boron_linkage_bonds(perception: COFPerception, **_options) -> list[Bond]:
    '''
    Find aryl attachments to boroxine, boronate ester and borazine rings.

    Check all atoms in each boron-containing ring so that N-aryl attachments in
    borazines are included. Only bonds leaving the ring system are cut; fused
    catechol rings remain intact. Free boronic acids are excluded.

    **parameters:**
        - perception: COFPerception

    **returns:**
        list of (ring atom, carbon) bonds to cut
    '''
    bonds = []
    seen = set()
    for ring in perception.rings:
        if not any(perception.symbols[i] == 'B' for i in ring):
            continue
        ring_set = set(ring)
        for index in ring:
            for neighbour in perception.heavy[index]:
                if neighbour in ring_set:
                    continue
                if perception.symbols[neighbour] != 'C':
                    continue
                if perception.share_ring(index, neighbour):
                    continue
                key = (min(index, neighbour), max(index, neighbour))
                if key in seen:
                    continue
                seen.add(key)
                bonds.append((index, neighbour))
    return bonds


def find_triazine_bonds(perception: COFPerception, **_options) -> list[Bond]:
    '''
    Find the C(triazine)-C(aryl) of a covalent triazine framework.

    Trimerising three nitriles builds a C3N3 ring, which is a three-connected
    node of the net, so the cut is made at the bonds that attach the monomers to
    it and never inside it.

    **parameters:**
        - perception: COFPerception

    **returns:**
        list of (ring carbon, aryl carbon) bonds to cut
    '''
    bonds = []
    for ring in perception.rings:
        if len(ring) != 6:
            continue
        counts = Counter(perception.symbols[i] for i in ring)
        if counts.get('C') != 3 or counts.get('N') != 3 or len(counts) != 2:
            continue

        ring_set = set(ring)
        carbons = [i for i in ring if perception.symbols[i] == 'C']
        alternating = all(
            all(
                perception.symbols[j] == 'N'
                for j in perception.heavy[i] if j in ring_set
            )
            for i in carbons
        )
        if not alternating:
            continue

        for carbon in carbons:
            for neighbour in perception.heavy[carbon]:
                if neighbour in ring_set:
                    continue
                if perception.symbols[neighbour] == 'C':
                    bonds.append((carbon, neighbour))
    return bonds


def find_olefin_bonds(
    perception: COFPerception,
    *,
    olefin_requires_nitrile: bool = True,
    **_options,
) -> list[Bond]:
    '''
    Find the C=C of an olefin-linked (sp2 carbon) COF.

    A Knoevenagel condensation joins an aryl acetonitrile to an aldehyde and
    leaves ``Ar-CH=C(CN)-Ar'``. The trouble is that a stilbene-cored monomer has
    exactly the same local environment as this linkage, and no amount of graph
    inspection can tell a linkage vinylene from a monomer vinylene. The nitrile
    is what disambiguates in practice, so it is required by default; set
    ``olefin_requires_nitrile=False`` for an aldol-derived olefin COF that
    carries no nitrile, accepting that a stilbene monomer will then be cut in
    half.

    **parameters:**
        - perception: COFPerception

        - olefin_requires_nitrile: bool
            Require one of the two carbons to carry a nitrile group.

    **returns:**
        list of (carbon, carbon) bonds to cut
    '''
    bonds = []
    for index, symbol in enumerate(perception.symbols):
        if symbol != 'C' or perception.in_ring(index):
            continue
        if perception.degree(index) != 3:
            continue
        for partner in perception.heavy[index]:
            if partner <= index or perception.symbols[partner] != 'C':
                continue
            if perception.in_ring(partner) or perception.degree(partner) != 3:
                continue
            if perception.distance(index, partner) > OLEFIN_MAX_BOND_LENGTH:
                continue
            if not (_has_ring_carbon_neighbour(perception, index, partner)
                    and _has_ring_carbon_neighbour(perception, partner, index)):
                continue
            if olefin_requires_nitrile and not (
                _bears_nitrile(perception, index)
                or _bears_nitrile(perception, partner)
            ):
                continue
            bonds.append((index, partner))
    return bonds


def _has_ring_carbon_neighbour(
    perception: COFPerception, index: int, exclude: int
) -> bool:
    '''True when the atom is attached to a ring carbon other than ``exclude``.'''
    return any(
        j != exclude
        and perception.symbols[j] == 'C'
        and perception.in_ring(j)
        for j in perception.heavy[index]
    )


def _bears_nitrile(perception: COFPerception, index: int) -> bool:
    '''True when the atom carries a nitrile substituent.'''
    for neighbour in perception.heavy[index]:
        if perception.symbols[neighbour] != 'C':
            continue
        others = [j for j in perception.heavy[neighbour] if j != index]
        if len(others) != 1:
            continue
        nitrogen = others[0]
        if (perception.symbols[nitrogen] == 'N'
                and len(perception.heavy[nitrogen]) == 1
                and perception.n_hydrogen[nitrogen] == 0):
            return True
    return False


def find_hydrazide_bonds(perception: COFPerception, **_options) -> list[Bond]:
    '''
    Find the C-N of a hydrazone linkage drawn without its double bond.

    A hydrazone COF condenses an aldehyde with a hydrazide to give

        Ar-CH=N-NH-CO-Ar'

    and `find_imine_bonds` already cuts that CH=N. It can only do so when the
    nitrogen carries no hydrogen, which is what an imine nitrogen looks like.
    Deposited structures are not always drawn that way: in TpTe-1-COF every
    nitrogen carries one hydrogen and the C-N and N-N distances are 1.51 and
    1.48 A, both single bonds, so the file describes the reduced hydrazide
    rather than the hydrazone. The connectivity is nonetheless unambiguous and
    the linkage is still acyclic, so the same cut applies.

    The N-N pair is what identifies it. One nitrogen carries the carbonyl and
    stays with the hydrazide monomer; the other carries the carbon that used
    to be the aldehyde, and that is the bond to cut. Requiring the partner
    nitrogen keeps ordinary amines and amides out, and requiring the cut
    carbon to be free of a terminal oxygen keeps the cut off the hydrazide
    side, where `find_amide_bonds` has already declined to act for the same
    reason.

    **parameters:**
        - perception: COFPerception

    **returns:**
        list of (nitrogen, carbon) bonds to cut
    '''
    bonds = []
    for index, symbol in enumerate(perception.symbols):
        if symbol != 'N' or perception.in_ring(index):
            continue
        # An imine nitrogen carries no hydrogen, and `find_imine_bonds`
        # already cuts it. This finder exists only for the reduced drawing
        # that finder cannot see, so the two must not both claim the same
        # bond: without this guard a properly drawn hydrazone is cut twice.
        if perception.n_hydrogen[index] == 0:
            continue
        neighbours = perception.heavy[index]
        if len(neighbours) != 2:
            continue
        partners = [j for j in neighbours if perception.symbols[j] == 'N']
        carbons = [j for j in neighbours if perception.symbols[j] == 'C']
        if len(partners) != 1 or len(carbons) != 1:
            continue
        if perception.in_ring(partners[0]):
            continue
        carbon = carbons[0]
        if perception.in_ring(carbon):
            continue
        # A terminal oxygen on the carbon means this is the hydrazide side,
        # whose C(=O)-N belongs inside the monomer.
        if any(perception.is_terminal_oxygen(j) for j in perception.heavy[carbon]):
            continue
        bonds.append((index, carbon))
    return bonds


def find_ester_bonds(perception: COFPerception, **_options) -> list[Bond]:
    '''
    Find the O-C(aryl) of an ester linkage.

    Condensing a carboxylic acid with a phenol gives Ar-CO-O-Ar'. COF-120 is
    built this way: its carbonyl carbon sits 1.221 A from a terminal oxygen
    and 1.403 A from a bridging one, the two distances of a carbonyl and a
    single C-O.

    Of the two bonds either side of the bridging oxygen, the cut is made at
    the one to the aromatic ring rather than at the one to the carbonyl. Both
    separate the same two skeletons and so give the same net, but they divide
    the atoms differently, and this division is the one that leaves whole
    monomers behind: the acid keeps both of its oxygens, as a carboxyl group
    that can be recognised as such, and the other vertex is the aryl ring
    itself. Cutting beside the carbonyl instead would hand one of the acid's
    oxygens to the phenol and leave neither fragment a recognisable
    functional group.

    The carbonyl carbon must be acyclic: inside a ring it is a lactone, where
    one cut opens the ring rather than separating two monomers, which is the
    ring-forming case this module refuses rather than guesses at.

    **parameters:**
        - perception: COFPerception

    **returns:**
        list of (oxygen, aryl carbon) bonds to cut
    '''
    bonds = []
    for index, symbol in enumerate(perception.symbols):
        if symbol != 'O' or perception.in_ring(index):
            continue
        neighbours = perception.heavy[index]
        if len(neighbours) != 2:
            continue

        carbonyls = [
            j for j in neighbours
            if perception.symbols[j] == 'C'
            and not perception.in_ring(j)
            and any(perception.is_terminal_oxygen(k)
                    for k in perception.heavy[j])
        ]
        aryls = [
            j for j in neighbours
            if perception.symbols[j] == 'C' and perception.in_ring(j)
        ]
        if len(carbonyls) != 1 or len(aryls) != 1:
            continue
        bonds.append((index, aryls[0]))
    return bonds


def ring_systems(perception: COFPerception) -> dict[int, int]:
    '''
    Group fused and directly bonded rings into cores.

    Directly bonded rings, such as biphenyl, share one core. Acyclic groups
    between distinct cores are considered separately by ``acyclic_bridges``.

    **parameters:**
        - perception: COFPerception

    **returns:**
        ring atom index -> id of the core it belongs to
    '''
    parent = {}

    def find(index):
        parent.setdefault(index, index)
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(first, second):
        first, second = find(first), find(second)
        if first != second:
            parent[first] = second

    for index in perception.graph:
        if perception.in_ring(index):
            find(index)
    for index in list(parent):
        for neighbour in perception.heavy[index]:
            if perception.in_ring(neighbour):
                union(index, neighbour)

    return {index: find(index) for index in parent}


def acyclic_bridges(perception: COFPerception) -> list[dict]:
    '''
    Find acyclic heavy-atom groups connecting at least two ring cores.

    Prune terminal substituents from each group to obtain its backbone.
    Hydrogens are excluded. Return only backbones with at least two atoms.

    **parameters:**
        - perception: COFPerception

    **returns:**
        list of dicts with keys ``atoms``, ``backbone`` and ``cores``
    '''
    core_of = ring_systems(perception)
    acyclic = [
        index for index, symbol in enumerate(perception.symbols)
        if symbol != 'H' and not perception.in_ring(index)
    ]
    remaining = set(acyclic)

    bridges = []
    while remaining:
        stack = [remaining.pop()]
        group = set(stack)
        while stack:
            current = stack.pop()
            for neighbour in perception.heavy[current]:
                if neighbour in remaining:
                    remaining.discard(neighbour)
                    group.add(neighbour)
                    stack.append(neighbour)

        cores = {
            core_of[neighbour]
            for index in group
            for neighbour in perception.heavy[index]
            if neighbour in core_of
        }
        if len(cores) < 2:
            continue

        # Remove terminal branches without removing connections to ring cores.
        backbone = set(group)
        while True:
            leaves = {
                index for index in backbone
                if sum(
                    1 for neighbour in perception.heavy[index]
                    if neighbour in backbone or neighbour in core_of
                ) < 2
            }
            if not leaves:
                break
            backbone -= leaves

        if len(backbone) < 2:
            continue
        bridges.append(
            {'atoms': sorted(group), 'backbone': sorted(backbone),
             'cores': sorted(cores)}
        )
    return bridges


def find_bridge_bonds(
    perception: COFPerception,
    *,
    claimed: frozenset = frozenset(),
    max_bridge_length: int = 4,
    **_options,
) -> list[Bond]:
    '''
    Find candidate cuts in bridges not claimed by a named linkage finder.

    Cut backbone bonds between different elements, up to ``max_bridge_length``
    atoms per backbone. Skip a bridge if a named finder has already claimed a
    bond incident to its backbone, to avoid additional fragments.

    C-C bonds are excluded because monomer spacers can resemble linkages.
    Neither cut site may have more than one hydrogen, which excludes glycol
    spacers and reduced imines with explicit hydrogens. This guard does not infer
    bond order or hybridisation when hydrogens are missing.

    These rules describe local connectivity, not a unique condensation reaction
    or assignment of starting monomers.

    **parameters:**
        - perception: COFPerception

        - claimed: bonds already found, as a frozenset of ordered (i, j) pairs

        - max_bridge_length: longest backbone, in atoms, still called a linkage

    **returns:**
        list of (i, j) bonds to cut
    '''
    bonds = []
    for bridge in acyclic_bridges(perception):
        backbone = bridge['backbone']
        if len(backbone) > max_bridge_length:
            continue
        backbone_set = set(backbone)
        if any(
            (min(i, j), max(i, j)) in claimed
            for i in backbone_set for j in perception.heavy[i]
        ):
            continue

        for first in backbone:
            for second in perception.heavy[first]:
                if second <= first or second not in backbone_set:
                    continue
                if perception.symbols[first] == perception.symbols[second]:
                    continue
                if max(perception.n_hydrogen[first],
                       perception.n_hydrogen[second]) > 1:
                    continue
                bonds.append((first, second))
    return bonds


LINKAGE_FINDERS = {
    'imine': find_imine_bonds,
    'ketoenamine': find_ketoenamine_bonds,
    'hydrazide': find_hydrazide_bonds,
    'azo': find_azo_bonds,
    'amide': find_amide_bonds,
    'ester': find_ester_bonds,
    'boron': find_boron_linkage_bonds,
    'triazine': find_triazine_bonds,
    'olefin': find_olefin_bonds,
    'bridge': find_bridge_bonds,
}

# Run fallback finders after named linkage finders.
FALLBACK_LINKAGES = ('bridge',)

DEFAULT_LINKAGES = tuple(LINKAGE_FINDERS)


def find_CN_double_bonds(ase_atom, graph=None) -> list[list[int]]:
    '''
    Find the C=N of an imine-type linkage.

    Kept for callers written against the earlier signature. ``graph`` is
    ignored; ring perception is needed as well, so the structure is perceived
    afresh. Prefer :func:`find_imine_bonds`.

    **parameters:**
        - ase_atom: ASE atoms object

        - graph: ignored

    **returns:**
        list of [carbon, nitrogen] pairs
    '''
    perception = perceive(ase_atom)
    return [[carbon, nitrogen]
            for nitrogen, carbon in find_imine_bonds(perception)]


def cof_linkage_bonds(
    ase_atom: Atoms,
    linkages: Sequence[str] | None = None,
    perception: COFPerception | None = None,
    **options,
) -> dict[str, list[Bond]]:
    '''
    Locate every COF linkage bond, grouped by the linkage that produced it.

    **parameters:**
        - ase_atom: ASE atoms object

        - linkages: names of the linkages to search for, default all of
            ``DEFAULT_LINKAGES``

        - perception: a COFPerception to reuse, computed here when not given

        - options: passed to the finders, e.g. ``olefin_requires_nitrile``

    **returns:**
        python dictionary of linkage name -> list of (i, j) bonds, holding only
        the linkages that were actually found
    '''
    if perception is None:
        perception = perceive(ase_atom)

    names = DEFAULT_LINKAGES if linkages is None else tuple(linkages)
    unknown = [name for name in names if name not in LINKAGE_FINDERS]
    if unknown:
        raise ValueError(
            f"Unknown COF linkage(s) {unknown}. "
            f"Available: {sorted(LINKAGE_FINDERS)}."
        )

    # Named finders claim bonds before the fallback runs.
    ordered = [name for name in names if name not in FALLBACK_LINKAGES]
    ordered += [name for name in names if name in FALLBACK_LINKAGES]

    found = {}
    claimed = set()
    for name in ordered:
        if name in FALLBACK_LINKAGES:
            bonds = LINKAGE_FINDERS[name](
                perception, claimed=frozenset(claimed), **options
            )
        else:
            bonds = LINKAGE_FINDERS[name](perception, **options)
        if bonds:
            found[name] = bonds
            claimed |= {(min(i, j), max(i, j)) for i, j in bonds}
    return found


def secondary_building_units(
    ase_atom: Atoms,
    linkages: Sequence[str] | None = None,
    perception: COFPerception | None = None,
    **options,
):
    '''
    Deconstruct a COF into its building units by cutting the linkage bonds.

    Mirrors ``mofdeconstructor.secondary_building_units`` so the same downstream
    code can consume either. The porphyrin list is returned for that reason and
    is empty unless the COF carries a metallated porphyrin.

    **parameters:**
        - ase_atom: ASE atoms object

        - linkages: names of the linkages to search for, default all

        - perception: a COFPerception to reuse, computed here when not given

        - options: passed to the finders, e.g. ``olefin_requires_nitrile``

    **returns:**
        list_of_connected_components: building units as lists of atom indices

        bonds_to_break: the cut bonds as [i, j]

        porphyrin_checker: indices of metals sitting in a porphyrin

        all_regions: region id -> indices of the components in that region

        breaking_pairs: the cut bonds as [i, j, sx, sy, sz], carrying the
            lattice offset each bond crosses
    '''
    if perception is None:
        perception = perceive(ase_atom)
    porphyrin_checker = mofdeconstructor.metal_in_porphyrin2(
        ase_atom, perception.graph
    )

    found = cof_linkage_bonds(
        ase_atom, linkages=linkages, perception=perception, **options
    )
    if found:
        logger.info(
            "cofstructure: linkages found -> %s",
            ", ".join(f"{name}: {len(bonds)}" for name, bonds in found.items()),
        )
    else:
        logger.info(
            "cofstructure: no supported linkage found. Linkages that close a "
            "new ring across both monomers (imide, dioxin, phenazine, "
            "benzoxazole, benzimidazole, benzothiazole, thiazole, quinoline) "
            "are not handled."
        )

    seen = set()
    bonds_to_break = []
    breaking_pairs = []
    for bonds in found.values():
        for first, second in bonds:
            key = (min(first, second), max(first, second))
            if key in seen:
                continue
            seen.add(key)
            bonds_to_break.append([first, second])
            shift = mofdeconstructor.get_bond_shift(
                first, second, perception.bond_offsets
            )
            breaking_pairs.append([first, second, shift[0], shift[1], shift[2]])

    bond_matrix = perception.bond_matrix.copy()
    for first, second in bonds_to_break:
        bond_matrix[first, second] = 0
        bond_matrix[second, first] = 0

    new_graph = mofdeconstructor.matrix2dict(bond_matrix)
    try:
        list_of_connected_components = \
            mofdeconstructor.connected_components(new_graph)
    except Exception:
        import networkx as nx
        graph = nx.from_dict_of_lists(new_graph)
        list_of_connected_components = [
            list(component) for component in nx.connected_components(graph)
        ]

    all_regions = {}
    all_pm_structures = [
        sorted(ase_atom[component].symbols)
        for component in list_of_connected_components
    ]
    for i in range(len(all_pm_structures)):
        temp = []
        for j in range(len(all_pm_structures)):
            if all_pm_structures[i] == all_pm_structures[j]:
                temp.append(j)
        if temp not in all_regions.values():
            all_regions[i] = temp

    return [
        list_of_connected_components,
        bonds_to_break,
        porphyrin_checker,
        all_regions,
        breaking_pairs,
    ]


def unique_building_units(list_of_connected_components,
                          bonds_to_break,
                          ase_atom,
                          porphyrin_checker,
                          all_regions,
                          wrap_system=True,
                          cheminfo=True,
                          add_dummy=False
                          ):
    '''
    Return the unique COF building units.

    A thin wrapper over ``mofdeconstructor.find_unique_building_units`` that
    keeps only the organic units, since a COF has no metal cluster to separate.

    **parameters:**
        the outputs of :func:`secondary_building_units`, plus the wrapping,
        cheminformatics and dummy-atom switches of the underlying function

    **returns:**
        list of ASE atoms objects, one per unique building unit
    '''
    _, cof_linkers, _ = mofdeconstructor.find_unique_building_units(
        list_of_connected_components,
        bonds_to_break,
        ase_atom,
        porphyrin_checker,
        all_regions,
        wrap_system=wrap_system,
        cheminfo=cheminfo,
        add_dummy=add_dummy,
    )
    return cof_linkers


def _prune_terminal_nodes(edges):
    '''
    Drop vertices that only ever make one connection.

    An unreacted end group or a fragment left dangling by a defect is a vertex
    of degree one, and no periodic net can carry one: it collides in the
    barycentric placement, so no canonical key exists. Removing one can leave its neighbour with a single
    connection in turn, so this runs to a fixed point.

    **parameters:**
        - edges: list of (u, v, sx, sy, sz)

    **returns:**
        edges: the surviving edges

        removed: set of node ids that were dropped
    '''
    removed = set()
    current = list(edges)
    while True:
        degree = Counter()
        for u, v, _, _, _ in current:
            degree[u] += 1
            degree[v] += 1
        terminal = {node for node, count in degree.items() if count < 2}
        if not terminal:
            return current, removed
        removed |= terminal
        current = [
            edge for edge in current
            if edge[0] not in terminal and edge[1] not in terminal
        ]
        if not current:
            return current, removed


def _renumber(edges):
    '''
    Renumber node ids to be contiguous from zero.

    ``periodic_graph_cgd`` writes each node id straight out, so a gap left by a
    pruned or spliced vertex would appear in the net as an isolated vertex.

    **parameters:**
        - edges: list of (u, v, sx, sy, sz)

    **returns:**
        edges: the same edges with contiguous ids

        mapping: old node id -> new node id
    '''
    nodes = sorted({node for u, v, *_ in edges for node in (u, v)})
    mapping = {node: index for index, node in enumerate(nodes)}
    return (
        [(mapping[u], mapping[v], sx, sy, sz) for u, v, sx, sy, sz in edges],
        mapping,
    )


def cof_topology_graph(
    ase_atom: Atoms,
    *,
    collapse_ditopic: bool = True,
    linkages: Sequence[str] | None = None,
    **options,
):
    '''
    Build the periodic net of a COF from its building units.

    Every building unit left by the linkage cuts is a vertex and every cut is an
    edge between the two units it separated. There is no metal to mark out the
    nodes, so the distinction is made on connectivity: a unit joined at three or
    more points is a branch point and stays a vertex, and one joined at exactly
    two points is a connection rather than a node, so it is spliced into a
    direct edge. Leaving it in would subdivide the net, and a subdivided net has
    no RCSR name -- the archive holds hcb, not its subdivision.

    A unit that is periodic within itself (a fused ribbon) contributes its own
    translations as self-edges, and contact translations are reduced modulo
    those, for the same reason a rod SBU needs it: without it the answer depends
    on where the traversal happened to start.

    **parameters:**
        - ase_atom: ASE atoms object, guest-free

        - collapse_ditopic: bool
            Splice two-connected units into edges. Set False to keep every
            building unit as a vertex.

        - linkages: names of the linkages to search for, default all

        - options: passed to the finders, e.g. ``olefin_requires_nitrile``

    **returns:**
        edges: list of (u, v, sx, sy, sz) with contiguous zero-based node ids

        node_atoms: node id -> the atom indices that node represents
    '''
    perception = perceive(ase_atom)
    components, _, _, _, breaking_pairs = secondary_building_units(
        ase_atom, linkages=linkages, perception=perception, **options
    )
    if not breaking_pairs:
        raise ValueError(
            "cof_topology_graph found no supported COF linkage to cut."
        )
    if len(components) < 2:
        raise ValueError(
            "cof_topology_graph cut the structure but it stayed in one piece, "
            "so there is no net. This is what a ring-forming linkage does: one "
            "cut leaves the monomers joined through the rest of the fused ring."
        )

    kept_graph, kept_offsets = generate_cgd.kept_bond_graph(
        ase_atom, breaking_pairs, perception.graph, perception.bond_offsets
    )
    base_edges = generate_cgd.base_edges_with_shifts(
        ase_atom, components, breaking_pairs, kept_graph, kept_offsets
    )
    self_translations = generate_cgd.component_self_translations(
        components, kept_graph, kept_offsets
    )

    pair_lattice = {}
    edges = []
    for u, v, sx, sy, sz in base_edges:
        if (u, v) not in pair_lattice:
            pair_lattice[(u, v)] = generate_cgd._lattice_basis(
                list(self_translations.get(u, []))
                + list(self_translations.get(v, []))
            )
        edges.append(
            (u, v) + generate_cgd._reduce_modulo_lattice(
                (sx, sy, sz), pair_lattice[(u, v)]
            )
        )

    for component_id, translations in self_translations.items():
        for sx, sy, sz in translations:
            edges.append((component_id, component_id, sx, sy, sz))

    edges = generate_cgd.dedup_periodic_edges(edges)
    edges, pruned = _prune_terminal_nodes(edges)
    if pruned:
        logger.info(
            "cofstructure: %d building unit(s) make a single connection "
            "(unreacted end group or defect) and cannot be vertices of a net.",
            len(pruned),
        )
    if not edges:
        raise ValueError(
            "cof_topology_graph found no building unit with two or more "
            "connections, so there is no net."
        )

    if collapse_ditopic:
        degree = Counter()
        for u, v, *_ in edges:
            degree[u] += 1
            degree[v] += 1
        # Splicing every vertex out would leave edges between vertices that are
        # themselves gone. That happens only when nothing branches, which means
        # the structure is a chain rather than a framework, so the net is left
        # subdivided and the caller can see what it is.
        if any(count > 2 for count in degree.values()):
            ditopic = {node for node, count in degree.items() if count == 2}
            edges, spliced = generate_cgd.splice_two_connected_nodes(
                edges, ditopic
            )
            pruned |= spliced
        else:
            logger.info(
                "cofstructure: no building unit is joined at more than two "
                "points, so there is no branch point to splice towards."
            )

    edges, mapping = _renumber(edges)
    node_atoms = {
        new: {int(a) for a in components[old]}
        for old, new in mapping.items()
    }
    return edges, node_atoms


def cgd_cof(
    ase_atom: Atoms,
    *,
    name: str = "net",
    collapse_ditopic: bool = True,
    linkages: Sequence[str] | None = None,
    **options,
) -> str:
    '''
    Return the CGD PERIODIC_GRAPH of a COF net.

    **parameters:**
        - ase_atom: ASE atoms object, guest-free

        - name: CGD graph ID

        - collapse_ditopic: splice two-connected building units into edges

        - linkages: names of the linkages to search for, default all

        - options: passed to the finders, e.g. ``olefin_requires_nitrile``

    **returns:**
        CGD content as a string
    '''
    edges, _ = cof_topology_graph(
        ase_atom,
        collapse_ditopic=collapse_ditopic,
        linkages=linkages,
        **options,
    )
    return generate_cgd.periodic_graph_cgd(edges, name)
