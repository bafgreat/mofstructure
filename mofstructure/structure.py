#!/usr/bin/python
'''
High level entry point to mofstructure.

MOFstructure wraps a single crystal structure and exposes the analyses that are
normally wanted together: guest removal, porosity, deconstruction into building
units, topology and open metal sites. Everything is computed on the guest free
system, so a structure only has to be read once.

    from mofstructure import structure

    mof = structure.MOFstructure(filename='UiO-66.cif')
    pores = mof.get_porosity()
    metal_sbus, organic_sbus = mof.get_sbu()

The individual algorithms live in mofdeconstructor, porosity and systre, and can
be called directly when more control is needed.
'''
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"
import os
from functools import reduce
import operator
import argparse
import pandas as pd
import shutil
import tempfile
import json
import numpy as np
import logging
from ase.io import read
from omsdetector_forked import MofCollection, mof
import mofstructure.mofdeconstructor as MOF_deconstructor
from mofstructure.porosity import zeo_calculation
import mofstructure.filetyper as read_write
from mofstructure.systre import SystreTopology
from mofstructure.generate_cgd import ligand_cluster_fingerprint

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


class MOFstructure(object):
    '''
    A single framework and the analyses that can be run on it.

    Give it either an ASE atoms object or a path to any file ASE can read. The
    methods return their results rather than writing them, so they can be mixed
    freely:

        mof = MOFstructure(filename='UiO-66.cif')
        pores = mof.get_porosity()
        metal_sbus, organic_sbus = mof.get_sbu()
        topology = mof.get_topology()

    Guest molecules are stripped internally before each analysis, so a structure
    that still contains solvent does not need cleaning up first.

    **parameters:**
        ase_atoms: ASE atoms object, used in preference to filename
        filename: path to a cif or any other ASE readable structure file
    '''

    def __init__(self,
                 ase_atoms=None,
                 filename=None
                 ):
        """
        A global method for extracting useful information
        """
        if ase_atoms is not None:
            self.ase_atoms = ase_atoms
        else:
            self.ase_atoms = read(filename)

    def remove_guest(self):
        '''
        Simple function to remove guest molecules in porous system.
        Note that this can work for any periodic system.

        **return:**
            ase_atoms (ase.Atoms): atom object with guest removed.
        '''
        index_non_guest = MOF_deconstructor.remove_unbound_guest(self.ase_atoms)
        return self.ase_atoms[index_non_guest]

    def get_sbu(self, wrap_system=True, cheminfo=True, add_dummy=False):
        """
        Extract the metal and linker secondary building units from the guest free system.

        **parameter**
            wrap_system (bool): If True, removes the effects of periodicity by merging the system.
            cheminfo (bool): If True, computes cheminformatic identifiers such as SMILES, InChI, and InChIKey.
            add_dummy (bool): If True, adds dummy atoms at the points of extension.

        **return:**
            metal_sbu (list): A list of unique metal secondary building units
            linker_sbu (list): A list of unique organic secondary building units
        """
        guest_free_atoms = self.remove_guest()
        connected_components, atoms_indices_at_breaking_point, porpyrin_checker, all_regions, _ = MOF_deconstructor.secondary_building_units(guest_free_atoms)
        if len(connected_components) > 0:
            metal_sbus, organic_sbus, _ = MOF_deconstructor.\
                find_unique_building_units(connected_components,
                                           atoms_indices_at_breaking_point,
                                           guest_free_atoms,
                                           porpyrin_checker,
                                           all_regions,
                                           wrap_system=wrap_system,
                                           cheminfo=cheminfo,
                                           add_dummy=add_dummy
                                           )
            return metal_sbus, organic_sbus
        else:
            logging.error(
                "Sorry, we were unable to successfully deconstruct your system."
                "Send us an email so that we can implement better rules for"
                "deconstructing your system."
            )
            return None

    def get_ligands(self, wrap_system=True, cheminfo=True, add_dummy=False):
        """
        Extract the metal and linker secondary building units from the guest free system.

        **parameter:**
            wrap_system (bool): If True, removes the effects of periodicity by merging the system.
            cheminfo (bool): If True, computes cheminformatic identifiers such as SMILES, InChI, and InChIKey.
            add_dummy (bool): If True, adds dummy atoms at the points of extension.

        **return:**
            metal_clusters (list): A list of unique metal atoms or clusters.
            organic_ligands (list): A list of unique organic ligands.
        """
        guest_free_atoms = self.remove_guest()
        connected_components, atoms_indices_at_breaking_point, porpyrin_checker, all_regions, _ = MOF_deconstructor.ligands_and_metal_clusters(guest_free_atoms)
        if len(connected_components) > 0:
            metal_clusters, organic_ligands, _ = MOF_deconstructor.\
                find_unique_building_units(connected_components,
                                           atoms_indices_at_breaking_point,
                                           guest_free_atoms,
                                           porpyrin_checker,
                                           all_regions,
                                           wrap_system=wrap_system,
                                           cheminfo=cheminfo,
                                           add_dummy=add_dummy
                                           )
            return metal_clusters, organic_ligands
        else:
            logging.error(
                "Sorry, we were unable to successfully deconstruct your system."
                "Send us an email so that we can implement better rules for"
                "deconstructing your system."
            )
            return None

    def get_porosity(self,
                     probe_radius=1.86,
                     number_of_steps=10000,
                     rad_file=None,
                     high_accuracy=True):
        '''
        A function to compute porosity data for a system.

        **parameters**

            probe_radius (float): Radius of the probe (default: 1.86).
            number_of_steps (int): Number of GCMC simulation cycles (default: 10000).
            high_accuracy (bool): If True, perform high-accuracy computations.
            rad_file: Optional file containing user defined atom radii. Must have the `.rad` extension

        **return:**
            pore (dict): A dictionary containing:
                - AV_Volume_fraction: Accessible volume void fraction.
                - AV_A^3: Accessible volume in A^2.
                - AV_cm^3/g: Accessible volume in cm³/g. This value is often infinite because it divides the computed
                    volume by Avogadro's number.
                - ASA_A^2: Accessible surface area in A^2.
                - ASA_m^2/cm^3: Accessible surface area in m2/cm3.
                - Number_of_channels: Number of channels (i.e., pores) present in the system.
                - LCD_A: The largest cavity diameter, defined as the diameter of the largest sphere that can be
                    inserted into the porous system without overlapping any atoms.
                - lfpd_A: The largest included sphere along the free sphere path, i.e., the largest sphere that can be
                    inserted into the pore.
                - PLD_A: The pore limiting diameter, defined as the largest sphere that can freely diffuse through the
                    porous network without overlapping any atoms.
        '''
        guest_free_atoms = self.remove_guest()
        if rad_file is None:
            pores = zeo_calculation(guest_free_atoms,
                                    probe_radius=probe_radius,
                                    number_of_steps=number_of_steps,
                                    high_accuracy=high_accuracy
                                    )
        else:
            pores = zeo_calculation(guest_free_atoms,
                                    probe_radius=probe_radius,
                                    number_of_steps=number_of_steps,
                                    high_accuracy=high_accuracy,
                                    rad_file=rad_file)
        return read_write.convert_numpy_types(pores)

    def get_topology(
        self,
        method="all_node",
        *,
        decimals=8,
        include_edge_centers=True,
        fallback_to_input_cgd=False,
    ):
        """
        Compute topology information for the guest-free system.

        **parameters:**
            method: str
                Topology extraction method passed to SystreTopology.

            decimals: int
                Number of decimal places used when hashing the relaxed topology payload.

            include_edge_centers: bool
                If True, include edge-center comments in the CRYSTAL2 text.

            fallback_to_input_cgd: bool
                If True, return a CRYSTAL2 wrapper from the input CGD when no
                relaxed component can be parsed from Systre output.

        **return:**
            python dictionary
                Mapping containing:
                    - topology
                    - dimension
                    - td10
                    - topology_hash
                    - cgd
        """
        guest_free_atoms = self.remove_guest()

        runner = SystreTopology(
            guest_free_atoms,
            method=method,
            name="net",
            keep_tmp=False,
        )

        result = runner.identify()
        comp = runner.best_component()

        if comp is None:
            return {
                "topology": result.topology,
                "dimension": None,
                "td10": None,
                "topology_hash": None,
                "cgd": runner.crystal2_text(
                    include_edge_centers=include_edge_centers,
                    fallback_to_input_cgd=fallback_to_input_cgd,
                ),
            }

        return {
            "topology": result.topology,
            "dimension": comp.dimension,
            "td10": comp.td10,
            "topology_hash": comp.topology_hash(decimals=decimals),
            "cgd": comp.crystal2_text(
                name=comp.rcsr_name or "net",
                include_edge_centers=include_edge_centers,
            ),
        }

    def get_ligand_cluster_fingerprint(self):
        """
        Describe how the ligands meet the metal clusters.

        Where ``get_topology`` asks Systre to name the net, this reads the
        deconstruction directly, so it answers for every framework, including
        the ones Systre leaves unnamed or refuses. It counts each ligand and
        cluster species, how many clusters each ligand bridges, and with what
        denticity, which is what makes it sensitive to defects: a missing
        linker lowers a cluster's connectivity, a linker hanging by one end is
        listed under ``terminal`` alongside its formula, and a carboxylate that
        has dropped from bridging to monodentate shows in the denticity
        histogram. It does not change when the atoms are listed in a different
        order, when the cell origin moves, or when the same crystal is given as
        a supercell.

        **return:**
            python dictionary
                Mapping containing:
                    - clusters
                    - ligands
                    - terminal
                    - refinement
                    - formula_units
                    - fingerprint_hash
        """
        return ligand_cluster_fingerprint(self.remove_guest())

    def draw_topology(self, method="all_node", supercell=(1, 1, 1),
                      filename=None, show=False, show_structure=True,
                      show_unit_cell=True, show_linker_sbu=True,
                      show_topology=False):
        '''
        Draw the complete topological network over the real framework geometry.

        Each node sits at the real centroid of the atoms it represents (metal
        atoms/clusters and the carboxyl and linker vertices), and edges follow
        the framework connectivity, including bonds that cross the cell. The
        figure is presented as a molecular structure rather than an axis-based
        plot: rotate, zoom and hover a node for its coordination.  Nodes in
        neighbouring periodic images are included whenever an edge reaches
        them, so every visible edge has a visible node at both ends.

        **parameters:**
            method: str
                Node definition: "sbus", "all_node", "single_node" or
                "ligand_cluster".

            supercell: tuple of three ints
                How many cells to draw along a, b, c. (1, 1, 1) draws one cell
                plus the edges leaving it.

            filename: str, optional
                If given, write the figure to this path. An .html file stays
                interactive; other extensions (.png, .pdf, ...) need the
                optional 'kaleido' package.

            show: bool
                If True, open the figure in a browser.

            show_structure: bool
                If True, draw the framework atoms and bonds behind the net.

            show_unit_cell: bool
                If True, draw the boundary of the displayed unit cells.

            show_linker_sbu: bool
                If True, show the centre-to-centre network produced by the
                selected topology method.  Its nodes and connections therefore
                change when the method changes.

            show_topology: bool
                If True, also show the abstract topology nodes and blue edges.
                It is False by default so the SBU–linker connectivity remains
                visually unambiguous.

        **returns:**
            plotly.graph_objects.Figure
        '''
        try:
            import plotly.graph_objects as go
        except ImportError as exc:
            raise ImportError(
                "draw_topology needs plotly. Install it with 'pip install plotly'."
            ) from exc

        from mofstructure.generate_cgd import net_geometry

        guest_free = self.remove_guest()
        method = method.lower().strip()
        positions, kinds, edges, cell = net_geometry(guest_free, method=method)

        degree = {nid: 0 for nid in positions}
        for u, v, *_ in edges:
            if u == v:
                degree[u] += 2
            else:
                degree[u] += 1
                degree[v] += 1

        reps = [
            (i, j, k)
            for i in range(int(supercell[0]))
            for j in range(int(supercell[1]))
            for k in range(int(supercell[2]))
        ]

        # Edge translations and node positions use the same periodic gauge.
        # Applying the stored translation preserves self-edges and distinct
        # connections to different images of the same node.
        edge_x, edge_y, edge_z = [], [], []
        # Keep every node instance touched by a displayed edge.  In particular,
        # an edge leaving the requested supercell must show its destination
        # node; otherwise the line appears to end in empty space and the
        # connectivity is ambiguous.
        visible_nodes = {
            (nid, tuple(cell_shift))
            for cell_shift in reps for nid in positions
        }
        for cell_shift in reps:
            base = np.array(cell_shift) @ cell
            for u, v, sx, sy, sz in edges:
                start = positions[u] + base
                edge_shift = np.array((sx, sy, sz), dtype=int)
                end = positions[v] + base + edge_shift @ cell
                if np.linalg.norm(end - start) < 1e-6:
                    continue
                visible_nodes.add((u, tuple(cell_shift)))
                visible_nodes.add((v, tuple(np.asarray(cell_shift) + edge_shift)))
                edge_x += [start[0], end[0], None]
                edge_y += [start[1], end[1], None]
                edge_z += [start[2], end[2], None]

        traces = []

        if show_structure:
            atom_x, atom_y, atom_z = [], [], []
            for translation in reps:
                shifted = guest_free.positions + np.array(translation) @ cell
                atom_x.extend(shifted[:, 0])
                atom_y.extend(shifted[:, 1])
                atom_z.extend(shifted[:, 2])
            traces.append(go.Scatter3d(
                x=atom_x, y=atom_y, z=atom_z, mode="markers",
                marker=dict(size=3, color="#64748b", opacity=0.20),
                hoverinfo="skip", name="framework atoms",
            ))

            atom_graph, _, atom_offsets = \
                MOF_deconstructor.compute_ase_neighbour_with_offsets(guest_free)
            bond_x, bond_y, bond_z = [], [], []
            for translation in reps:
                seen_bonds = set()
                base = np.array(translation) @ cell
                for i, neighbours in atom_graph.items():
                    for j in neighbours:
                        j = int(j)
                        for shift in atom_offsets.get((int(i), j), []):
                            shift = tuple(int(x) for x in shift)
                            reverse = (j, int(i), -shift[0], -shift[1], -shift[2])
                            key = (int(i), j, *shift)
                            canonical = min(key, reverse)
                            if canonical in seen_bonds:
                                continue
                            seen_bonds.add(canonical)
                            start = guest_free.positions[int(i)] + base
                            end = guest_free.positions[j] + base + np.array(shift) @ cell
                            bond_x += [start[0], end[0], None]
                            bond_y += [start[1], end[1], None]
                            bond_z += [start[2], end[2], None]
            traces.append(go.Scatter3d(
                x=bond_x, y=bond_y, z=bond_z, mode="lines",
                line=dict(color="#94a3b8", width=1), opacity=0.32,
                hoverinfo="skip",
                name="framework bonds",
            ))

        if show_unit_cell:
            cell_x, cell_y, cell_z = [], [], []
            corners = [(i, j, k) for i in (0, 1) for j in (0, 1) for k in (0, 1)]
            cell_edges = []
            for corner in corners:
                for axis in range(3):
                    neighbour = list(corner)
                    if neighbour[axis] == 0:
                        neighbour[axis] = 1
                        cell_edges.append((corner, tuple(neighbour)))
            for translation in reps:
                for a, b in cell_edges:
                    start = (np.array(translation) + np.array(a)) @ cell
                    end = (np.array(translation) + np.array(b)) @ cell
                    cell_x += [start[0], end[0], None]
                    cell_y += [start[1], end[1], None]
                    cell_z += [start[2], end[2], None]
            traces.append(go.Scatter3d(
                x=cell_x, y=cell_y, z=cell_z, mode="lines",
                line=dict(color="#94a3b8", width=2, dash="dot"),
                opacity=0.55,
                hoverinfo="skip", name="unit cell",
            ))

        if show_topology:
            traces.append(go.Scatter3d(
                x=edge_x, y=edge_y, z=edge_z, mode="lines",
                line=dict(color="#2563eb", width=5),
                hoverinfo="none", name="topology edges",
            ))

        if show_linker_sbu:
            # Use the selected method's graph.  A method-independent
            # ligand_cluster overlay made all four drawings look alike and hid
            # the contractions that define each topological representation.
            map_positions, map_kinds, map_edges = positions, kinds, edges

            contact_x, contact_y, contact_z = [], [], []
            visible_map_nodes = {
                (node, tuple(translation))
                for translation in reps for node in map_positions
            }
            contact_degree = {node: 0 for node in map_positions}
            for translation in reps:
                base = np.array(translation) @ cell
                for u, v, sx, sy, sz in map_edges:
                    edge_shift = np.array((sx, sy, sz), dtype=int)
                    destination = np.asarray(translation) + edge_shift
                    contact_degree[u] += 1 if translation == reps[0] else 0
                    contact_degree[v] += 1 if translation == reps[0] else 0
                    start = map_positions[u] + base
                    end = map_positions[v] + base + edge_shift @ cell
                    visible_map_nodes.add((u, tuple(translation)))
                    visible_map_nodes.add((v, tuple(destination)))
                    contact_x += [start[0], end[0], None]
                    contact_y += [start[1], end[1], None]
                    contact_z += [start[2], end[2], None]
            traces.append(go.Scatter3d(
                x=contact_x, y=contact_y, z=contact_z, mode="lines",
                line=dict(color="#008000", width=7, dash="solid"),
                hoverinfo="none", name=f"{method} connections",
            ))

            mapping_styles = {
                "metal": ("SBU / metal centres", "square", "#dc2626", 8),
                "organic": (
                    "organic / linker centres", "diamond", "#9333ea", 7
                ),
            }
            for kind, (label, symbol, color, size) in mapping_styles.items():
                mapped = sorted(
                    (node, shift)
                    for node, shift in visible_map_nodes
                    if map_kinds[node] == kind
                )
                if not mapped:
                    continue
                coords = np.array([
                    map_positions[node] + np.array(shift) @ cell
                    for node, shift in mapped
                ])
                text = [
                    f"{method}: {label[:-1]} {node}<br>connections "
                    f"{contact_degree[node]}"
                    for node, _ in mapped
                ]
                traces.append(go.Scatter3d(
                    x=coords[:, 0], y=coords[:, 1], z=coords[:, 2],
                    mode="markers",
                    marker=dict(
                        size=size, symbol=symbol, color=color,
                        line=dict(color="#ffffff", width=1),
                    ),
                    text=text, hoverinfo="text", name=label,
                ))

        if show_topology:
            palette = {"metal": "#fb8500", "organic": "#14213d"}
            for kind in ("metal", "organic"):
                nodes = [nid for nid in positions if kinds[nid] == kind]
                if not nodes:
                    continue
                instances = sorted(
                    (nid, shift) for nid, shift in visible_nodes if nid in nodes
                )
                coords = np.array([
                    positions[nid] + np.array(shift) @ cell
                    for nid, shift in instances
                ])
                text = [
                    f"{kind} node {nid}<br>coordination {degree[nid]}"
                    for nid, _ in instances
                ]
                traces.append(go.Scatter3d(
                    x=coords[:, 0], y=coords[:, 1], z=coords[:, 2],
                    mode="markers",
                    marker=dict(
                        size=10 if kind == "metal" else 7,
                        color=palette[kind],
                        line=dict(color="#000000", width=1),
                    ),
                    text=text, hoverinfo="text",
                    name=f"{kind} nodes",
                ))

        title = f"{method} net"
        topo = self.get_topology(method=method).get("topology")
        if topo:
            title = f"{title} — {topo}"

        fig = go.Figure(data=traces)
        fig.update_layout(
            title=dict(
                text=title,
                x=0.5,
                xanchor="center",
                font=dict(size=21, color="#17231b"),
            ),
            showlegend=True,
            paper_bgcolor="#f4f7f5",
            plot_bgcolor="#f4f7f5",
            font=dict(
                family="Inter, Avenir, Helvetica, Arial, sans-serif",
                color="#33413a",
            ),
            margin=dict(l=8, r=8, b=8, t=58),
            hoverlabel=dict(
                bgcolor="#ffffff",
                bordercolor="#008000",
                font=dict(color="#17231b", size=13),
            ),
            legend=dict(
                x=0.99,
                y=0.99,
                xanchor="right",
                yanchor="top",
                bgcolor="rgba(255,255,255,0.82)",
                bordercolor="rgba(23,35,27,0.18)",
                borderwidth=1,
                font=dict(size=12),
                itemsizing="constant",
            ),
            scene=dict(
                # A topology is a structure, not a chart.  Suppress axes,
                # ticks and grid planes while retaining Plotly's useful 3-D
                # rotation and hover interaction.
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                zaxis=dict(visible=False),
                bgcolor="#f4f7f5",
                aspectmode="data",
                camera=dict(
                    eye=dict(x=1.45, y=1.45, z=1.15),
                    projection=dict(type="orthographic"),
                ),
            ),
        )

        if filename:
            if str(filename).lower().endswith(".html"):
                fig.write_html(filename)
            else:
                fig.write_image(filename)
        if show:
            fig.show()
        return fig

    def get_oms(self):
        """
        Function to compute open metal sites
        """
        general_info = {}
        overlap = MOF_deconstructor.inter_atomic_distance_check(self.ase_atoms)
        general_info['has_overlapping_atoms'] = not overlap

        if len(self.ase_atoms) > 3000:
            print(f'The system size is very large: {len(self.ase_atoms)} atoms')
            print('Hence oms computation will take a while ')
            print('Thanks for your patience!!!')

        lattice = self.ase_atoms.get_cell().tolist()
        species = self.ase_atoms.get_chemical_symbols()
        coords = self.ase_atoms.get_positions()
        oms_structure = mof.MofStructure(lattice=lattice,
                                         species=species,
                                         coords=coords,
                                         coords_are_cartesian=True,
                                         name='oms_system'
                                         )

        coordination_sphere = oms_structure.metal_coord_spheres
        tmp_dir = tempfile.mkdtemp()
        try:
            oms_structure.analyze_metals(tmp_dir)

            json_files = [f for f in os.listdir(tmp_dir) if f.endswith('.json')]
            if not json_files:
                raise FileNotFoundError(f"No JSON file found in temporary directory: {tmp_dir}")


            json_file = os.path.join(tmp_dir, json_files[0])
            data = read_write.load_data(json_file)


            general_info["metals"] = data.get("metal_species")
            general_info["has_oms"] = data.get("has_oms")
            general_info["density"] = data.get("density")
            general_info["uc_volume"] = data.get("uc_volume")
            general_info["error_in_systems"] = data.get("problematic")


            metal_sites_dic = data.get("metal_sites")
            if metal_sites_dic is not None:
                unique_tuples = {
                    (d["metal"], d["number_of_linkers"], d["is_open"])
                    for d in metal_sites_dic
                }
                unique_metal = [
                    {
                        "metal": metal,
                        "coordination_number": number,
                        "is_open": is_open,
                        "enviroment": get_coordination_environment(coordination_sphere, metal, number)
                    }
                    for metal, number, is_open in unique_tuples
                ]
                oms_metal = [i.get("metal") for i in unique_metal if i["is_open"]]
                general_info["metal_info"] = unique_metal
                general_info["open_metals"] = oms_metal
        finally:
            shutil.rmtree(tmp_dir)


        return general_info



def get_coordination_environment(coordination_spheres, metal_element, coordination_number):
    """
    Returns the coordination environment (neighbors) of a specified metal element
    only if the number of neighbors equals the given coordination number.

    **parameters:**
        coordination_spheres (list): A list of metal coordination spheres (MetalSite objects).
        metal_element (str): The metal element symbol (e.g., 'Ag').
        coordination_number (int): The desired number of neighbors (coordination number).

    Returns:
        List[str]: A list of species symbols (neighbors) bonded to the metal center.
                   Returns an empty list if no matching environment is found.
    """
    metal_element = str(metal_element)

    for metal_site in coordination_spheres:
        # Convert the species list to strings.
        species_str = [str(s) for s in metal_site.species]

        # Early check: if the total number of atoms is not metal + neighbors, skip.
        if len(species_str) != coordination_number + 1:
            continue

        try:
            center_idx = species_str.index(metal_element)
        except ValueError:
            continue

        # Get the neighbors by removing the metal element.
        neighbors = species_str[:center_idx] + species_str[center_idx+1:]
        if len(neighbors) == coordination_number:
            return neighbors

    return []
