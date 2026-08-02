#!/usr/bin/python
'''
Stacking analysis of layered covalent organic frameworks.

Layered COFs are described by how consecutive sheets sit relative to one
another, which is measured here as the lateral offset between adjacent layers
and the interlayer spacing along the stacking direction. That separates
eclipsed AA stacking from the various staggered arrangements.
'''
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"
import argparse
import os
import glob
import numpy as np
from mofstructure.filetyper import write_json
from ase.io import read
import mofstructure.mofdeconstructor as MOF_deconstructor


def compute_cof_stacking(ase_atom):
    """
    A a simple function to compute the stacking pattern of COFs or
    layered materials like graphene

    **parameter:**
        ase_atom : ASE Atoms object
    **returns**
        layers : list of list wherei each list correspond to a layar
        lateral_offsets : list of list where each list contains the lateral offsets between two layers
        interlayer_height : list of list where each list contains the interlayer heights between two layers
    """
    indices = MOF_deconstructor.remove_unbound_guest(ase_atom)

    ase_atom = ase_atom[indices]

    graph, _ = MOF_deconstructor.compute_ase_neighbour(ase_atom)
    layers = []
    lateral_offsets = []
    interlayer_height = []
    components = MOF_deconstructor.connected_components(graph)
    if len(components) > 1:
        for i in range(len(components)-1):
            for j in range(1, len(components)):
                if i != j:
                    layers.append([i+1, j+1])

                    layer1_indices = components[i]
                    layer2_indices = components[j]

                    layer1_positions = ase_atom[layer1_indices].get_positions()
                    layer2_positions = ase_atom[layer2_indices].get_positions()
                    center_1 = np.mean(layer1_positions, axis=0)
                    center_2 = np.mean(layer2_positions, axis=0)

                    slip_x = float(round(abs(center_1[0]-center_2[0]), 2))
                    slip_y = float(round(abs(center_1[1]-center_2[1]), 2))
                    lateral_offsets.append([slip_x, slip_y])
                    interlayer_height.append(float(round(abs(center_1[2]-center_2[2]), 2)))
        return layers, lateral_offsets, interlayer_height


def Main():
    """
    Main function to compute the stacking pattern of COFs or layered materials like graphene
    if input is a file it reads directly from the file and if input is a directory
    it reads all the files in the directory and computes the stacking pattern for each file
    For directory input it also writes the output to a json file in the same directory
    """
    parser = argparse.ArgumentParser(description='Compute the stacking pattern of COFs or layered materials like graphene')
    parser.add_argument('input', type=str, help='Input file or directory')
    parser.add_argument('-o', '--output', type=str, default='cof_stacking_output.json', help='Output file for directory input (default: cof_stacking_output.json)')
    parser.add_argument('-v', '--verbose', action='store_true', help='print verbose output')
    args = parser.parse_args()
    input_path = args.input
    if os.path.isfile(input_path):
        ase_atom = read(input_path)
        result = compute_cof_stacking(ase_atom)
        if result is None:
            print(f"{input_path} is not a layered structure; "
                  "no stacking pattern could be computed.")
            return
        layers, lateral_offsets, interlayer_height = result
        print("Layers:", layers)
        print("Lateral Offsets:", lateral_offsets)
        print("Interlayer Heights:", interlayer_height)
    elif os.path.isdir(input_path):
        output_data = {}
        for file in glob.glob(os.path.join(input_path, '*.cif*')):
            try:
                print(f"Processing file: {file}")
                ase_atom = read(file)
                layers, lateral_offsets, interlayer_height = compute_cof_stacking(ase_atom)

                basename = os.path.basename(file).split('.')[0]
                output_data[basename] = {
                    "Layers": layers,
                    "Lateral Offsets": lateral_offsets,
                    "Interlayer Heights": interlayer_height
                }
                print(f"File: {file}")
                print("Layers:", layers)
                print("Lateral Offsets:", lateral_offsets)
                print("Interlayer Heights:", interlayer_height)
            except Exception as e:
                print(f"Error processing file {file}: {e}")
        output_file = os.path.join(".", args.output)
        write_json(output_data, output_file)
        print(f"Output written to {output_file}")
    else:
        print("Invalid input. Please provide a valid file or directory path.")

if __name__ == "__main__":
    Main()

