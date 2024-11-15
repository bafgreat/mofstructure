#!/usr/bin/python
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"
import numpy as np
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
                    center_1 =  np.mean(layer1_positions, axis=0)
                    center_2 =  np.mean(layer2_positions, axis=0)

                    slip_x = round(abs(center_1[0]-center_2[0]), 2)
                    slip_y = round(abs(center_1[1]-center_2[1]), 2)
                    lateral_offsets.append([slip_x, slip_y])
                    interlayer_height.append(round(abs(center_1[2]-center_2[2]), 2))
        return layers, lateral_offsets, interlayer_height
