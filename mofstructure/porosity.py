#!/usr/bin/python
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke"
__version__ = '0.1.0'
__status__ = "production"
import os
import re
import numpy as np
from pyzeo.netstorage import AtomNetwork
from pyzeo.area_volume import volume, surface_area
from pymatgen.io.ase import AseAtomsAdaptor
import mofstructure.filetyper as File_typer



def zeo_calculation(ase_atom, probe_radius=1.86, number_of_steps=5000):
    '''
    Main script to compute geometric structure of porous systems.
    The focus here is on MOF, but the script can run on any periodic
    system.
    The script computes the accesible surface area, accessible volume
    and the pore geometry.There are many more outputs which can be extracted
    from ,vol_str and sa_str.
    More there are also other computation that can be done. Check out the
    test directory in dependencies/pyzeo/test. Else contact bafgreat@gmail.com
    if you need more output and can't figure it out.
    Main parameter:
        ase atom object
    return
        python dictionary containing
        1) Accessible volume fraction
        2) Accessible volume (A^3)
        3) Accessible surface area (A^2)
        4) List of the surface area along identified pockets (A^2)
        5) Largest cavity diameter:The largest sphere that can
        be inserted system without overlapping with any of the atoms
        6) Free sphere
        7) Pore limiting diameter: The largest sphere that can freely
        diffuse through the porous network without overlapping with
        any of the atoms
    '''
    tmp_cssr = 'tmp.cssr'
    tmp_out = 'tmp.res'
    tmp = ase_to_zeoobject(ase_atom)
    File_typer.put_contents(tmp_cssr, tmp)
    parameters = {}
    atmnet = AtomNetwork.read_from_CSSR(tmp_cssr)
    vol_str = volume(
        atmnet, probe_radius, probe_radius, number_of_steps, high_accuracy=False)
    vol_str = vol_str.decode("utf-8").split()
    parameters['AV_Volume_fraction'] = np.float64(vol_str[10])
    parameters['AV'] = np.float64(vol_str[8])
    sa_str = surface_area(atmnet, probe_radius, probe_radius, number_of_steps, high_accuracy=False)
    sa_str = sa_str.decode("utf-8").split()
    parameters['ASA'] = np.float64(sa_str[8])
    parameters['Number_of_channels'] = np.int64(sa_str[20])
    atmnet.calculate_free_sphere_parameters(tmp_out)
    outlines = File_typer.get_contents(tmp_out)
    data = outlines[0].split()
    parameters['LCD'] = np.float64(data[1])
    parameters['lfpd'] = np.float64(data[3])
    parameters['PLD'] = np.float64(data[2])
    if os.path.exists(tmp_cssr):
        os.remove(tmp_cssr)
    if os.path.exists(tmp_out):
        os.remove(tmp_out)
    return parameters


def ase_to_zeoobject(ase_atom):
    '''
    Converts an ase atom type to a zeo++ Cssr object
    In zeo++ the xyz coordinate system is rotated to a zyx format.

    Args:
        ase atom
    Return:
        zeo++ cssr object
    '''
    pymol = AseAtomsAdaptor.get_structure(ase_atom)
    a, b, c = ase_atom.cell.lengths()
    alpha, beta, gama = ase_atom.cell.angles()
    load = [
        f"{c:.4f} {b:.4f} {a:.4f}",
        f"{gama:.2f} {beta:.2f} {alpha:.2f} SPGR =  1 P 1    OPT = 1",
        f"{len(ase_atom)} 0",
        f"{pymol.formula}"
    ]
    for index, atom in enumerate(ase_atom):
        charge = pymol[index].charge if hasattr(pymol[index], "charge") else 0
        element = atom.symbol
        position = ase_atom.get_scaled_positions()[index]
        load.append(
            f"{index+1} {element} { position[2]:.4f} {position[1]:.4f}  {position[0]:.4f} 0 0 0 0 0 0 0 0 {charge:.4f}")
    return "\n".join(load)
