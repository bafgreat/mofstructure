#!/usr/bin/python
from __future__ import print_function

__author__ = "Dr. Dinga Wonanke"
__status__ = "production"
import os
import numpy as np
from pyzeo.netstorage import AtomNetwork
from pyzeo.area_volume import volume, surface_area
from pymatgen.io.ase import AseAtomsAdaptor
import mofstructure.filetyper as File_typer


def zeo_calculation(ase_atom, probe_radius=1.86, number_of_steps=10000, high_accuracy=True, rad_file=None):
    '''
    Main script to compute geometric structure of porous systems.
    The focus here is on MOF, but the script can run on any porous periodic
    system. The script computes the accesible surface area, accessible volume
    and the pore geometry. There are many more outputs which can be extracted
    from ,vol_str and sa_str. Moreover there are also other computation that can be done.
    Check out the test directory in dependencies/pyzeo/test. Else contact bafgreat@gmail.com
    if you need more output and can't figure it out.
    Main parameter:
        ase_atom: ase atom object
        probe_radius: The radius of the probe. Here 1.86 is used as default
        number_of_steps: Number of GCMC simulation cycles
        high_accuracy: key to determine where to perform high accuracy computation

    return
        python dictionary containing
        1) AV_Volume_fraction: Accessible volume void fraction
        2) AV_A^3: Accessible volume in A^2
        3) AV_cm^3/g: Accessible volume in cm^3/g. This values is often infinity because
                      it literatly divides the given value by the avogadro's number
        4) ASA_A^2: Accessible surface area A^2
        5) ASA_m^2/cm^3: Accessible surface area in m^2/cm^3
        6) Number_of_channels: Number of channels present in the porous system, which correspond to the number of
                               pores within the system
        7) LCD_A: The largest cavity diameter is the largest sphere that can be inserted in a porous
                  system without overlapping with any of the atoms in the system.
        8) lfpd_A:The largest included sphere along free sphere path is
                  largest sphere that can be inserted in the pore
        9)PLD_A:The pore limiting diameter is the largest sphere that can freely
                 diffuse through the porous network without overlapping with any
                 of the atoms in the system
    '''
    tmp_cssr = 'tmp.cssr'
    tmp_out = 'tmp.res'
    tmp = ase_to_zeoobject(ase_atom)
    File_typer.put_contents(tmp_cssr, tmp)
    parameters = {}
    atmnet = AtomNetwork.read_from_CSSR(tmp_cssr)
    if rad_file is not None:
        try:
            atmnet = AtomNetwork.read_from_CSSR(tmp_cssr, rad_file)
        except Exception as e:
            print ("please edit your rad file. In the meantime, default radii will be used.")
            atmnet = AtomNetwork.read_from_CSSR(tmp_cssr)

    vol_str = volume(
        atmnet, probe_radius, probe_radius, number_of_steps, high_accuracy=high_accuracy)
    if high_accuracy is True:
        vol_str = vol_str[0].decode("utf-8").split()
    else:
        vol_str = vol_str.decode("utf-8").split()
    parameters['AV_Volume_fraction'] = np.float64(vol_str[10])
    parameters['AV_A^3'] = np.float64(vol_str[8])
    # parameters['AV_cm^3/g'] = np.float64(vol_str[12])
    sa_str = surface_area(atmnet, probe_radius, probe_radius,
                          number_of_steps, high_accuracy=high_accuracy)
    if high_accuracy is True:
        sa_str = sa_str[0].decode("utf-8").split()
    else:
        sa_str = sa_str.decode("utf-8").split()
    parameters['ASA_A^2'] = np.float64(sa_str[8])
    parameters['ASA_m^2/cm^3'] = np.float64(sa_str[10])
    parameters['Number_of_channels'] = np.int64(sa_str[20])
    atmnet.calculate_free_sphere_parameters(tmp_out)
    outlines = File_typer.get_contents(tmp_out)
    data = outlines[0].split()
    parameters['LCD_A'] = np.float64(data[1])
    parameters['lfpd_A'] = np.float64(data[3])
    parameters['PLD_A'] = np.float64(data[2])
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
    a_axis, b_axis, c_axis = ase_atom.cell.lengths()
    alpha, beta, gama = ase_atom.cell.angles()
    load = [
        f"{c_axis:.4f} {b_axis:.4f} {a_axis:.4f}",
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