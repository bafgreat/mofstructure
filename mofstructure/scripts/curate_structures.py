#!/usr/bin/python
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"
import os
import argparse
from mofstructure import structure
import mofstructure.mofdeconstructor as MOF_deconstructor
import mofstructure.filetyper as read_write



def remove_guest(ase_atom):
    '''
    Simple function to remove guest molecules in porous system
    '''
    index_non_guest = MOF_deconstructor.remove_unbound_guest(ase_atom)
    return ase_atom[index_non_guest]



def curate_data(cif_files, result_folder, verbose=False):
    '''
    A workflow to quickly curate a crystal structure. At the moment the
    main fucntionality is the removal of unbound guest molecules. The
    code takes a file, reades the coordinates and removes any unbound
    guest molecule.

    The function by default creates a new folder called Curated_data
    and write the curated data in exactly the same ase format as the input.

    The function can act both on a single file or list of files found in
    as folder/directory.

    **Use Case**
    on the command line simply run

    mofstructure_curate cif_file
    or
    mofstructure_curate cif_folder

    **Out look**
    In the future the code should correct for missing hydrogens and
    curate systems with overlapping atoms.

    Parameters
    ----------
    cif_file : a cif file or any ase readable file containing a MOF.
    result_folder : path to output folder
    '''

    if not os.path.exists(result_folder):
        os.makedirs(result_folder)
    path_2_general = 'general_information.json'
    if os.path.exists(path_2_general):
        general_info = read_write.load_data(path_2_general)
    else:
        general_info = {}

    if  isinstance(cif_files, str) and os.path.isfile(cif_files):
        basename = os.path.basename(cif_files)
        general_info[basename] = {}
        mof_object = structure.MOFstructure(filename=cif_files)
        overlap = MOF_deconstructor.inter_atomic_distance_check(mof_object.ase_atoms)
        general_info[basename]['has_overlapping_atoms'] = not overlap
        ase_atom = mof_object.remove_guest()
        ase_atom.write(f'{result_folder}/{basename}')
    elif isinstance(cif_files, list):
        try:
            seen_data = [f for f in os.listdir(result_folder) if f.endswith('.cif')]
            for cif_file in cif_files:
                basename = os.path.basename(cif_file)
                general_info[basename] = {}
                if basename not in seen_data:
                    print(f'Processing {basename}')
                    mof_object = structure.MOFstructure(filename=cif_file)
                    overlap = MOF_deconstructor.inter_atomic_distance_check(mof_object.ase_atoms)
                    general_info[basename]['has_overlapping_atoms'] = not overlap
                    ase_atom = mof_object.remove_guest()

                    ase_atom.write(f'{result_folder}/{basename}')
        except Exception as e:
            print(f"Error occurred when processing {cif_file}: {e}")
            pass

    read_write.append_json(general_info, path_2_general)
    if verbose:
        print(f"Saved results to {result_folder}")
    return


def main():
    '''
    CLI for cleaning crytstal structures.
    '''
    parser = argparse.ArgumentParser(
        description='Run work_flow function with optional verbose output')
    parser.add_argument('cif_folder', type=str,
                        help='single file or a directory. ')

    parser.add_argument('-s', '--save_dir', type=str,
                        default='Curated_data', help='directory to save output files')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='print verbose output')
    args = parser.parse_args()
    if os.path.isfile(args.cif_folder):
        cif_files = args.cif_folder
    elif os.path.isdir(args.cif_folder):
        cif_files = [os.path.join(args.cif_folder, f)
        for f in os.listdir(args.cif_folder)
        if f.endswith('.cif') and not f.startswith('._')]

    curate_data(cif_files, args.save_dir, args.verbose)
