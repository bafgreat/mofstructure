#!/usr/bin/python
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"
import os
import pickle
import csv
import json
import codecs
from zipfile import ZipFile
import numpy as np
from ase import Atoms
import ase
import pandas as pd


class AtomsEncoder(json.JSONEncoder):
    '''
    ASE atom type encorder for json to enable serialising
    ase atom object.
    '''

    def default(self, encorder_obj):
        '''
        define different encoder to serialise ase atom objects
        '''
        if isinstance(encorder_obj, Atoms):
            coded = dict(positions=[list(pos) for pos in encorder_obj.get_positions()], lattice_vectors=[
                         list(c) for c in encorder_obj.get_cell()], labels=list(encorder_obj.get_chemical_symbols()))
            if len(encorder_obj.get_cell()) == 3:
                coded['periodic'] = ['True', 'True', 'True']
            coded['n_atoms'] = len(list(encorder_obj.get_chemical_symbols()))
            coded['atomic_numbers'] = encorder_obj.get_atomic_numbers().tolist()
            keys = list(encorder_obj.info.keys())
            if 'atom_indices_mapping' in keys:
                info = encorder_obj.info
                coded.update(info)
            return coded
        if isinstance(encorder_obj, ase.spacegroup.Spacegroup):
            return encorder_obj.todict()
        return json.JSONEncoder.default(self, encorder_obj)


def numpy_to_json(ndarray, file_name):
    '''
    Serialise a numpy object
    '''
    json.dump(ndarray.tolist(), codecs.open(file_name, 'w',
              encoding='utf-8'), separators=(',', ':'), sort_keys=True)
    return


def write_json(json_obj, file_name):
    '''
    write a python dictionary object to json
    '''
    # Serializing json
    json_object = json.dumps(json_obj, indent=4, sort_keys=True)
    with open(file_name, "w", encoding='utf-8') as outfile:
        outfile.write(json_object)

def json_to_numpy(json_file):
    '''
    serialised a numpy array to json
    '''
    json_reader = codecs.open(json_file, 'r', encoding='utf-8').read()
    json_reader = np.array(json.loads(json_reader))
    return read_json


def json_to_ase_atom(data,  encoder, filename):
    '''
    serialise an ase atom type and write as json
    '''
    with open(filename, 'w', encoding='utf-8') as f_obj:
        json.dump(data, f_obj, indent=4, sort_keys=False, cls=encoder)
    return

def append_json_atom(data,  encoder, filename):
    '''
    append a data containing an ase atom object
    '''
    if not os.path.exists(filename):
        with open(filename, 'w', encoding='utf-8') as f_obj:
            f_obj.write('{}')
    elif os.path.getsize(filename) == 0:
        with open(filename, 'w', encoding='utf-8') as f_obj:
            f_obj.write('{}')
    with open(filename, 'r+', encoding='utf-8') as f_obj:
        # First we load existing data into a dict.
        file_data = json.load(f_obj)
        # Join new_data with file_data inside emp_details
        file_data.update(data)
        # Sets file's current position at offset.
        f_obj.seek(0)
        # convert back to json.

        json.dump(data, f_obj, indent=4, sort_keys=False, cls=encoder)



def append_json(new_data, filename):
    '''
    append a new data in an existing json file
    '''
    if not os.path.exists(filename):
        with open(filename, 'w', encoding='utf-8') as file:
            file.write('{}')
    elif os.path.getsize(filename) == 0:
        with open(filename, 'w', encoding='utf-8') as file:
            file.write('{}')
    with open(filename, 'r+', encoding='utf-8') as file:
        # First we load existing data into a dict.
        file_data = json.load(file)
        # Overwrite existing keys with new_data
        file_data.update(new_data)
        # Sets file's current position at offset.
        file.seek(0)
        # convert back to json.
        json.dump(file_data, file, indent=4, sort_keys=True)


def read_json(file_name):
    '''
    load a json file
    '''
    with open(file_name, 'r', encoding='utf-8') as f_obj:
        data = json.load(f_obj)

    return data


def csv_read(csv_file):
    '''
    Read a csv file
    '''
    f_obj = open(csv_file, 'r', encoding='utf-8')
    data = csv.reader(f_obj)
    return data


def get_contents(filename):
    '''
    Read a file and return a list content
    '''
    with open(filename, 'r', encoding='utf-8') as f_obj:
        contents = f_obj.readlines()
    return contents


def put_contents(filename, output):
    '''
    write a list object into a file
    '''
    with open(filename, 'w', encoding='utf-8') as f_obj:
        f_obj.writelines(output)
    return


def append_contents(filename, output):
    '''
    append contents into a file
    '''
    with open(filename, 'a', encoding='utf-8') as f_obj:
        f_obj.writelines(output)
    return


def pickle_load(filename):
    '''
    load a pickle file
    '''
    data = open(filename, 'rb')
    data = pickle.load(data)
    return data


def read_zip(zip_file):
    '''
    read a zip file
    '''
    content = ZipFile(zip_file, 'r')
    content.extractall(zip_file)
    content.close()
    return content


def load_data(filename):
    '''
    function that recognises file extenion and chooses the correction
    function to load the data.
    '''
    file_ext = filename[filename.rindex('.')+1:]
    if file_ext == 'json':
        data = read_json(filename)
    elif file_ext == 'csv':
        data = pd.read_csv(filename)
    elif file_ext == 'p':
        data = pickle_load(filename)
    elif file_ext == 'xlsx':
        data = pd.read_excel(filename)
    else:
        data = get_contents(filename)
    return data
