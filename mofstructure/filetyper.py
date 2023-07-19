#!/usr/bin/python
from __future__ import print_function
__author__ = "Dr. Dinga Wonanke"
__version__ = '0.1.0'
__status__ = "production"
import sys
import subprocess
import os
import csv
import json, codecs
import numpy as np
import pickle
from ase import Atoms
import ase 

    
    
class AtomsEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, Atoms):
            coded = dict(positions=[list(pos) for pos in obj.get_positions()],lattice_vectors=[list(c) for c in obj.get_cell()],labels=list(obj.get_chemical_symbols()))
            if len(obj.get_cell()) == 3:
                coded['periodic'] = ['True', 'True', 'True']
            coded['n_atoms'] = len(list(obj.get_chemical_symbols()))
            coded['atomic_numbers'] = obj.get_atomic_numbers().tolist()
                
            keys = list(obj.info.keys())
            if 'atom_indices_mapping' in keys:
                info = obj.info
                coded.update(info)
            return coded
        if isinstance(obj, ase.spacegroup.Spacegroup):
            return obj.todict()
        return json.JSONEncoder.default(self, obj)

def Numpy_to_Json(Array, file_name):
    json.dump(Array.tolist(), codecs.open(file_name, 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=True)
    return
    
def Write_Json(list,file_name):
    json.dump(list, codecs.open(file_name, 'w', encoding='utf-8'))
    
    
def Json_to_Numpy(json_file):
    read_json = codecs.open(json_file, 'r', encoding='utf-8').read()
    read_json = np.array(json.loads(read_json))
    return read_json
    
def Write_JSON_ATOM(data,  AtomsEncoder, filename):
    with open(filename, 'w') as f:
        json.dump(data, f, indent=4, sort_keys=False, cls=AtomsEncoder)
    return

def Append_JSON_ATOM(data,  AtomsEncoder, filename):
    with open(filename, 'r+') as f:
        # First we load existing data into a dict.
        file_data = json.load(f)
        # Join new_data with file_data inside emp_details
        file_data.update(data)
        # Sets file's current position at offset.
        f.seek(0)
        # convert back to json.

        json.dump(data, f, indent=4, sort_keys=False, cls=AtomsEncoder)
  

def append_json(new_data, filename):
    with open(filename,'r+') as file:
          # First we load existing data into a dict.
        file_data = json.load(file)
        # Join new_data with file_data inside emp_details
        file_data.update(new_data)
        # Sets file's current position at offset.
        file.seek(0)
        # convert back to json.
        json.dump(file_data, file, indent = 4, sort_keys=True)
    
def Read_Json(file_name):
    f = open(file_name)
    data = json.load(f)
    f.close()
    return data
    

def CSV_reader(file):
    f = open(file, 'r')
    data = csv.reader(f)
    return data

def get_contents(filename):
    with open(filename, 'r') as f:
        contents = f.readlines()
    return contents

def put_contents(filename, output):
    with open(filename, 'w') as f:
        f.writelines(output)
    return


def append_contents(filename, output):
    with open(filename, 'a') as f:
        f.writelines(output)
    return


def Pickle_load(filename):
    data = open(filename, 'rb')
    data = pickle.load(data)
    return data
    
def Read_zip(zip_file):
    from zipfile import ZipFile
    content = ZipFile(zip_file, 'r')
    content.extractall(zip_file)
    content.close()
    return content
