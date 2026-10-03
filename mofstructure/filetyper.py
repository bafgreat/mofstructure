#!/usr/bin/python
'''
Reading and writing of the file formats used across mofstructure.

load_data dispatches on the file extension, so json, csv, pickle, xlsx and
msgpack can all be read through one call, with anything unrecognised returned
as raw text. The writers cover the same formats plus the append helpers the
command line tools use to grow a database file structure by structure.
'''
__author__ = "Dr. Dinga Wonanke"
__status__ = "production"
import os
import tempfile
import pickle
import csv
import json
import codecs
from zipfile import ZipFile
import numpy as np
from ase import Atoms
import ase
import pandas as pd
import msgpack
from importlib.resources import files

# Folder the command line tools write their database into. mofstructure
# handles MOFs, COFs and zeolites alike, so the name follows the package
# rather than one class of material.
DEFAULT_SAVE_DIR = "MOFstructureDB"
# Subfolders of the save directory. Every command must agree on these or a
# run stops adding to one database and starts building a second.
STRUCTURE_DATA = "Structure_Data"
XYZ_DB = "XYZ_DB"


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


def write_json_atomic(data, filename, encoder=None, indent=4, sort_keys=True):
    '''
    Write a json file so that it is never left half written.

    The data go to a temporary file in the same directory, which is flushed
    to disk and then renamed over the target. A rename within one file system
    is atomic, so a job killed at any moment leaves either the old file or the
    new one, never a truncated file that cannot be loaded.

    **parameters:**
        - data: json-serialisable object

        - filename: str
            Destination.

        - encoder: json.JSONEncoder subclass, optional

        - indent: int or None

        - sort_keys: bool
    '''
    directory = os.path.dirname(os.path.abspath(filename)) or "."
    os.makedirs(directory, exist_ok=True)
    handle, temporary = tempfile.mkstemp(
        prefix=".tmp_", suffix=".json", dir=directory)
    try:
        with os.fdopen(handle, "w", encoding="utf-8") as file:
            json.dump(data, file, indent=indent, sort_keys=sort_keys,
                      cls=encoder)
            file.flush()
            os.fsync(file.fileno())
        os.replace(temporary, filename)
    except BaseException:
        if os.path.exists(temporary):
            os.remove(temporary)
        raise


def _load_json_or_empty(filename):
    '''
    Load a json mapping for updating, treating a missing or empty file as {}.

    A file that exists but cannot be parsed is not silently replaced, since
    that would discard every record it holds. The error names the file so it
    can be repaired or moved aside.
    '''
    if not os.path.exists(filename) or os.path.getsize(filename) == 0:
        return {}
    try:
        with open(filename, encoding="utf-8") as file:
            return json.load(file)
    except json.JSONDecodeError as error:
        raise ValueError(
            f"{filename} is not valid json ({error}). It was probably "
            "truncated by an interrupted run with an older version of "
            "mofstructure. Move it aside and rebuild it with "
            "mofstructure_merge, or repair it by hand."
        ) from error


def append_json_atom(data, encoder, filename):
    '''
    Add records containing ASE atoms objects to a json file.

    **parameters:**
        - data: python dictionary
            New records; existing keys are overwritten.

        - encoder: json.JSONEncoder subclass able to serialise the atoms.

        - filename: str
    '''
    file_data = _load_json_or_empty(filename)
    file_data.update(data)
    write_json_atomic(file_data, filename, encoder=encoder, sort_keys=False)


def summary_frame(records, drop=()):
    '''
    Build the csv summary of a per structure result dictionary.

    A structure that failed outright is recorded as None rather than as a
    dictionary of results. Such an entry has no fields to become columns, and
    stacking it raises, so it is left out of the summary rather than allowed
    to take down the write at the end of a long run.

    **parameters:**
        - records: mapping
            Structure name to a dictionary of results, or None.

        - drop: iterable of str
            Columns to leave out, used for bulky text fields.

    **returns:**
        pandas.DataFrame
            One row per structure that produced results, indexed by name.
    '''
    rows = {name: result for name, result in records.items()
            if isinstance(result, dict) and result}
    data_f = pd.DataFrame.from_dict(rows, orient='index')
    data_f.index.name = 'mof_names'
    return data_f.drop(columns=list(drop), errors='ignore')


def append_json(new_data, filename):
    '''
    Add records to a json file, overwriting existing keys.

    The file is rewritten atomically, so an interrupted job cannot leave it
    truncated. Rewriting a large file after every structure is slow, which is
    why the batch commands record progress in json-lines files and only call
    this when they consolidate.

    **parameters:**
        - new_data: python dictionary

        - filename: str
    '''
    file_data = _load_json_or_empty(filename)
    file_data.update(new_data)
    write_json_atomic(file_data, filename)


def read_json(file_name):
    '''
    load a json file
    '''
    with open(file_name, encoding='utf-8') as f_obj:
        data = json.load(f_obj)

    return data


def csv_read(csv_file):
    '''
    Read a csv file
    '''
    f_obj = open(csv_file, encoding='utf-8')
    data = csv.reader(f_obj)
    return data


def get_contents(filename):
    '''
    Read a file and return a list content
    '''
    with open(filename, encoding='utf-8') as f_obj:
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


def save_dict_msgpack(data: dict, filename: str) -> None:
    """Save a dictionary to a file using MessagePack."""
    with open(filename, "wb") as f:
        msgpack.pack(data, f, use_bin_type=True)


def load_dict_msgpack(filename: str) -> dict:
    """Load a dictionary from a MessagePack file."""
    with open(filename, "rb") as f:
        return msgpack.unpack(f, raw=False, strict_map_key=False)


def convert_numpy_types(data):
    '''
    A function that jsonifies a dictionary by removing numpy data types.

    **parameter:**
        data (dict): data to jsonify

    **returns:**
        data (dic): jsonified data
    '''
    if isinstance(data, dict):
        return {key: convert_numpy_types(value) for key, value in data.items()}
    elif isinstance(data, list):
        return [convert_numpy_types(element) for element in data]
    elif isinstance(data, np.integer):
        return int(data)
    elif isinstance(data, np.floating):
        return float(data)
    else:
        return data


def load_data(filename):
    '''
    function that recognises file extenion and chooses the correction
    function to load the data.
    '''
    basename = os.path.basename(filename)
    file_ext = basename.split('.')[-1]

    if file_ext == 'json':
        data = read_json(filename)
    elif file_ext == 'csv':
        data = pd.read_csv(filename)
    elif file_ext == 'p':
        data = pickle_load(filename)
    elif file_ext == 'xlsx':
        data = pd.read_excel(filename)
    elif file_ext == 'msgpack':
        data = load_dict_msgpack(filename)
    else:
        data = get_contents(filename)
    return data


def load_iupac_names():
    '''
    Load the ligand name database that ships with the package.

    The keys are the identifiers produced by
    mofdeconstructor.name_lookup_keys, meaning a full InChIKey, a canonical
    SMILES and an InChIKey connectivity block for every molecule. Query it
    through mofdeconstructor.lookup_iupac_name rather than by direct indexing,
    since a ligand taken out of a framework has to be normalised the same way
    the keys were before it will match.

    **returns:**
        dictionary mapping a lookup key to an IUPAC name.
    '''
    msgpack_path = files("mofstructure").joinpath("db/iupacname_smiles.msgpack")
    return load_data(msgpack_path)
