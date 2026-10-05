#!/usr/bin/python
'''
Geometric pore analysis through zeo++.

Accessible volume, accessible surface area and the pore diameters (LCD, PLD and
the largest free sphere along the percolation path) are computed by handing the
structure to zeo++ as a CSSR file.

zeo++ fails in two ways that a plain try/except cannot handle. It aborts the
process outright when the Voronoi decomposition fails its internal volume
check, which raises SIGABRT rather than a python exception, and on a large or
awkward cell it can run for hours without finishing. zeo_calculation therefore
runs the calculation in a child interpreter under a timeout, so neither an
abort nor a stall can take down a batch job.

Every call returns the same keys whatever happens: a structure that fails
gives the same record with None in place of each number and a
`porosity_status` saying why. A directory of structures therefore yields rows
of one shape that go straight into a DataFrame, rather than a mixture of
records and blanks that has to be repaired before it can be used.
'''

__author__ = "Dr. Dinga Wonanke"
__status__ = "production"
import math
import os
import pickle
import shutil
import signal
import subprocess
import sys
import tempfile
from pyzeo.netstorage import AtomNetwork
from pyzeo.area_volume import volume, surface_area
from pymatgen.io.ase import AseAtomsAdaptor
import mofstructure.filetyper as File_typer

# Field names are lower case, ASCII and free of the ^ and / that make a
# column awkward to address in pandas, SQL or parquet. Order is fixed so the
# columns of a table built structure by structure do not depend on insertion
# order.
POROSITY_FIELDS = (
    "av_volume_fraction",
    "av_a3",
    "asa_a2",
    "asa_m2_per_cm3",
    "number_of_channels",
    "lcd_a",
    "pld_a",
    "lfpd_a",
)

# Names used up to 0.1.9.0, kept here so an older table can be migrated.
LEGACY_FIELD_NAMES = {
    "AV_Volume_fraction": "av_volume_fraction",
    "AV_A^3": "av_a3",
    "ASA_A^2": "asa_a2",
    "ASA_m^2/cm^3": "asa_m2_per_cm3",
    "Number_of_channels": "number_of_channels",
    "LCD_A": "lcd_a",
    "PLD_A": "pld_a",
    "lfpd_A": "lfpd_a",
}

# Seconds a single structure may take before the child is killed. zeo++ needs
# seconds for an ordinary MOF; a cell large enough to run past this is one
# that would otherwise stall the batch indefinitely.
DEFAULT_TIMEOUT = 1800


def empty_porosity_record(status):
    '''
    Build the record returned when zeo++ could not analyse a structure.

    The keys match a successful record exactly, so a failure adds a row of
    missing values to a table rather than a hole in its shape.

    **parameters:**
        - status: str
            Why the calculation produced no numbers.

    **returns:**
        python dictionary
            Every field in POROSITY_FIELDS set to None, plus
            `porosity_status`.
    '''
    record = {field: None for field in POROSITY_FIELDS}
    record["porosity_status"] = status
    return record


def zeo_calculation(ase_atom, probe_radius=1.86, number_of_steps=10000,
                    high_accuracy=True, rad_file=None,
                    timeout=DEFAULT_TIMEOUT):
    '''
    Compute the pore geometry of a periodic system in a separate interpreter.

    zeo++ validates its Voronoi decomposition against the cell volume and calls
    abort() when the check fails, which raises SIGABRT rather than a python
    exception. No try/except can intercept that, so an unlucky structure takes
    the whole process down with it. A large cell can also leave zeo++ running
    far longer than the rest of the batch is worth waiting for. Running the
    calculation in a child under a timeout keeps both to that one structure.

    **parameters:**
        - ase_atom: ase.Atoms
            The periodic structure to analyse.

        - probe_radius: float
            Radius of the probe in Angstrom.

        - number_of_steps: int
            Monte Carlo samples used for the volume and area.

        - high_accuracy: bool
            Use the more expensive Voronoi decomposition.

        - rad_file: str, optional
            File of user defined atomic radii, `.rad` extension.

        - timeout: float, optional
            Seconds to allow the structure before the child is killed. None
            waits indefinitely, which risks stalling a batch.

    **returns:**
        python dictionary
            The fields of POROSITY_FIELDS and `porosity_status`. A structure
            zeo++ cannot handle returns the same keys with None in place of
            each number, so the caller always receives one row of one shape.
    '''
    workdir = tempfile.mkdtemp(prefix='mofstructure_zeo_')
    try:
        with open(os.path.join(workdir, 'input.pkl'), 'wb') as file_handle:
            pickle.dump({'ase_atom': ase_atom,
                         'probe_radius': probe_radius,
                         'number_of_steps': number_of_steps,
                         'high_accuracy': high_accuracy,
                         'rad_file': rad_file,
                         'timeout': timeout}, file_handle)

        # The child must import this copy of mofstructure. Started from the
        # workdir it would otherwise find whichever copy is installed, which
        # in a source checkout can be an older release.
        package_root = os.path.dirname(
            os.path.dirname(os.path.abspath(__file__)))
        env = dict(os.environ)
        env['PYTHONPATH'] = os.pathsep.join(
            filter(None, [package_root, env.get('PYTHONPATH')]))

        # cwd is the workdir so the tmp.cssr and tmp.res scratch files land
        # there and go away with it, even when the child aborts
        try:
            completed = subprocess.run(
                [sys.executable, '-m', 'mofstructure.porosity', workdir],
                cwd=workdir, check=False, timeout=timeout, env=env)
        except subprocess.TimeoutExpired:
            # run() has already killed the child by this point.
            print(f'zeo++ exceeded {timeout} s on this structure, '
                  f'skipping porosity')
            return empty_porosity_record('timeout')

        output_file = os.path.join(workdir, 'output.pkl')
        if completed.returncode != 0 or not os.path.exists(output_file):
            print(f'zeo++ could not analyse this structure '
                  f'(exit {completed.returncode}), skipping porosity')
            return empty_porosity_record(f'failed:{completed.returncode}')

        with open(output_file, 'rb') as file_handle:
            return pickle.load(file_handle)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def compute_zeo_parameters(ase_atom, probe_radius=1.86, number_of_steps=10000, high_accuracy=True, rad_file=None):
    '''
    Main script to compute geometric structure of porous systems.
    The focus here is on MOF, but the script can run on any porous periodic
    system. The script computes the accesible surface area, accessible volume
    and the pore geometry. There are many more outputs which can be extracted
    from ,vol_str and sa_str. Moreover there are also other computation that can be done.
    Check out the test directory in dependencies/pyzeo/test.
    Else contact bafgreat@gmail.com. if you need more output and can't figure it out.

    **parameter:**
        ase_atom: ASE atoms object representing the porous structure.
        probe_radius (float): Radius of the probe (default: 1.86).
        number_of_steps (int): Number of GCMC simulation cycles (default: 10000).
        high_accuracy (bool): If True, perform high-accuracy computations.
        rad_file: Optional file containing user defined atom radii. Must have the `.rad` extension

    **returns:**
        **python dictionary containing**
        1) av_volume_fraction: Accessible volume void fraction
        2) av_a3: Accessible volume in A^3
        3) asa_a2: Accessible surface area in A^2
        4) asa_m2_per_cm3: Accessible surface area in m^2/cm^3
        5) number_of_channels: Number of channels present in the porous system, which correspond to the number of pores within the system
        6) lcd_a: The largest cavity diameter is the largest sphere that can be inserted in a porous system without overlapping with any of the atoms in the system.
        7) pld_a: The pore limiting diameter is the largest sphere that can freely diffuse through the porous network without overlapping with any of the atoms in the system
        8) lfpd_a: The largest included sphere along free sphere path is largest sphere that can be inserted in the pore
        9) porosity_status: 'ok' when the numbers above were computed

        Values are plain python floats and ints rather than numpy scalars, so
        the record serialises to JSON without a custom encoder.
    '''
    # Scratch files go in a directory of their own. zeo++ writes them by name
    # in the working directory, so two calculations sharing a directory would
    # otherwise overwrite each other's input.
    scratch = tempfile.mkdtemp(prefix='mofstructure_zeo_run_')
    tmp_cssr = os.path.join(scratch, 'tmp.cssr')
    tmp_out = os.path.join(scratch, 'tmp.res')
    tmp = ase_to_zeoobject(ase_atom)
    parameters = {}
    try:
        File_typer.put_contents(tmp_cssr, tmp)
        if rad_file is not None:
            try:
                atmnet = AtomNetwork.read_from_CSSR(
                    tmp_cssr, rad_flag=True, rad_file=rad_file)
            except Exception:
                print("please edit your rad file. In the meantime, "
                      "default radii will be used.")
                atmnet = AtomNetwork.read_from_CSSR(tmp_cssr)
        else:
            atmnet = AtomNetwork.read_from_CSSR(tmp_cssr)

        vol_str = volume(atmnet, probe_radius, probe_radius,
                         number_of_steps, high_accuracy=high_accuracy)
        if high_accuracy is True:
            vol_str = vol_str[0].decode("utf-8").split()
        else:
            vol_str = vol_str.decode("utf-8").split()
        parameters['av_volume_fraction'] = float(vol_str[10])
        parameters['av_a3'] = float(vol_str[8])
        sa_str = surface_area(atmnet, probe_radius, probe_radius,
                              number_of_steps, high_accuracy=high_accuracy)
        if high_accuracy is True:
            sa_str = sa_str[0].decode("utf-8").split()
        else:
            sa_str = sa_str.decode("utf-8").split()
        parameters['asa_a2'] = float(sa_str[8])
        parameters['asa_m2_per_cm3'] = float(sa_str[10])
        parameters['number_of_channels'] = int(sa_str[20])
        atmnet.calculate_free_sphere_parameters(tmp_out)
        outlines = File_typer.get_contents(tmp_out)
        data = outlines[0].split()
        parameters['lcd_a'] = float(data[1])
        parameters['pld_a'] = float(data[2])
        parameters['lfpd_a'] = float(data[3])
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    # Fixed order, so the columns of a table do not depend on the caller.
    ordered = {field: parameters[field] for field in POROSITY_FIELDS}
    ordered['porosity_status'] = 'ok'
    return ordered


def _run_as_child():
    '''
    Entry point used by zeo_calculation, which runs this module as a script so
    that an abort() inside zeo++ only takes down the child interpreter.
    '''
    workdir = sys.argv[1]
    # Scratch files belong inside the directory the parent cleans up. A child
    # killed on timeout never runs its own cleanup, so anything written
    # elsewhere would be left behind on every slow structure.
    tempfile.tempdir = workdir
    with open(os.path.join(workdir, 'input.pkl'), 'rb') as file_handle:
        arguments = pickle.load(file_handle)

    # The parent enforces the timeout, but a parent that is itself killed
    # (a batch worker past its limit, a cancelled job) cannot, and zeo++
    # would then run on with nothing to stop it. SIGALRM ends the child even
    # inside compiled code. The margin lets the parent's own timeout fire
    # first, so a stall is still reported as a timeout.
    timeout = arguments.pop('timeout', None)
    if timeout and hasattr(signal, 'alarm'):
        signal.signal(signal.SIGALRM, signal.SIG_DFL)
        signal.alarm(int(math.ceil(timeout)) + 60)

    parameters = compute_zeo_parameters(**arguments)

    with open(os.path.join(workdir, 'output.pkl'), 'wb') as file_handle:
        pickle.dump(parameters, file_handle)


def ase_to_zeoobject(ase_atom):
    '''
    Converts an ase atom type to a zeo++ Cssr object
    In zeo++ the xyz coordinate system is rotated to a zyx format.

    **parameter:**
        ase_atom: ase atom object

    **returns:**
        cssr_object: string representing zeo++ cssr object
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


if __name__ == '__main__':
    _run_as_child()
