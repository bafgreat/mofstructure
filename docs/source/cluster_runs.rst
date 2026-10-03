Large runs and clusters
=======================

``mofstructure_database``, ``mofstructure_topology`` and
``mofstructure_porosity`` are built to run over tens of thousands of
structures without stopping. This page describes how they do that, how to
split a run over a cluster array and what the output looks like.

What protects a long run
------------------------

* **One worker process per structure.** Each structure is analysed in a
  worker process. If compiled code crashes (a segmentation fault or abort in
  Open Babel, zeo++ or numpy) only the worker dies. The structure is recorded
  as ``crashed`` and a new worker carries on with the next structure.
* **A hard time limit.** ``--max-time`` is a wall-clock limit for the whole
  analysis of one structure, including reading, guest removal,
  deconstruction, the canonical form and the embedding. A structure that runs
  past it is killed and recorded as ``timeout``.
* **An optional memory limit.** ``--memory-limit`` caps the memory of each
  worker in GB (Linux). A structure that needs more fails on its own instead
  of the scheduler killing the whole job.
* **Results are saved as they finish.** Each result is appended as one line
  to a checkpoint file and flushed to disk. A job killed by the scheduler
  loses at most the structures in progress.
* **Runs resume.** Repeating the same command skips every structure already
  recorded. ``--retry-failed`` runs the failures again.
* **Database files are never left half written.** They are written to a
  temporary file and renamed into place.

Options shared by the three commands
------------------------------------

========================  ====================================================
``-j, --workers N``        worker processes (default 1)
``--max-time S``           seconds per structure before it is killed
                           (database 7200, topology 1800, porosity 2400;
                           0 for no limit)
``--memory-limit GB``      memory per worker, Linux only
``--shard i/N``            process part *i* of *N* (0-based); ``slurm/N``
                           reads ``SLURM_ARRAY_TASK_ID``
``--retry-failed``         run again structures recorded as error, timeout
                           or crashed
``--recycle N``            structures per worker before it is replaced
                           (default 50), which releases memory
``--merge-every MIN``      rebuild the json and csv files during the run
                           every MIN minutes (default 60, 0 only at the end;
                           not used with ``--shard``)
========================  ====================================================

On one node
-----------

.. code-block:: bash

   mofstructure_database cif_folder -s MOFdb -j 32 --oms --memory-limit 8

If the job is killed, submit the same command again.

On a SLURM array
----------------

Each array task processes one shard and writes its own checkpoint part.
Join the parts once all tasks have finished.

.. code-block:: bash

   #!/bin/bash
   #SBATCH --job-name=mofstructure
   #SBATCH --array=0-31
   #SBATCH --cpus-per-task=16
   #SBATCH --mem=64G
   #SBATCH --time=24:00:00
   #SBATCH --output=logs/mofstructure_%A_%a.log

   mofstructure_database cif_folder -s MOFdb \
       -j ${SLURM_CPUS_PER_TASK} --memory-limit 3.5 \
       --shard slurm/32 --oms

.. code-block:: bash

   # after the array has finished, or at any time to inspect progress
   mofstructure_merge MOFdb

A task that hit its time limit is resubmitted with the same array index and
continues from its checkpoint. All shards must use the same ``N``.

Topology, porosity and open metal sites as separate jobs
--------------------------------------------------------

The three analyses can be submitted as independent jobs on the same folder
and the same save directory. Each writes its own checkpoints, so they do not
interfere, and each rebuilds only the files it produces.

.. code-block:: bash

   # topology.sh
   #SBATCH --cpus-per-task=32 --mem=128G --time=48:00:00
   mofstructure_topology cif_folder -s MOFdb --method all_node \
       -j ${SLURM_CPUS_PER_TASK} --timeout 600 --max-time 1800 \
       --memory-limit 3.5 --quiet

   # porosity.sh
   #SBATCH --cpus-per-task=32 --mem=128G --time=48:00:00
   mofstructure_porosity cif_folder -s MOFdb \
       -j ${SLURM_CPUS_PER_TASK} --timeout 1800 --memory-limit 3.5

   # oms.sh
   #SBATCH --cpus-per-task=32 --mem=128G --time=24:00:00
   mofstructure_oms cif_folder -s MOFdb \
       -j ${SLURM_CPUS_PER_TASK} --memory-limit 3.5

When a job reaches its time limit, submit the same script again and it
continues. When all three have finished, ``mofstructure_merge MOFdb`` writes
the final files and a ``run_status.csv`` covering all three.

``mofstructure_topology`` names a structure by its file name without the
last suffix (``ABAVIJ.MOF_subset`` for ``ABAVIJ.MOF_subset.cif``), while the
porosity, OMS and database commands keep the text before the first dot
(``ABAVIJ``). Normalise the names before joining the tables.

Output
------

.. code-block:: text

   MOFdb/
     Structure_Data/
       _progress/                     checkpoints, one line per structure
         database[.part-iiii-of-nnnn].jsonl
         topology.<method>[.part-...].jsonl
         porosity[.part-...].jsonl
       sbus_and_linkers.json/.csv
       ligands_data.json/.csv
       porosity_data.json/.csv
       structure_oms_and_general_info.json/.csv
       topology_data.json/.csv
       fingerprint_data.json/.csv
       run_status.csv
     XYZ_DB/                          building units as xyz files

The checkpoints grow by one line per structure throughout the run. The json
and csv files have the same names and layout as in earlier versions; a json
file is one object that has to be rewritten whole, so they are rebuilt from
the checkpoints every ``--merge-every`` minutes, at the end of a run that is
not sharded, and by ``mofstructure_merge``, which can be run at any time,
also while a run is in progress. Records already in the json files
from an earlier version are kept.

``run_status.csv`` has one row per structure and command:

============  =============================================================
``status``    ``ok``; ``error`` (the analysis raised, see ``detail``);
              ``timeout`` (killed after ``--max-time``); ``crashed``
              (the worker died, ``detail`` names the signal)
``elapsed_s`` wall-clock seconds for the structure
``analysis_errors``  for ``mofstructure_database``, the analyses that
              failed within a structure that otherwise completed, for
              example ``{"topology": "..."}``
============  =============================================================

A structure that timed out or crashed has a row in ``porosity_data.json``
with ``porosity_status`` set to ``timeout`` or ``crashed`` and missing values
elsewhere, so the tables keep one row per input structure.

Each checkpoint line is a json object:

.. code-block:: json

   {"name": "ABAVIJ", "file": "cif_folder/ABAVIJ.cif", "status": "ok",
    "elapsed_s": 1.83, "data": {"topology": {...}, "fingerprint": {...},
    "sbu": {...}, "ligand": {...}, "porosity": {...}}}

Debugging one structure
-----------------------

``-j 0`` runs everything in the calling process without isolation or time
limit, so a debugger or print statements work as usual.
