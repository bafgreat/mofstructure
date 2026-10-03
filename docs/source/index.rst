mofstructure documentation
==========================

Introduction
------------

`mofstructure` takes a crystal structure and answers the questions that
usually follow: what it is built from, how porous it is, and what net it
forms. It works on metal-organic frameworks, covalent organic frameworks and
zeolites, reading CIF or any other format ASE handles.

.. raw:: html

   <video width="600" height="400" autoplay loop muted>
      <source src="_static/movie3.mp4" type="video/mp4">
      Your browser does not support the video tag.

   </video>

Key Features of `mofstructure`
------------------------------

The `mofstructure` module includes a variety of features that simplify common operations and enhance the workflow for users working with MOFs and similar materials. Some of the key functionalities include:

1. **Identify the net of a framework:**

   - The module names the underlying net from an archive of RCSR symbols, IZA
     zeolite framework-type codes and EPINET nets, and computes a canonical key
     that identifies the net whether or not any archive names it. The key is
     unchanged by supercell, atom order or choice of origin, so it can be
     stored as a database handle.

2. **Computation of Geometric Properties of MOFs:**

   - `mofstructure` integrates seamlessly with the `zeo++` software in the background to enable quick and accurate computation of all porosity-related properties. Users can easily obtain essential metrics such as Pore Limiting Diameter (PLD), Largest Cavity Diameter (LCD), Accessible Surface Area (ASA), and other geometric characteristics critical to the analysis of MOFs.

3. **Automated Removal of Unbound Guest Molecules:**

   - The module offers an automated process for identifying and removing unbound guest molecules from the framework. This feature is particularly useful when preparing structures for simulations or other computational analyses where the presence of unbound molecules could skew results.

4. **Deconstruction of Metal-Organic Frameworks into Building Units:**

   - `mofstructure` allows users to deconstruct MOFs into their constituent building units, including organic ligands, metal clusters, organic secondary building units (SBUs), and metal SBUs. For each building unit, the module computes important cheminformatic identifiers such as SMILES strings, InChI, and InChIKey. Additionally, it identifies the type of metal SBU and determines the coordination number of the central metal atom, which is crucial for understanding the structural properties of the framework.

5. **Determination of open metal sites (OMS) in MOFs:**

   - The module can identify and characterize open metal sites within MOF structures. This is particularly important for applications such as catalysis, where the presence of OMS can significantly influence the material's performance.

6. **Wrapping Systems Around Unit Cells to Remove the Effect of Periodic Boundary Conditions (PBC):**

   - When visualizing CIF files or converting CIF files to XYZ format, systems may appear uncoordinated due to the effects of periodic boundary conditions. `mofstructure` provides a solution by wrapping systems around their unit cells, ensuring a more accurate and visually coherent representation of the structure.

7. **Separation of Building Units into Regions:**

   - This feature is essential for users who need to substitute specific ligands or building units within a framework. By separating building units into distinct regions, `mofstructure` enables targeted modifications, allowing for precise customization of the framework's properties.

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   installation
   usage
   cluster_runs
   examples
   api_reference
   updates

Support
=======

The module does more than this guide covers. If you are stuck, or need a
quantity that is not yet available, please open an issue on GitHub or email
bafgreat@gmail.com.

Roadmap
=======

1. Substitute building units in a MOF to enable framework functionalisation
2. Topological analysis of metal-organic cages and other discrete assemblies

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`