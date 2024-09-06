.. mofstructure documentation master file, created by
   sphinx-quickstart on Sat Aug 31 17:58:25 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

mofstructure documentation
========================================
Introduction
============

`mofstructure` is a powerful and user-friendly Python module designed for the manipulation and analysis of metal-organic frameworks (MOFs) and other porous materials, including covalent organic frameworks (COFs) and zeolites. The module offers a range of functions that streamline complex tasks related to the study and modification of these structures. Whether you are a researcher, scientist, or student working in the field of materials science, `mofstructure` provides the tools you need to efficiently analyze and manipulate these intricate frameworks.

Key Features of `mofstructure`
------------------------------

The `mofstructure` module includes a variety of features that simplify common operations and enhance the workflow for users working with MOFs and similar materials. Some of the key functionalities include:

1. **Computation of Geometric Properties of MOFs:**
   - `mofstructure` integrates seamlessly with the `zeo++` software in the background to enable quick and accurate computation of all porosity-related properties. Users can easily obtain essential metrics such as Pore Limiting Diameter (PLD), Largest Cavity Diameter (LCD), Accessible Surface Area (ASA), and other geometric characteristics critical to the analysis of MOFs.

2. **Automated Removal of Unbound Guest Molecules:**
   - The module offers an automated process for identifying and removing unbound guest molecules from the framework. This feature is particularly useful when preparing structures for simulations or other computational analyses where the presence of unbound molecules could skew results.

3. **Deconstruction of Metal-Organic Frameworks into Building Units:**
   - `mofstructure` allows users to deconstruct MOFs into their constituent building units, including organic ligands, metal clusters, organic secondary building units (SBUs), and metal SBUs. For each building unit, the module computes important cheminformatic identifiers such as SMILES strings, InChI, and InChIKey. Additionally, it identifies the type of metal SBU and determines the coordination number of the central metal atom, which is crucial for understanding the structural properties of the framework.

4. **Wrapping Systems Around Unit Cells to Remove the Effect of Periodic Boundary Conditions (PBC):**
   - When visualizing CIF files or converting CIF files to XYZ format, systems may appear uncoordinated due to the effects of periodic boundary conditions. `mofstructure` provides a solution by wrapping systems around their unit cells, ensuring a more accurate and visually coherent representation of the structure.

5. **Separation of Building Units into Regions:**
   - This feature is essential for users who need to substitute specific ligands or building units within a framework. By separating building units into distinct regions, `mofstructure` enables targeted modifications, allowing for precise customization of the framework's properties.

.. image:: images/Rotation.gif
   :alt: Generalities

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   installation
   usage
   examples
   api_reference
   updates

Support
=======

The module contains much more functionalities. If you are struggling or wish to compute a new quantity that is not yet present, feel free to send me an email: bafgreat@gamil.com.

Tutorial
=========
.. raw:: html

   <iframe src="./doc/how-to-doc.html" width="100%" height="600">
     <p>Click <a href="./doc/how-to-doc.html">here</a> for a how-to tutorial.</p>
   </iframe>

Roadmap
=======

In the future, the code should be able to:

1. Compute RCSR topological code
2. Substitute building units in a MOF to enable framework functionalization
3. Automatically curate CIFs
4. Deconstruct COFs into their building units

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`