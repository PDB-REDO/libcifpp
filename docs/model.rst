Molecular Model
===============

Theoretically it is possible to get along with only the classes *cif::file*, *cif::datablock* and *cif::category*. But to keep your data complete and valid you then have to update lots of categories for all but the simplest manipulations. For this *libcifpp* comes with a higher level API modelling atoms, residues, monomers, polymers and complete structures in their respective classes.

Note that these classes only work properly if you are using *mmCIF* files and have an mmcif_pdbx dictionary available, either compiled in using `mrc <https://github.com/mhekkel/mrc.git>`_ or installed in the proper location.

.. note::

	This part of *libcifpp* is the least developed part. What is available should work but functionality should eventually be extended.

Atom
----

The :cpp:class:`cif::mm::atom` is a lightweight proxy class giving access to the data stored in *atom_site* and *atom_site_anisotrop*. It only caches the most often used item data and every modification is directly written back into the *mmCIF* categories.

Atoms can be copied by value with low cost. The atom class only contains a pointer to an implementation that is reference counted.

Residue, Monomer and Polymer
----------------------------

The :cpp:class:`cif::mm::residue`, :cpp:class:`cif::mm::monomer` and :cpp:class:`cif::mm::polymer` implement what you'd expect. A monomer is a residue that is part of a polymer and thus has a sequence number and siblings.

Sugars & Branches
-----------------

There are also classes for modelling sugars and sugar branches. You can create sugar branches

Structure
---------

The :cpp:class:`cif::mm::structure` can be used to load one of the models from an *mmCIF* file. By default the first model is loaded. (Multiple models are often only available files containing structures defined using NMR).

A structure holds a reference to a *cif::datablock* and retrieves its data from this datablock and writes any modification back into that datablock.

One of the most useful parts of the structure class is the ability to create and modify residues. This updates related *chem_comp* and *entity* categories as well.