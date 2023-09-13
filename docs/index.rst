Introduction
============

Information on 3D structures of proteins originally came formatted in `PDB <http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html>`_ files. Although the specification for this format had some real restrictions like a mandatory HEADER and CRYST line, many programs implemented this very poorly often writing out only ATOM records. And users became used to this.

The legacy PDB format has some severe limitations rendering it useless for all but very small protein structures. A new format called `mmCIF <https://mmcif.wwpdb.org/>`_ has been around for decades and now is the default format for the Protein Data Bank.

The software developed in the `PDB-REDO <https://pdb-redo.eu/>`_ project aims at improving 3D models based on original experimental data. For this, the tools need to be able to work with both legacy PDB and mmCIF files. A decision was made to make mmCIF leading internally in all programs and convert legacy PDB directly into mmCIF before processing the data. A robust conversion had to be developed to make this possible since, as noted above, files can come with more or less information making it sometimes needed to do a sequence alignment to find out the exact residue numbers.

And so libcif++ came to life, a library to work with mmCIF files. Work on this library started early 2017 and has developed quite a bit since then. To reduce dependency on other libraries, some functionality was added that is not strictly related to reading and writing mmCIF files but may be useful nonetheless. This is mostly code that is used in 3D calculations and symmetry operations.

Design
------

The main part of the library is a set of classes that work with mmCIF files. They are:

* :cpp:class:`cif::file`
* :cpp:class:`cif::datablock`
* :cpp:class:`cif::category`

The :cpp:class:`cif::file` class encapsulates the contents of a mmCIF file. In such a file there are one or more :cpp:class:`cif::datablock` objects and each datablock contains one or more :cpp:class:`cif::category` objects.

Synopsis
--------

Using *libcifpp* is easy, if you are familiar with modern C++:

.. literalinclude:: ../README.md
	:language: c++
	:start-after: ```c++
	:end-before: ```

.. toctree::
   :maxdepth: 2
   :caption: Contents

   self
   basics.rst
   compound.rst
   model.rst
   resources.rst
   symmetry.rst
   bitsandpieces.rst
   api/library_root.rst
   genindex.rst

