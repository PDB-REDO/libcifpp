PDB files
=========

Although *libcifpp* uses *mmCIF* internally, it can read legacy PDB files and write *mmCIF* datablocks back to PDB format. Compressed files (gzip, xz/lzma) are detected automatically.

Reading
-------

Use :cpp:func:`cif::pdb::read` to read a file. It detects whether the content is PDB or mmCIF and returns a :cpp:class:`cif::file`:

.. code-block:: cpp

    cif::file f = cif::pdb::read("1cbs.pdb.gz");
    auto &db = f.front();

    // From here on the data is used like any other mmCIF datablock
    auto &atom_site = db["atom_site"];

The overload taking a ``std::istream`` works the same way and can be used for data that is already in memory.

Writing
-------

Use :cpp:func:`cif::pdb::write` to write a datablock. The path overload chooses the output format based on the filename extension; the stream overload writes legacy PDB:

.. code-block:: cpp

    cif::pdb::write("out.pdb", db);
    cif::pdb::write("out.cif", db);

    // Or write PDB data to an existing stream
    cif::pdb::write(std::cout, db);

The conversion is best effort: the legacy format cannot represent all *mmCIF* data, so a round trip may lose information.

Fixing and reconstructing
-------------------------

PDB files often lack categories that are required for a valid *mmCIF* file. Two helpers are available:

* :cpp:func:`cif::pdb::fixup_pdbx` performs a quick fix for the most common missing categories.
* :cpp:func:`cif::pdb::reconstruct_pdbx` does a more thorough reconstruction, including entities, sequences and polymer information.

.. code-block:: cpp

    cif::file f = cif::pdb::read("input.pdb");

    if (not cif::pdb::is_valid_pdbx_file(f))
        cif::pdb::reconstruct_pdbx(f);

    cif::pdb::write("output.cif", f.front());
