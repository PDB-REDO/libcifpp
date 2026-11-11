CQL
===

The structure of CIF files (even of STAR files, from which this format is derived) looks suspiciously like a relational database. When you consider categories to be tables and items to be columns you're almost there. The only problem is linking tables, in common CIF files this is done based on multiple columns and the rules are a bit fuzzy allowing for empty columns to still match related columns that do have a value.

An early version of the tool *mmCQL* contained a SQL-like language interpreter to SELECT and UPDATE values in CIF files. This functionality has been expanded by implementing a full SQL interface using the `SQLite <https://sqlite.org>`_ library. libcifpp categories are exposed as virtual tables in a SQLite environment and can be queried and manipulated using SQL syntax.

Creating or altering categories through SQL is not supported yet. The virtual table schema is derived from the datablock and its dictionary, so statements like CREATE TABLE and ALTER TABLE would have to be intercepted by a preprocessor. That is on the to-do list.

The *mmcql* tools in `cif-tools <https://github.com/PDB-REDO/cif-tools>`_ use this backend and form a command line application similar to the *sqlite* or *psql* tools for regular SQLite files and PostgreSQL databases respectively.

Synopsis
--------

.. literalinclude:: ../examples/example-cql.cpp
	:language: c++
	:start-at: #include <cif++/cif++.hpp>

Usage
-----

To start using CQL, you will first have to create a :cpp:class:`cif::cql::connection` to a :cpp:class:`cif::datablock`. Using this connection you can create a :cpp:class:`cif::cql::transaction`. With the transaction you can execute SQL statements using :cpp:func:`cif::cql::transaction::exec`, or use :cpp:func:`cif::cql::transaction::stream` to directly pull values from the result.

The result of :cpp:func:`cif::cql::transaction::exec` is a :cpp:class:`cif::cql::result` class which uses a :cpp:class:`cif::category` as storage class.

Implementation Details
----------------------

When the datablock contains a validator (i.e. you loaded a dictionary) the SQL engine knows about all possible items/columns that are allowed. It also knows about links/relations between categories, just like the regular libcifpp query mechanism. So, updating and deleting will cascade automatically.

Another point is data types. CIF files can have numbers, strings or NULL values. Same goes for SQLite. However, when a file was loaded without a dictionary, the type of an item is dependent on its content. If something was parsed as being a number, the type will be numeric. If however, the file does contain a dictionary/validator, the type is determined by this dictionary. So, even if it looks like a number, it still might be a string internally. Good example is the ID field in atom_site, or the auth_seq_id/auth_seq_num fields. In the WHERE clause this may have unexpected results, so you may have to fall back to using `CAST <https://sqlite.org/lang_expr.html#castexpr>`_.

Each row has a stable rowid that is not reused while the category exists. You can use it in UPDATE and DELETE statements, and a rowid that no longer exists results in an error rather than affecting a different row.

The API for this functionality is a bit new, there may be room for improvement. Ideas are welcome.