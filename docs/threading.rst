Thread safety
=============

Most of the data structures in *libcifpp* are not internally synchronised. A :cpp:class:`cif::file`, :cpp:class:`cif::datablock` or :cpp:class:`cif::category` should not be modified by multiple threads at the same time, and concurrent reads are only safe when no thread modifies the data. Use one instance per thread, or protect access with your own mutex.

Global state
------------

A few parts of the library do maintain global state:

* :cpp:class:`cif::compound_factory` is a singleton. Call :cpp:func:`cif::compound_factory::init` before the first use and before starting threads that use the factory. Do not call :cpp:func:`cif::compound_factory::clear` while another thread is still using the instance. If each thread needs its own set of compound definitions, initialise the factory with ``useThreadLocalInstanceOnly = true``.
* Resource registration through :cpp:func:`cif::add_file_resource` and :cpp:func:`cif::add_data_directory` is synchronised, but the simplest and safest pattern is still to register all resources before starting worker threads.
* :cpp:class:`cif::validator_factory` is thread-safe and returns shared, immutable validators.

CQL
---

A :cpp:class:`cif::cql::connection` and its :cpp:class:`cif::cql::transaction` objects are not thread-safe. Use one connection per thread. The underlying datablock must also not be modified through another path while a query is running.

Progress bars
-------------

The fancy progress bar is safe to update from multiple threads: its progress counter is atomic and its output is protected by an internal mutex. The simple progress bar does not start a thread and is intended for single-threaded use. Writing output from multiple progress bars at the same time will interleave on the terminal, so use a single progress bar per output stream.
