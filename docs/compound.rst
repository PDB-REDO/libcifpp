Chemical Compounds
==================

The data in *CIF* and *mmCIF* files often describes the structure of some chemical compounds. The structure is recorded in the categories *atom_site* and friends. Records in these categories refer to chemical compounds using a compound ID. This compound ID is the ID field of the *chem_comp* category. For all of the known compounds in the PDB there is an entry in the Chemical Component Dictionary or `CCD <https://www.wwpdb.org/data/ccd>`_. If *libcifpp* was properly installed you have a copy of this file somewhere on your disk. And if you have installed the update scripts, a fresh version of this file will be retrieved weekly.

Information about compounds is captured in the :cpp:class:`cif::compound`. An instance of a compound object for a certain compound ID can be obtained by using the singleton :cpp:class:`cif::compound_factory`.

If the compound you want to use is not available in the CCD, you can add that information yourself. For this you can use the method :cpp:func:`cif::compound_factory::push_dictionary`. There is also :cpp:class:`cif::compound_source`, a stack based helper that adds a dictionary for the lifetime of the object.

So, given that we have the CCD and user defined compound definitions, what will you get when you try to retrieve such a compound by ID? The factory has a stack of compound generators. The base of the stack is the CCD file (*components.cif*), if it can be found. Dictionaries you add with *push_dictionary* or *compound_source* are pushed on top of that stack. The generators are searched from the top of the stack down, so the most recently added dictionary is searched first and the CCD is searched last. If no compound was created at all, nullptr is returned.

Updating CCD
------------

The CCD data is stored in a single file called *components.cif* and can be downloaded from `CCD <https://www.wwpdb.org/data/ccd>`_.

As can be read in the section on resources (:doc:`/resources`) files in libcifpp are loaded in a specific order. If the CCD datafile was downloaded during installation, a copy can be found in the directory */usr/share/libcifpp/* (if you installed in */usr*). This is a static file and will not be updated until the next installation of libcifpp.

When configuring libcifpp, the *CIFPP_INSTALL_UPDATE_SCRIPT* option controls whether an update script is installed. It is enabled by default on Linux and FreeBSD and can be disabled with:

.. code-block:: console

	cmake -S . -B build -DCIFPP_INSTALL_UPDATE_SCRIPT=OFF # ... more options?

When enabled, a script named *update-libcifpp-data* is installed in the *cron.weekly* directory below *CMAKE_INSTALL_SYSCONFDIR*. For a system-wide installation with prefix */usr* this is */etc/cron.weekly*; for a local installation it will be e.g. *$HOME/.local/etc/cron.weekly*, which your system may not run automatically.

The script reads a config file named *libcifpp.conf* from *CMAKE_INSTALL_SYSCONFDIR* (usually */etc/libcifpp.conf*). In this file the following line needs to be uncommented:

.. code-block:: console

	update=true

After that, the update script will weekly download the latest *components.cif* file to the directory configured by *CIFPP_CACHE_DIR*. The default is */var/cache/libcifpp* for a */usr/local* prefix and *CMAKE_INSTALL_FULL_LOCALSTATEDIR/cache/libcifpp* otherwise.
