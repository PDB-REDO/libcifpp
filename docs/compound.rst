Chemical Compounds
==================

The data in *CIF* and *mmCIF* files often describes the structure of some chemical compounds. The structure is recorded in the categories *atom_site* and friends. Records in these categories refer to chemical compounds using a compound ID. This compound ID is the ID field of the *chem_comp* category. For all of the known compounds in the PDB there is an entry in the Chemical Compounds Dictionary or `CCD <https://www.wwpdb.org/data/ccd>`_. If *libcifpp* was properly installed you have a copy of this file somewhere on your disk. And if you have installed the update scripts, a fresh version of this file will be retrieved weekly.

As an alternative to CCD there are the monomer library files from `CCP4 <https://www.ccp4.ac.uk/>`_. These contain somewhat different data but the overlap is good enough for usage in *libcifpp*.

Information about compounds is captured in the :cpp:class:`cif::compound`. An instance of a compound object for a certain compound ID can be obtained by using the singleton :cpp:class:`cif::compound_factory`.

If the compound you want to use is not available in the CCD or in CCP4, you can add that information yourself. For this you can use the method :cpp:func:`cif::compound_factory::push_dictionary`.

So, given that we have CCD, CCP4 monomer library and used defined compound definitions, what will you get when you try to retrieve such a compound by ID? The answer is, the factory has a stack of compound generators. The first thrown on the stack is the one for a CCD file (*components.cif*) if it can be found. Then, if the *CLIBD_MON* environmental variable is defined, a generator for monomer library files is added to the stack. And then all generators for files you added using *push_dictionary* are added in order. The generators are searched in the reverse order in which they were added to see if it creates a compound object for the ID. If no compound was created at all, nullptr is returned.

Updating CCD
------------

The CCD data is stored in a single file called *components.cif* and can be downloaded from `CCD <https://www.wwpdb.org/data/ccd>`_. 

As can be read in the section on resources (:doc:`/resources`) files in libcifpp are loaded in a specific order. If the CCD datafile was downloaded during installation, a copy can be found in the directory */usr/share/libcifpp/* (if you installed in */usr*). This is a static file and will not be updated until the next installation of libcifpp.

When configuring libcifpp, you can specify the *CIFPP_INSTALL_UPDATE_SCRIPT* option, as in:

.. code-block:: console

	cmake -S . -B build -DCIFPP_INSTALL_UPDATE_SCRIPT=ON # ... more options?

This will install a script named *update-libcifpp-data* in */etc/cron.weekly* or */etc/periodic/weekly*. This file uses a config file named */etc/libcifpp.conf* which you then need to edit. In this config file the following line needs to be uncommented:

.. code-block:: console

	# update=true

After that, the update script will weekly download the latest components.cif file to */var/cache/libcifpp*.
