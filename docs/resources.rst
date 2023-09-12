Resources
=========

Programs using libcifpp often need access to common data files. E.g. CIF dictionary files, CCP4 monomer restraints files or the CCD data file. In libcifpp these files are called resources. These files are often also based on external sources that are updated on a regular basis.

Loading Resources
-----------------

No matter where the resource is located, you should always use the single libcifpp API call :cpp:func:`cif::load_resource` to load them. This function returns a *std::istream* wrapped inside a *std::unique_ptr*. 