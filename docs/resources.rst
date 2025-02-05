Resources
=========

Programs using libcifpp often need access to common data files. E.g. CIF dictionary files, CCP4 monomer restraints files or the CCD data file. In libcifpp these files are called resources. These files are often also based on external sources that are updated on a regular basis.

Resources can be compiled into the executable so that the resulting
application can be made portable to other machines. For this you
need to use `mrc <https://github.com/mhekkel/mrc.git>`_ which only works
on Un*x like systems using the ELF executable format or on MS Windows

But resources may also be located as files on the filesytem at
specific locations. And you can specify your own location for
files (a directory) or even override named resources with your
own data.

Loading Resources
-----------------

No matter where the resource is located, you should always use the single libcifpp API call :cpp:func:`cif::load_resource` to load them. This function returns a *std::istream* wrapped inside a *std::unique_ptr*. 

The order in which resources are searched for is:

* Use the resource that was defined by calling :cpp:func:`cif::add_file_resource`
  for this name.

* Search the paths specified by :cpp:func:`cif::add_data_directory`, last one
  added is searched first

* Search the so-called *CACHE_DIR*. This location is defined
  at compile time and based on the installation directory of
  libcifpp. Usually it is */var/cache/libcifpp*.
  It is in this directory where the cron job for libcifpp will
  put the updated files weekly.

* If the *CCP4* environment is available, the
  *$ENV{CCP4}/share/libcifpp* is searched.

* If the environment variable *LIBCIFPP_DATA_DIR* is set it
  is searched

* The *DATA_DIR* is searched, this is also a variable defined
  at compile time, also based on the installation directory
  of libcifpp. It usually is */usr/share/libcifpp*

* As a last resort an attempt is made to load the data from
  resources compiled by `mrc <https://github.com/mhekkel/mrc.git>`_.

