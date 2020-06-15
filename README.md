libcif++
========

This is the README file for libcif++.

Requirements
------------

The code for this library was written in C++14. You therefore need a
recent compiler to build it. For the development gcc 7.4 and clang 6.0
have been used.

Other requirements are:

- A recent CCP4 environment (with the clipper libraries that have
  electrion scattering support).
- GNU make version 4.1 or higher.
- Boost libraries, the current version was developed using version 1.65
- [mrc](https://github.com/mhekkel/mrc), a resource compiler that
  allows including data files into the executable making them easier to
  install. Strictly this is optional, but at the expense of a lot of
  functionality.
- [newuoa-cpp](https://github.com/elsid/newuoa-cpp), required to
  calculate atom radii.

Building
--------

Make sure you install the libraries and tools in list above first
before building. You don't have to install them in system locations,
paths can be set as described in the next section.

There are two makefiles, one located in the directory libcif++ and one
in tools.

Both makefiles will include a *make.config* file (which will be
generated if it doesn't exist). This configuration file can be used to
override local settings, e.g. the location of certain libraries.

Before running make, first `source` the ccp4 environment. This will
take care of setting up the make.config file correctly. You then only
have to edit the `ZEEP_DIR` variable to point to the correct
directory in case you did not install libzeep, and if you use
libzeep version 4, you will have to change `ZEEP_INCL` and `ZEEP_LIB`
to be `$(ZEEP_DIR)/include` and `$(ZEEP_DIR)/lib` respectively.
