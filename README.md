libcif++
========

This library contains code to work with mmCIF and PDB files.


Requirements
------------

The code for this library was written in C++17. You therefore need a
recent compiler to build it. For the development gcc 9.3 and clang 9.0
have been used.

Other requirements are:

- GNU make version 4.1 or higher.
- Boost libraries, at least version 1.71
- [mrc](https://github.com/mhekkel/mrc), a resource compiler that
  allows including data files into the executable making them easier to
  install. Strictly this is optional, but at the expense of a lot of
  functionality.

Building
--------

Simply configure, make and make install.

There's one configure flag that might be of interest: if you specify
DEBUG=1 the make will create a debug version by default. Otherwise,
you can run make with DEBUG=1 to create a debug version.
