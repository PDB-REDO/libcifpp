libcif++ and pdb-redo tools
===========================

This is the README file for libcif++ and the suite of PDB-REDO tools using this library.

### Requirements

The code for this library was written in C++14. You therefore need a recent compiler to be able to build it.

Other requirements are:

- Boost libraries, the current version was developed using version 1.65
- [mrc](https://github.com/mhekkel/mrc), a resource compiler that allows including data files into the executable making them easier to install. Strictly this is optional, but at the expense of a lot of functionality.
- [newuoa-cpp](https://github.com/elsid/newuoa-cpp), required to calculate atom radii.
- [libzeep](https://github.com/mhekkel/libzeep), a library that contains a full validating XML parser as well as a complete HTTP and SOAP server implementation.
- [nlohmann/json](https://github.com/nlohmann/json), a header only library to parse and write JSON. This requirement will be removed in a future release (once libzeep 4 is out of beta).

### Building

There are two makefiles, one located in the directory libcif++ and one in tools.

Both makefiles will include a *make.config* file (which will be generated if it doesn't exist). This configuration file can be used to override local settings, e.g. the location of certain libraries.