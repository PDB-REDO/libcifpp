libcifpp
========

This library contains code to work with mmCIF and PDB files.

Requirements
------------

The code for this library was written in C++17. You therefore need a
recent compiler to build it. For the development gcc 9.3 and clang 9.0
have been used as well as MSVC version 2019.

Other requirements are:

- Boost libraries, at least version 1.70
- [mrc](https://github.com/mhekkel/mrc), a resource compiler that
  allows including data files into the executable making them easier to
  install. Strictly this is optional, but at the expense of functionality.

Building
--------

This library uses [cmake](https://cmake.org). The usual way of building
and installing is to create a `build` directory and run cmake there.

On linux e.g. you would issue the following commands:

```
	git clone https://github.com/PDB-REDO/libcifpp.git
	cd libcifpp
	mkdir build
	cd build
	cmake ..
	cmake --build . --config Release
	ctest -C Release
	cmake --install .
```
This checks out the source code from github, creates a new directory
where cmake stores its files. Run a configure, build the code and run
tests. And then it installs the library and auxiliary files.

The default is to install everything in `$HOME/.local` on Linux and
`%LOCALAPPDATA%` on Windows (the AppData/Local folder in your home directory).

