libcifpp
========

This library contains code to work with mmCIF and PDB files.

Requirements
------------

The code for this library was written in C++17. You therefore need a
recent compiler to build it. For the development gcc 9.4 and clang 9.0
have been used as well as MSVC version 2019.

Other requirements are:

- [mrc](https://github.com/mhekkel/mrc), a resource compiler that
  allows including data files into the executable making them easier to
  install. Strictly speaking this is optional, but at the expense of
  functionality.
- [libeigen](https://eigen.tuxfamily.org/index.php?title=Main_Page), a
  library to do amongst others matrix calculations. This usually can be
  installed using your package manager, in Debian/Ubuntu it is called
  `libeigen3-dev`
- zlib, the development version of this library. On Debian/Ubuntu this
  is the package `zlib1g-dev`.
- [boost](https://www.boost.org). The boost libraries are only needed if
  you want to build the testing code.

When building using MS Visual Studio, you will also need [libzeep](https://github.com/mhekkel/libzeep)
since MSVC does not yet provide a C++ template required by libcifpp.

Building
--------

This library uses [cmake](https://cmake.org). The usual way of building
and installing is to create a `build` directory and run cmake there.

On linux e.g. you would issue the following commands to build and install
libcifpp in your `$HOME/.local` folder:

```bash
 git clone https://github.com/PDB-REDO/libcifpp.git --recurse-submodules
 cd libcifpp
 cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$HOME/.local -DCMAKE_BUILD_TYPE=Release
 cmake --build build
 cmake --install build
```

This checks out the source code from github, creates a new directory
where cmake stores its files. Run a configure, build the code and then
it installs the library and auxiliary files.

If you want to run the tests before installing, you should add `-DENABLE_TESTING=ON`
to the first cmake command.
