# libcifpp

This library contains code to work with mmCIF and legacy PDB files.

## Documentation

The documentation can be found at https://www.hekkelman.com/libcifpp-doc/

## Synopsis

```c++
// A simple program counting residues with an OXT atom

#include <filesystem>
#include <iostream>

#include <cif++.hpp>

namespace fs = std::filesystem;

int main(int argc, char *argv[])
{
    if (argc != 2)
        exit(1);

    // Read file, can be PDB or mmCIF and can even be compressed with gzip.
    cif::file file = cif::pdb::read(argv[1]);

    if (file.empty())
    {
        std::cerr << "Empty file" << std::endl;
        exit(1);
    }

    // Take the first datablock in the file
    auto &db = file.front();

    // Use the atom_site category
    auto &atom_site = db["atom_site"];

    // Count the atoms with atom-id "OXT"
    auto n = atom_site.count(cif::key("label_atom_id") == "OXT");

    std::cout << "File contains " << atom_site.size() << " atoms of which "
              << n << (n == 1 ? " is" : " are") << " OXT" << std::endl
              << "residues with an OXT are:" << std::endl;

    // Loop over all atoms with atom-id "OXT" and print out some info.
    // That info is extracted using structured binding in C++
    for (const auto &[asym, comp, seqnr] :
            atom_site.find<std::string, std::string, int>(
                cif::key("label_atom_id") == "OXT",
                "label_asym_id", "label_comp_id", "label_seq_id"))
    {
        std::cout << asym << ' ' << comp << ' ' << seqnr << std::endl;
    }

    return 0;
}
```

## Installation

You might be able to use libcifpp from a package manager used by your
OS distribution. But most likely this package will be out-of-date.
Therefore it is recommended to build *libcifpp* from code. It is not
hard to do.

### Requirements

The code for this library was written in C++17. You therefore need a
recent compiler to build it. For the development gcc 9.4 and clang 9.0
have been used as well as MSVC version 2019.

Other requirements are:

- [cmake](https://cmake.org) A build tool.
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
- [boost](https://www.boost.org).

When building using MS Visual Studio, you will also need [libzeep](https://github.com/mhekkel/libzeep)
since MSVC does not yet provide a C++ template required by libcifpp.

The Boost libraries are only needed in case you want to build the test
code or if you are using GCC. That last condition is due to a long
standing bug in the implementation of std::regex. It simply crashes
on the regular expressions used in the mmcif_pdbx dictionary and so
we use the boost regex implementation instead.

### Building

Building the code is as simple as typing:

```console
 git clone https://github.com/PDB-REDO/libcifpp.git --recurse-submodules
 cd libcifpp
 cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$HOME/.local -DCMAKE_BUILD_TYPE=Release
 cmake --build build
 cmake --install build
```

This checks out the source code from github, creates a new directory
where cmake stores its files. Run a configure, build the code and then
it installs the library and auxiliary files.

Tests are created by default, and to test the code you can run:

```console
 cmake --build build
 ctest --test-dir build
```
