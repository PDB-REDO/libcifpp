[![github CI](https://github.com/pdb-redo/libcifpp/actions/workflows/cmake-multi-platform.yml/badge.svg)](https://github.com/pdb-redo/libcifpp/actions)
[![GitHub License](https://img.shields.io/github/license/pdb-redo/libcifpp)](https://github.com/pdb-redo/libcifpp/LICENSE)

# libcifpp

As the name implies, this library was originally written to work with mmCIF files
using C++ as programming language. The design of this library leanes heavily on
the structure of CIF files. These files can be thought of as a text dump of a
relational databank with, often but not always, a very strict schema describing
the data. These schema's are called dictionaries.

Using information from the content of a mmCIF file and an optional schema,
libcifpp allows you to access the data in the file as a collection of datablock
each containing a collection of categories with rows of data. The categories can
be searched for data using queries written in regular C++ syntax. When a dictionary
was specified, inserted data is checked for validity. Likewise removal of data
may result in cascaded removal of linked data in other categories using
parent/child relationship information.

Since there were still many programs using the legacy PDB format at the time
development started, a layer was added that converts data to and from PDB format
into mmCIF format. This means you can manipulate PDB files as if they were
normal mmCIF files.

Apart from this basic functionality, libcifpp also offers code to help with
symmetry calculations, 3d manipulations and obtaining information from the CCD
[Chemical Component Dictionary](https://www.wwpdb.org/data/ccd).

## Documentation

The documentation can be found at [github.io](https://pdb-redo.github.io/libcifpp/)

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
        std::cerr << "Empty file\n";
        exit(1);
    }

    // Take the first datablock in the file
    auto &db = file.front();

    // Use the atom_site category
    auto &atom_site = db["atom_site"];

    // Count the atoms with atom-id "OXT"
    auto n = atom_site.count(cif::key("label_atom_id") == "OXT");

    std::cout << "File contains " << atom_site.size() << " atoms of which "
              << n << (n == 1 ? " is" : " are") << " OXT\n"
              << "residues with an OXT are:\n";

    // Loop over all atoms with atom-id "OXT" and print out some info.
    // That info is extracted using structured binding in C++
    for (const auto &[asym, comp, seqnr] :
            atom_site.find<std::string, std::string, int>(
                cif::key("label_atom_id") == "OXT",
                "label_asym_id", "label_comp_id", "label_seq_id"))
    {
        std::cout << asym << ' ' << comp << ' ' << seqnr << '\n';
    }

    return 0;
}
```

## Installation

You might be able to use libcifpp from a package manager used by your
OS distribution. But most likely this package will be out-of-date.
Therefore it is recommended to build *libcifpp* from code. It is not
hard to do. But it is recommended to read the following instructions
carefully.

### Requirements

The code for this library was written in C++17. You therefore need a
recent compiler to build it. For the development gcc >= 9.4 and clang >= 9.0
have been used as well as MSVC version 2019.

The other requirement you really need to have installed on your computer
is a version of [CMake](https://cmake.org). For now the minimum version
is 3.16 but that may soon change into a higher version. You should also
install the gui version of CMake to set build options easily, on Debian
I prefer to use the curses version installed with `cmake-curses-gui`.

It is very useful to have [mrc](https://github.com/mhekkel/mrc) available.
However, this is only an option if you use Windows or an operating system
using the ELF executable format (i.e. Linux or FreeBSD). MRC is a resource
compiler that allows including data files into the executable making them
easier to install.

Other libraries you might want to install beforehand are:

- [libeigen](https://eigen.tuxfamily.org/index.php?title=Main_Page), a
  library to do amongst others matrix calculations. This usually can be
  installed using your package manager, in Debian/Ubuntu it is called
  `libeigen3-dev`
- [zlib](https://github.com/madler/zlib), the development version of this
  library. On Debian/Ubuntu this is the package `zlib1g-dev`.
- [boost](https://www.boost.org), in Debian/Ubuntu this is `libboost-dev`.
  
  The Boost libraries are only needed in case you are using GCC due to a long
  standing bug in GNU's implementation of std::regex. It simply crashes
  on the regular expressions used in the mmcif_pdbx dictionary and so
  we use the boost regex implementation instead.

### Building

First you need to download the code:

```console
 git clone https://github.com/PDB-REDO/libcifpp.git
 cd libcifpp
```

You should start by considering where to install libcifpp. If you have
sufficient permissions on your computer you perhaps should use the
default but libcifpp can be configured to be installed anywhere
including e.g. *$HOME/.local*.

Next step is to configure, for this use the CMake gui application. If you
installed the curses version of cmake you can type `ccmake`. On Windows
you can use `cmake-gui.exe`.

To install in the default location:

```console
ccmake -S . -B build
```

To install elsewhere, e.g. *$HOME/.local*:

```console
ccmake -S . -B build -DCMAKE_INSTALL_PREFIX=$HOME/.local
```

In the cmake window, start the configure command (use button or press 'c').
After the first configure step you will see a list of settable options.
Alter these to match your preferences. Most options are self explaining
and contain a description. Some may need a bit more explanation:

- CIFPP_DATA_DIR, this directory will be used to store initial versions
  of the mmcif_pdbx dictionary as well as the optional CCD file.

- CIFPP_DOWNLOAD_CCD

  The CCD file is huge and perhaps you think you don't
  need it. In that case you can leave this OFF. But that will limit the
  use cases.

- CIFPP_INSTALL_UPDATE_SCRIPT
  
  The files in CIFPP_DATA_DIR are quickly becoming out of date. On
  FreeBSD and Linux you can install a script that updates these files
  on a weekly basis.

- CIFPP_CRON_DIR
  
  The directory where the update script is to be installed.

- CIFPP_ETC_DIR
  
  The update script will only work if the file called *libcifpp.conf*
  in this *etc* directory will contain an uncommented line with

```console
update=true
```

- CIFPP_CACHE_DIR

  When you installed and enabled the update script, new files are
  written to this directory.

- CIFPP_RECREATE_SYMOP_DATA
  
  If you had CCP4 sourced into your environment, this option allows
  you to recreate the symop data file.

- BUILD_FOR_CCP4
  
  Build a special version of libcifpp to be installed in the CCP4
  environment.

After setting these options you can run the configure step again and
then use generate to create the makefiles.

Building and installing is then as simple as:

```console
cmake --build build
cmake --install build
```

If this fails due to lack of permissions, you can try:

```console
sudo cmake --install build
```

Tests are created by default, and to test the code you can run:

```console
ctest --test-dir build
```
