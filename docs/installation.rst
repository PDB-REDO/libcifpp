Installation
============

This is a short overview. The full instructions, including all CMake options, are in the `README <https://github.com/PDB-REDO/libcifpp#installation>`_.

Requirements
------------

* A C++23 compiler (gcc >= 15, clang >= 20 or MSVC 2022)
* CMake 3.23 or newer
* zlib and liblzma development files
* PCRE2 development files (optional; a bundled fallback is used when not found)
* Eigen 3.4 or newer (optional, can be downloaded by CMake)
* SQLite 3 development files (optional; otherwise the bundled amalgamation is used)

The gxrio library is included as a git submodule, so clone with ``--recursive`` or run ``git submodule update --init``. If gxrio is already installed it will be used instead.

Building
--------

.. code-block:: console

    cmake -S . -B build
    cmake --build build
    ctest --test-dir build
    cmake --install build

The CCD data file can be downloaded during installation with ``-DCIFPP_DOWNLOAD_CCD=ON``. An update script can be installed with ``-DCIFPP_INSTALL_UPDATE_SCRIPT=ON``; see :doc:`/compound` for its configuration.
