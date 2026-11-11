# SPDX-License-Identifier: BSD-2-Clause
# 
# Copyright (c) 2026 NKI/AVL, Netherlands Cancer Institute
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

find_package(Doxygen REQUIRED)
find_package(Sphinx REQUIRED)

# Find all the public headers
# get_target_property(CIFPP_PUBLIC_HEADER_DIR libCIFPP INTERFACE_INCLUDE_DIRECTORIES)
set(CIFPP_PUBLIC_HEADER_DIR ${PROJECT_SOURCE_DIR}/include)
file(GLOB_RECURSE CIFPP_PUBLIC_HEADERS ${CIFPP_PUBLIC_HEADER_DIR}/*.hpp)

set(DOXYGEN_INPUT_DIR ${CIFPP_PUBLIC_HEADER_DIR})
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/xml)
set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/index.xml)
set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

# Replace variables inside @@ with the current values
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

add_custom_command(
    OUTPUT ${DOXYGEN_OUTPUT_DIR}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${DOXYGEN_OUTPUT_DIR})

add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
    DEPENDS ${CIFPP_PUBLIC_HEADERS} ${DOXYGEN_OUTPUT_DIR} ${DOXYFILE_OUT}
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
    MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
    COMMENT "Generating docs")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in ${CMAKE_CURRENT_SOURCE_DIR}/conf.py @ONLY)

set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)

add_custom_target("Sphinx-cifpp" ALL
    COMMAND ${CMAKE_COMMAND} -E rm -rf ${CMAKE_CURRENT_SOURCE_DIR}/api
    COMMAND ${SPHINX_EXECUTABLE} -b html
    -Dbreathe_projects.cifpp=${DOXYGEN_OUTPUT_DIR}
    ${SPHINX_SOURCE} ${SPHINX_BUILD}
    DEPENDS ${DOXYGEN_INDEX_FILE}
    BYPRODUCTS ${CMAKE_CURRENT_SOURCE_DIR}/api
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating documentation with Sphinx")

set_property(TARGET Sphinx-cifpp APPEND PROPERTY ADDITIONAL_CLEAN_FILES
    ${DOXYGEN_OUTPUT_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/api)

install(
    DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/sphinx/
    DESTINATION ${CMAKE_INSTALL_DOCDIR}
    PATTERN .doctrees EXCLUDE
    PATTERN .buildinfo EXCLUDE)
