# SPDX-License-Identifier: BSD-2-Clause

# Copyright (c) 2021 NKI/AVL, Netherlands Cancer Institute

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.

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

cmake_minimum_required(VERSION 3.16)

# set the project name
project(libcifpp VERSION 5.1.2 LANGUAGES CXX)

list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(GNUInstallDirs)
include(CheckFunctionExists)
include(CheckIncludeFiles)
include(CheckLibraryExists)
include(CMakePackageConfigHelpers)
include(CheckCXXSourceCompiles)
include(Dart)
include(GenerateExportHeader)

set(CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set(CMAKE_CXX_VISIBILITY_PRESET hidden)
# set(CMAKE_VISIBILITY_INLINES_HIDDEN 1)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter -Wno-missing-field-initializers")
elseif(MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
endif()

# Building shared libraries?
option(BUILD_SHARED_LIBS "Build a shared library instead of a static one" OFF)

# We do not want to write an export file for all our symbols...
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

# Optionally build a version to be installed inside CCP4
option(BUILD_FOR_CCP4 "Build a version to be installed in CCP4" OFF)

# Lots of code depend on the availability of the components.cif file
option(CIFPP_DOWNLOAD_CCD "Download the CCD file components.cif during installation" ON)

# An optional cron script can be installed to keep the data files up-to-date
if(UNIX)
	option(CIFPP_INSTALL_UPDATE_SCRIPT "Install the script to update CCD and dictionary files" ON)
endif()

# When CCP4 is sourced in the environment, we can recreate the symmetry operations table
if(EXISTS "$ENV{CCP4}")
	if(EXISTS "$ENV{CLIBD}/syminfo.lib")
		option(CIFPP_RECREATE_SYMOP_DATA "Recreate SymOp data table in case it is out of date" ON)
	else()
		set(CIFPP_RECREATE_SYMOP_DATA OFF)
		message(WARNING "Symop data table recreation requested, but file syminfo.lib was not found in $ENV{CLIBD}")
	endif()
else()
	set(CIFPP_RECREATE_SYMOP_DATA OFF)
	message("Not trying to recreate symop_table_data.hpp since CCP4 is not defined")
endif()

# Unit tests
option(ENABLE_TESTING "Build test exectuables" OFF)

if(BUILD_FOR_CCP4)
	if("$ENV{CCP4}" STREQUAL "" OR NOT EXISTS $ENV{CCP4})
		message(FATAL_ERROR "A CCP4 built was requested but CCP4 was not sourced")
	else()
		list(APPEND CMAKE_MODULE_PATH "$ENV{CCP4}")
		list(APPEND CMAKE_PREFIX_PATH "$ENV{CCP4}")
		set(CMAKE_INSTALL_PREFIX "$ENV{CCP4}")
		set(CMAKE_INSTALL_FULL_DATADIR "${CMAKE_INSTALL_PREFIX}/share/libcifpp")

		if(WIN32)
			set(BUILD_SHARED_LIBS ON)
		endif()
	endif()
endif()

if(WIN32)
    if(${CMAKE_SYSTEM_VERSION} GREATER_EQUAL 10) # Windows 10
        add_definitions(-D _WIN32_WINNT=0x0A00)
    elseif(${CMAKE_SYSTEM_VERSION} EQUAL 6.3) # Windows 8.1
        add_definitions(-D _WIN32_WINNT=0x0603)
    elseif(${CMAKE_SYSTEM_VERSION} EQUAL 6.2) # Windows 8
        add_definitions(-D _WIN32_WINNT=0x0602)
    elseif(${CMAKE_SYSTEM_VERSION} EQUAL 6.1) # Windows 7
        add_definitions(-D _WIN32_WINNT=0x0601)
    elseif(${CMAKE_SYSTEM_VERSION} EQUAL 6.0) # Windows Vista
        add_definitions(-D _WIN32_WINNT=0x0600)
    else() # Windows XP (5.1)
        add_definitions(-D _WIN32_WINNT=0x0501)
    endif()

	add_definitions(-DNOMINMAX)
endif()

if(MSVC)
    # make msvc standards compliant...
    add_compile_options(/permissive-)

	if(BUILD_SHARED_LIBS)
		set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
	else()
		set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
	endif()
endif()

# Libraries

# Start by finding out if std:regex is usable. Note that the current
# implementation in GCC is not acceptable, it crashes on long lines.
# The implementation in libc++ (clang) and MSVC seem to be OK.
check_cxx_source_compiles("
#include <iostream>
#ifndef __GLIBCXX__
#error
#endif
int main(int argc, char *argv[]) { return 0; }" GXX_LIBSTDCPP)

if(GXX_LIBSTDCPP)
	message(STATUS "Testing for known regex bug, since you're using GNU libstdc++")

	try_run(STD_REGEX_RUNNING STD_REGEX_COMPILING
		${CMAKE_CURRENT_BINARY_DIR}/test ${PROJECT_SOURCE_DIR}/cmake/test-rx.cpp)

	if(STD_REGEX_RUNNING STREQUAL FAILED_TO_RUN)
		message(STATUS "You are probably trying to compile using the g++ standard library which contains a crashing std::regex implementation. Will use boost::regex instead")

		add_subdirectory(regex EXCLUDE_FROM_ALL)
		set(BOOST_REGEX ON)
	endif()
endif()

set(CMAKE_THREAD_PREFER_PTHREAD)
set(THREADS_PREFER_PTHREAD_FLAG)
find_package(Threads)

if(MSVC)
	# Avoid linking the shared library of zlib
	# Search ZLIB_ROOT first if it is set.
	if(ZLIB_ROOT)
		set(_ZLIB_SEARCH_ROOT PATHS ${ZLIB_ROOT} NO_DEFAULT_PATH)
		list(APPEND _ZLIB_SEARCHES _ZLIB_SEARCH_ROOT)
	endif()

	# Normal search.
	set(_ZLIB_x86 "(x86)")
	set(_ZLIB_SEARCH_NORMAL
		PATHS "[HKEY_LOCAL_MACHINE\\SOFTWARE\\GnuWin32\\Zlib;InstallPath]"
		"$ENV{ProgramFiles}/zlib"
		"$ENV{ProgramFiles${_ZLIB_x86}}/zlib")
	unset(_ZLIB_x86)
	list(APPEND _ZLIB_SEARCHES _ZLIB_SEARCH_NORMAL)

	foreach(search ${_ZLIB_SEARCHES})
		find_library(ZLIB_LIBRARY NAMES zlibstatic NAMES_PER_DIR ${${search}} PATH_SUFFIXES lib)
	endforeach()
endif()

find_package(ZLIB REQUIRED)

find_package(Eigen3 REQUIRED)

include(FindFilesystem)
list(APPEND CIFPP_REQUIRED_LIBRARIES ${STDCPPFS_LIBRARY})

include(FindAtomic)
list(APPEND CIFPP_REQUIRED_LIBRARIES ${STDCPPATOMIC_LIBRARY})

if(MSVC)
	# this dependency can go once MSVC supports std::experimental::is_detected
	find_package(zeep 5.1.8 REQUIRED)
	list(APPEND CIFPP_REQUIRED_LIBRARIES zeep::zeep)
endif()

# Create a revision file, containing the current git version info
include(VersionString)
write_version_header(${PROJECT_SOURCE_DIR}/src/ LIB_NAME "LibCIFPP")

# SymOp data table
if(CIFPP_RECREATE_SYMOP_DATA)
	# The tool to create the table
	add_executable(symop-map-generator "${PROJECT_SOURCE_DIR}/src/symop-map-generator.cpp")

	add_custom_command(
		OUTPUT ${PROJECT_SOURCE_DIR}/src/symop_table_data.hpp
		COMMAND $<TARGET_FILE:symop-map-generator> $ENV{CLIBD}/syminfo.lib $ENV{CLIBD}/symop.lib ${PROJECT_SOURCE_DIR}/src/symop_table_data.hpp
	)

	add_custom_target(
		OUTPUT ${PROJECT_SOURCE_DIR}/src/symop_table_data.hpp
		DEPENDS symop-map-generator "$ENV{CLIBD}/syminfo.lib" "$ENV{CLIBD}/symop.lib"
	)
endif()

# Sources
set(project_sources
	${PROJECT_SOURCE_DIR}/src/category.cpp
	${PROJECT_SOURCE_DIR}/src/condition.cpp
	${PROJECT_SOURCE_DIR}/src/datablock.cpp
	${PROJECT_SOURCE_DIR}/src/dictionary_parser.cpp
	${PROJECT_SOURCE_DIR}/src/file.cpp
	${PROJECT_SOURCE_DIR}/src/item.cpp
	${PROJECT_SOURCE_DIR}/src/parser.cpp
	${PROJECT_SOURCE_DIR}/src/row.cpp
	${PROJECT_SOURCE_DIR}/src/validate.cpp
	${PROJECT_SOURCE_DIR}/src/text.cpp
	${PROJECT_SOURCE_DIR}/src/utilities.cpp

	${PROJECT_SOURCE_DIR}/src/atom_type.cpp
	${PROJECT_SOURCE_DIR}/src/compound.cpp
	${PROJECT_SOURCE_DIR}/src/point.cpp
	${PROJECT_SOURCE_DIR}/src/symmetry.cpp

	${PROJECT_SOURCE_DIR}/src/model.cpp

	${PROJECT_SOURCE_DIR}/src/pdb/cif2pdb.cpp
	${PROJECT_SOURCE_DIR}/src/pdb/pdb2cif.cpp
	${PROJECT_SOURCE_DIR}/src/pdb/pdb2cif_remark_3.cpp
	${PROJECT_SOURCE_DIR}/src/pdb/tls.cpp
)

set(project_headers
	${PROJECT_SOURCE_DIR}/include/cif++.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/utilities.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/item.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/datablock.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/file.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/validate.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/iterator.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/parser.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/forward_decl.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/dictionary_parser.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/condition.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/category.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/row.hpp

	${PROJECT_SOURCE_DIR}/include/cif++/atom_type.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/compound.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/point.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/symmetry.hpp

	${PROJECT_SOURCE_DIR}/include/cif++/model.hpp

	${PROJECT_SOURCE_DIR}/include/cif++/pdb/cif2pdb.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/pdb/io.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/pdb/pdb2cif.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/pdb/pdb2cif_remark_3.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/pdb/tls.hpp
)

add_library(cifpp ${project_sources} ${project_headers} ${PROJECT_SOURCE_DIR}/src/symop_table_data.hpp)
add_library(cifpp::cifpp ALIAS cifpp)
generate_export_header(cifpp EXPORT_FILE_NAME cif++/exports.hpp)

if(BOOST_REGEX)
	target_compile_definitions(cifpp PRIVATE USE_BOOST_REGEX=1 BOOST_REGEX_STANDALONE=1)
	target_include_directories(cifpp PRIVATE regex/include)
endif()

if(MSVC)
	target_compile_definitions(cifpp PUBLIC NOMINMAX=1)
endif()

set_target_properties(cifpp PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(cifpp
	PUBLIC
	"$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include;${PROJECT_BINARY_DIR};${EIGEN3_INCLUDE_DIR}>"
	"$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

target_link_libraries(cifpp PUBLIC Threads::Threads ZLIB::ZLIB Eigen3::Eigen ${CIFPP_REQUIRED_LIBRARIES})

if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
	target_link_options(cifpp PRIVATE -undefined dynamic_lookup)
endif(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")

if(CIFPP_DOWNLOAD_CCD)
	# download the components.cif file from CCD
	set(COMPONENTS_CIF ${PROJECT_SOURCE_DIR}/data/components.cif)

	if(NOT EXISTS ${COMPONENTS_CIF})
		if(NOT EXISTS ${PROJECT_SOURCE_DIR}/data)
			file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/data/)
		endif()

		file(DOWNLOAD https://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif ${COMPONENTS_CIF} SHOW_PROGRESS)
	endif()

	add_custom_target(COMPONENTS ALL DEPENDS ${COMPONENTS_CIF})
endif()

# Installation directories
set(CIFPP_DATA_DIR "${CMAKE_INSTALL_FULL_DATADIR}/libcifpp")
target_compile_definitions(cifpp PUBLIC DATA_DIR="${CIFPP_DATA_DIR}")

if(UNIX)
	set(CIFPP_CACHE_DIR "${CMAKE_INSTALL_FULL_LOCALSTATEDIR}/cache/libcifpp")
	target_compile_definitions(cifpp PUBLIC CACHE_DIR="${CIFPP_CACHE_DIR}")

	set(CIFPP_ETC_DIR "${CMAKE_INSTALL_FULL_SYSCONFDIR}")
endif()

# Install rules
install(TARGETS cifpp
	EXPORT cifppTargets
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(MSVC AND BUILD_SHARED_LIBS)
	install(
		FILES $<TARGET_PDB_FILE:${PROJECT_NAME}>
		DESTINATION ${CMAKE_INSTALL_LIBDIR}
		OPTIONAL)
endif()

install(EXPORT cifppTargets
	FILE "cifppTargets.cmake"
	NAMESPACE cifpp::
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cifpp
)

install(
	DIRECTORY include/cif++
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	COMPONENT Devel
)

install(
	FILES include/cif++.hpp
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
	COMPONENT Devel
)

install(
	FILES ${PROJECT_BINARY_DIR}/cif++/exports.hpp
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cif++
	COMPONENT Devel
)

install(FILES
	${PROJECT_SOURCE_DIR}/rsrc/mmcif_ddl.dic
	${PROJECT_SOURCE_DIR}/rsrc/mmcif_pdbx.dic
	${PROJECT_SOURCE_DIR}/rsrc/mmcif_ma.dic
	${COMPONENTS_CIF}
	DESTINATION ${CIFPP_DATA_DIR}
)

if(${CIFPP_CACHE_DIR})
	install(FILES
		${PROJECT_SOURCE_DIR}/rsrc/mmcif_ddl.dic
		${PROJECT_SOURCE_DIR}/rsrc/mmcif_pdbx.dic
		${PROJECT_SOURCE_DIR}/rsrc/mmcif_ma.dic
		${COMPONENTS_CIF}
		DESTINATION ${CIFPP_CACHE_DIR}
	)
endif()

set(CONFIG_TEMPLATE_FILE ${PROJECT_SOURCE_DIR}/cmake/cifppConfig.cmake.in)

configure_package_config_file(
	${CONFIG_TEMPLATE_FILE}
	${CMAKE_CURRENT_BINARY_DIR}/cifpp/cifppConfig.cmake
	INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cifpp
	PATH_VARS CIFPP_DATA_DIR
)

install(FILES
	"${CMAKE_CURRENT_BINARY_DIR}/cifpp/cifppConfig.cmake"
	"${CMAKE_CURRENT_BINARY_DIR}/cifpp/cifppConfigVersion.cmake"
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cifpp
	COMPONENT Devel
)

set(cifpp_MAJOR_VERSION ${CMAKE_PROJECT_VERSION_MAJOR})
set_target_properties(cifpp PROPERTIES
	VERSION ${PROJECT_VERSION}
	SOVERSION ${cifpp_MAJOR_VERSION}
	INTERFACE_cifpp_MAJOR_VERSION ${cifpp_MAJOR_VERSION})

set_property(TARGET cifpp APPEND PROPERTY
	COMPATIBLE_INTERFACE_STRING cifpp_MAJOR_VERSION
)

write_basic_package_version_file(
	"${CMAKE_CURRENT_BINARY_DIR}/cifpp/cifppConfigVersion.cmake"
	VERSION ${PROJECT_VERSION}
	COMPATIBILITY AnyNewerVersion
)

# pkgconfig support
set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix ${CMAKE_INSTALL_PREFIX})
set(libdir ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(includedir ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libcifpp.pc.in
	${CMAKE_CURRENT_BINARY_DIR}/libcifpp.pc.in @ONLY)
file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/libcifpp.pc
	INPUT ${CMAKE_CURRENT_BINARY_DIR}/libcifpp.pc.in)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libcifpp.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

if(ENABLE_TESTING)
	enable_testing()

	find_package(Boost REQUIRED)

	list(APPEND CIFPP_tests unit-v2 unit-3d format model rename-compound sugar spinner)

	foreach(CIFPP_TEST IN LISTS CIFPP_tests)
		set(CIFPP_TEST "${CIFPP_TEST}-test")
		set(CIFPP_TEST_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/test/${CIFPP_TEST}.cpp")

		add_executable(${CIFPP_TEST} ${CIFPP_TEST_SOURCE})

		target_link_libraries(${CIFPP_TEST} PRIVATE Threads::Threads cifpp::cifpp Boost::boost)

		if(MSVC)
			# Specify unwind semantics so that MSVC knowns how to handle exceptions
			target_compile_options(${CIFPP_TEST} PRIVATE /EHsc)
		endif()

		add_custom_target("run-${CIFPP_TEST}" DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/Run${CIFPP_TEST}.touch ${CIFPP_TEST})

		add_custom_command(
			OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/Run${CIFPP_TEST}.touch
			COMMAND $<TARGET_FILE:${CIFPP_TEST}> -- ${CMAKE_CURRENT_SOURCE_DIR}/test)

		add_test(NAME ${CIFPP_TEST}
			COMMAND $<TARGET_FILE:${CIFPP_TEST}> -- ${CMAKE_CURRENT_SOURCE_DIR}/test)
	endforeach()
endif()

# Optionally install the update scripts for CCD and dictionary files
if(CIFPP_INSTALL_UPDATE_SCRIPT)
	if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux" OR ${CMAKE_SYSTEM_NAME} STREQUAL "GNU")
		set(CIFPP_CRON_DIR "${CIFPP_ETC_DIR}/cron.weekly")
	elseif(UNIX) # assume all others are like FreeBSD...
		set(CIFPP_CRON_DIR "${CIFPP_ETC_DIR}/periodic/weekly")
	else()
		message(FATAL_ERROR "Don't know where to install the update script")
	endif()

	configure_file(${PROJECT_SOURCE_DIR}/tools/update-libcifpp-data.in update-libcifpp-data @ONLY)
	install(
		FILES ${CMAKE_CURRENT_BINARY_DIR}/update-libcifpp-data
		DESTINATION ${CIFPP_CRON_DIR}
		PERMISSIONS OWNER_EXECUTE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ
	)

	install(DIRECTORY DESTINATION ${CIFPP_CACHE_DIR})
	install(DIRECTORY DESTINATION "${CIFPP_ETC_DIR}/libcifpp/cache-update.d")

	# a config to, to make it complete
	if(NOT EXISTS "${CIFPP_ETC_DIR}/libcifpp.conf")
		file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/libcifpp.conf [[# Uncomment the next line to enable automatic updates
# update=true
]])
		install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libcifpp.conf DESTINATION "${CIFPP_ETC_DIR}")
		install(CODE "message(\"A configuration file has been written to ${CIFPP_ETC_DIR}/libcifpp.conf, please edit this file to enable automatic updates\")")
	endif()

	target_compile_definitions(cifpp PUBLIC CACHE_DIR="${CIFPP_CACHE_DIR}")
endif()

set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_SOURCE_TGZ ON)
set(CPACK_SOURCE_TBZ2 OFF)
set(CPACK_SOURCE_TXZ OFF)
set(CPACK_SOURCE_TZ OFF)
set(CPACK_SOURCE_IGNORE_FILES "/data/components.cif;/build;/.vscode;/.git")
set(CPACK_PACKAGE_FILE_NAME "${PROJECT_NAME}-${PROJECT_VERSION}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME ${CPACK_PACKAGE_FILE_NAME})
include(CPack)
