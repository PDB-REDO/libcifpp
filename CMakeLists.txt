# SPDX-License-Identifier: BSD-2-Clause

# Copyright (c) 2021 NKI/AVL, Netherlands Cancer Institute

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

cmake_minimum_required(VERSION 3.16)

# set the project name
project(
	libcifpp
	VERSION 7.0.3
	LANGUAGES CXX)

list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(CheckFunctionExists)
include(CheckIncludeFiles)
include(CheckLibraryExists)
include(CMakePackageConfigHelpers)
include(CheckCXXSourceCompiles)
include(GenerateExportHeader)
include(CTest)
include(FetchContent)
include(ExternalProject)

set(CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# When building with ninja-multiconfig, build both debug and release by default
if(CMAKE_GENERATOR STREQUAL "Ninja Multi-Config")
	set(CMAKE_CROSS_CONFIGS "Debug;Release")
	set(CMAKE_DEFAULT_CONFIGS "Debug;Release")
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS
		"${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter -Wno-missing-field-initializers"
	)
elseif(MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
endif()

# Build documentation?
option(BUILD_DOCUMENTATION "Build the documentation" OFF)

# Optionally build a version to be installed inside CCP4
option(BUILD_FOR_CCP4 "Build a version to be installed in CCP4")

# Building shared libraries?
if(NOT(BUILD_FOR_CCP4 AND WIN32))
	option(BUILD_SHARED_LIBS "Build a shared library instead of a static one" OFF)
endif()

if(BUILD_FOR_CCP4)
	unset(CIFPP_DOWNLOAD_CCD)
	unset(CIFPP_INSTALL_UPDATE_SCRIPT)
else()
	# Lots of code depend on the availability of the components.cif file
	option(CIFPP_DOWNLOAD_CCD
		"Download the CCD file components.cif during installation" ON)

	# An optional cron script can be installed to keep the data files up-to-date
	if(UNIX AND NOT APPLE)
		option(CIFPP_INSTALL_UPDATE_SCRIPT
			"Install the script to update CCD and dictionary files" ON)
	endif()
endif()

# When CCP4 is sourced in the environment, we can recreate the symmetry
# operations table
if(EXISTS "$ENV{CCP4}/lib/data/syminfo.lib")
	option(CIFPP_RECREATE_SYMOP_DATA
		"Recreate SymOp data table in case it is out of date" ON)
endif()

# CCP4 build
if(BUILD_FOR_CCP4)
	if("$ENV{CCP4}" STREQUAL "" OR NOT EXISTS $ENV{CCP4})
		message(FATAL_ERROR "A CCP4 built was requested but CCP4 was not sourced")
	else()
		list(PREPEND CMAKE_MODULE_PATH "$ENV{CCP4}")
		list(PREPEND CMAKE_PREFIX_PATH "$ENV{CCP4}")
		set(CMAKE_INSTALL_PREFIX "$ENV{CCP4}")

		if(WIN32)
			set(BUILD_SHARED_LIBS ON)
		endif()
	endif()
endif()

# Now include the GNUInstallDirs module
include(GNUInstallDirs)

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

	# We do not want to write an export file for all our symbols...
	set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

if(MSVC)
	# make msvc standards compliant...
	add_compile_options(/permissive- /bigobj)
	add_link_options(/NODEFAULTLIB:library)

	if(BUILD_SHARED_LIBS)
		set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
	else()
		set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
	endif()
endif()

# Libraries

# Start by finding out if std:regex is usable. Note that the current
# implementation in GCC is not acceptable, it crashes on long lines. The
# implementation in libc++ (clang) and MSVC seem to be OK.
check_cxx_source_compiles(
	"
#include <iostream>
#ifndef __GLIBCXX__
#error
#endif
int main(int argc, char *argv[]) { return 0; }"
	GXX_LIBSTDCPP)

if(GXX_LIBSTDCPP)
	message(
		STATUS "Testing for known regex bug, since you're using GNU libstdc++")

	try_run(STD_REGEX_RUNNING STD_REGEX_COMPILING
		${CMAKE_CURRENT_BINARY_DIR}/test
		${CMAKE_CURRENT_SOURCE_DIR}/cmake/test-rx.cpp)

	if(STD_REGEX_RUNNING STREQUAL FAILED_TO_RUN)
		message(
			STATUS
			"You are probably trying to compile using the g++ standard library which contains a crashing std::regex implementation. Will use boost::regex instead"
		)

		find_package(Boost 1.80 QUIET COMPONENTS regex)

		if(NOT Boost_FOUND)
			set(BOOST_REGEX_STANDALONE ON)

			FetchContent_Declare(
				boost-rx
				GIT_REPOSITORY https://github.com/boostorg/regex
				GIT_TAG boost-1.83.0)

			FetchContent_MakeAvailable(boost-rx)
		endif()

		set(BOOST_REGEX ON)
	endif()
endif()

set(CMAKE_THREAD_PREFER_PTHREAD)
set(THREADS_PREFER_PTHREAD_FLAG)
find_package(Threads)

if(MSVC)
	# Avoid linking the shared library of zlib Search ZLIB_ROOT first if it is
	# set.
	if(ZLIB_ROOT)
		set(_ZLIB_SEARCH_ROOT PATHS ${ZLIB_ROOT} NO_DEFAULT_PATH)
		list(APPEND _ZLIB_SEARCHES _ZLIB_SEARCH_ROOT)
	endif()

	# Normal search.
	set(_ZLIB_x86 "(x86)")
	set(_ZLIB_SEARCH_NORMAL
		PATHS "[HKEY_LOCAL_MACHINE\\SOFTWARE\\GnuWin32\\Zlib;InstallPath]"
		"$ENV{ProgramFiles}/zlib" "$ENV{ProgramFiles${_ZLIB_x86}}/zlib")
	unset(_ZLIB_x86)
	list(APPEND _ZLIB_SEARCHES _ZLIB_SEARCH_NORMAL)

	if(BUILD_FOR_CCP4)
		list(PREPEND _ZLIB_SEARCHES "$ENV{CCP4}/lib")
	endif()

	foreach(search ${_ZLIB_SEARCHES})
		find_library(
			ZLIB_LIBRARY
			NAMES zlibstatic NAMES_PER_DIR ${${search}}
			PATH_SUFFIXES lib)
	endforeach()
endif()

find_package(ZLIB REQUIRED)

# Using Eigen3 is a bit of a thing. We don't want to build it completely since
# we only need a couple of header files. Nothing special. But often, eigen3 is
# already installed and then we prefer that.
find_package(Eigen3 3.4 QUIET)

if(Eigen3_FOUND AND TARGET Eigen3::Eigen)
	get_target_property(EIGEN_INCLUDE_DIR Eigen3::Eigen
		INTERFACE_INCLUDE_DIRECTORIES)
else()
	# Create a private copy of eigen3 and populate it only, no need to build
	FetchContent_Declare(
		my-eigen3
		GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
		GIT_TAG 3.4.0)

	FetchContent_GetProperties(my-eigen3)

	if(NOT my-eigen3_POPULATED)
		FetchContent_Populate(my-eigen3)
	endif()

	set(EIGEN_INCLUDE_DIR ${my-eigen3_SOURCE_DIR})
endif()

include(FindFilesystem)
list(APPEND CIFPP_REQUIRED_LIBRARIES ${STDCPPFS_LIBRARY})

include(FindAtomic)
list(APPEND CIFPP_REQUIRED_LIBRARIES ${STDCPPATOMIC_LIBRARY})

# Create a revision file, containing the current git version info
include(VersionString)
write_version_header(${CMAKE_CURRENT_SOURCE_DIR}/src/ LIB_NAME "LibCIFPP")

# SymOp data table
if(CIFPP_RECREATE_SYMOP_DATA)
	# The tool to create the table
	add_executable(symop-map-generator
		"${CMAKE_CURRENT_SOURCE_DIR}/src/symop-map-generator.cpp")

	add_custom_command(
		OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/src/symop_table_data.hpp
		COMMAND
		$<TARGET_FILE:symop-map-generator> $ENV{CLIBD}/syminfo.lib
		$ENV{CLIBD}/symop.lib ${CMAKE_CURRENT_SOURCE_DIR}/src/symop_table_data.hpp)

	add_custom_target(
		OUTPUT
		${CMAKE_CURRENT_SOURCE_DIR}/src/symop_table_data.hpp
		DEPENDS symop-map-generator "$ENV{CLIBD}/syminfo.lib"
		"$ENV{CLIBD}/symop.lib")
endif()

# Sources
set(project_sources
	${CMAKE_CURRENT_SOURCE_DIR}/src/category.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/datablock.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/dictionary_parser.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/file.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/item.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/parser.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/row.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/validate.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/text.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/utilities.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/atom_type.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/compound.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/point.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/symmetry.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/model.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/pdb/cif2pdb.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/pdb/pdb2cif.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/pdb/pdb_record.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/pdb/pdb2cif_remark_3.hpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/pdb/pdb2cif_remark_3.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/pdb/reconstruct.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/pdb/validate-pdbx.cpp
)

set(project_headers
	include/cif++.hpp
	include/cif++/utilities.hpp
	include/cif++/item.hpp
	include/cif++/datablock.hpp
	include/cif++/exports.hpp
	include/cif++/file.hpp
	include/cif++/validate.hpp
	include/cif++/iterator.hpp
	include/cif++/parser.hpp
	include/cif++/forward_decl.hpp
	include/cif++/dictionary_parser.hpp
	include/cif++/condition.hpp
	include/cif++/category.hpp
	include/cif++/row.hpp
	include/cif++/atom_type.hpp
	include/cif++/compound.hpp
	include/cif++/point.hpp
	include/cif++/symmetry.hpp
	include/cif++/model.hpp
	include/cif++/pdb.hpp
	include/cif++/pdb/cif2pdb.hpp
	include/cif++/pdb/io.hpp
	include/cif++/pdb/pdb2cif.hpp
	include/cif++/pdb/tls.hpp)

add_library(cifpp STATIC)
add_library(cifpp::cifpp ALIAS cifpp)

target_sources(cifpp
	PRIVATE ${project_sources}
	${CMAKE_CURRENT_SOURCE_DIR}/src/symop_table_data.hpp
	PUBLIC
	FILE_SET cifpp_headers TYPE HEADERS
	BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include
	FILES ${project_headers}
)

# The code now really requires C++20
target_compile_features(cifpp PUBLIC cxx_std_20)

set(CMAKE_DEBUG_POSTFIX d)
set_target_properties(cifpp PROPERTIES DEBUG_POSTFIX "d")

generate_export_header(cifpp EXPORT_FILE_NAME
	${CMAKE_CURRENT_SOURCE_DIR}/include/cif++/exports.hpp)

if(BOOST_REGEX)
	target_compile_definitions(cifpp PRIVATE USE_BOOST_REGEX=1
		BOOST_REGEX_STANDALONE=1)
	get_target_property(BOOST_REGEX_INCLUDE_DIR Boost::regex
		INTERFACE_INCLUDE_DIRECTORIES)
endif()

if(MSVC)
	target_compile_definitions(cifpp PUBLIC NOMINMAX=1)
endif()

set_target_properties(cifpp PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(
	cifpp
	PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
	"$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
	PRIVATE "${BOOST_REGEX_INCLUDE_DIR}" "${EIGEN_INCLUDE_DIR}")

target_link_libraries(cifpp PUBLIC Threads::Threads ZLIB::ZLIB
	${CIFPP_REQUIRED_LIBRARIES})

if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
	target_link_options(cifpp PRIVATE -undefined dynamic_lookup)
endif(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")

if(CIFPP_DOWNLOAD_CCD)
	# download the components.cif file from CCD
	set(COMPONENTS_CIF ${CMAKE_CURRENT_SOURCE_DIR}/rsrc/components.cif)

	if(EXISTS ${COMPONENTS_CIF})
		file(SIZE ${COMPONENTS_CIF} CCD_FILE_SIZE)

		if(CCD_FILE_SIZE EQUAL 0)
			message(STATUS "Removing empty ${COMPONENTS_CIF} file")
			file(REMOVE "${COMPONENTS_CIF}")
		endif()
	endif()

	if(NOT EXISTS ${COMPONENTS_CIF})
		# Since the file(DOWNLOAD) command in cmake does not use compression, we try
		# to download the gzipped version and decompress it ourselves.
		find_program(GUNZIP gunzip)

		if(WIN32 OR GUNZIP STREQUAL "GUNZIP-NOTFOUND")
			file(
				DOWNLOAD https://files.wwpdb.org/pub/pdb/data/monomers/components.cif
				${COMPONENTS_CIF}
				SHOW_PROGRESS
				STATUS CCD_FETCH_STATUS)
		else()
			if(NOT EXISTS "${COMPONENTS_CIF}.gz")
				file(
					DOWNLOAD
					https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz
					${COMPONENTS_CIF}.gz
					SHOW_PROGRESS
					STATUS CCD_FETCH_STATUS)
			endif()

			add_custom_command(
				OUTPUT ${COMPONENTS_CIF}
				COMMAND "${GUNZIP}" ${COMPONENTS_CIF}.gz
				WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/rsrc/)

			add_custom_target(COMPONENTS ALL DEPENDS ${COMPONENTS_CIF})
		endif()

		# Do not continue if downloading went wrong
		list(POP_FRONT CCD_FETCH_STATUS CCD_FETCH_STATUS_CODE)

		if(CCD_FETCH_STATUS_CODE)
			message(
				FATAL_ERROR "Error trying to download CCD file: ${CCD_FETCH_STATUS}")
		endif()
	endif()
endif()

# Installation directories
if(BUILD_FOR_CCP4)
	set(CIFPP_DATA_DIR
		"$ENV{CCP4}/share/libcifpp"
		CACHE PATH "Directory where dictionary and other static data is stored")
else()
	set(CIFPP_DATA_DIR
		"${CMAKE_INSTALL_FULL_DATADIR}/libcifpp"
		CACHE PATH "Directory where dictionary and other static data is stored")
endif()

if(CIFPP_DATA_DIR)
	target_compile_definitions(cifpp PUBLIC DATA_DIR="${CIFPP_DATA_DIR}")
endif()

if(UNIX AND NOT BUILD_FOR_CCP4)
	if("${CMAKE_INSTALL_PREFIX}" STREQUAL "/usr/local")
		set(CIFPP_CACHE_DIR
			"/var/cache/libcifpp"
			CACHE PATH "The directory where downloaded data files are stored")
	else()
		set(CIFPP_CACHE_DIR
			"${CMAKE_INSTALL_FULL_LOCALSTATEDIR}/cache/libcifpp"
			CACHE PATH "The directory where downloaded data files are stored")
	endif()

	target_compile_definitions(cifpp PUBLIC CACHE_DIR="${CIFPP_CACHE_DIR}")

	set(CIFPP_ETC_DIR
		"${CMAKE_INSTALL_FULL_SYSCONFDIR}"
		CACHE PATH "The directory where the update configuration file is stored")
else()
	unset(CIFPP_CACHE_DIR)
endif()

# Install rules
install(TARGETS cifpp
	EXPORT cifpp
	FILE_SET cifpp_headers DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(MSVC AND BUILD_SHARED_LIBS)
	install(
		FILES $<TARGET_PDB_FILE:cifpp>
		DESTINATION ${CMAKE_INSTALL_LIBDIR}
		OPTIONAL)
endif()

# Clean up old config files (with old names)
file(GLOB OLD_CONFIG_FILES
	${CMAKE_INSTALL_FULL_LIBDIR}/cmake/cifpp/cifppConfig*.cmake
	${CMAKE_INSTALL_FULL_LIBDIR}/cmake/cifpp/cifppTargets*.cmake)

if(OLD_CONFIG_FILES)
	message(
		STATUS "Installation will remove old config files: ${OLD_CONFIG_FILES}")
	install(CODE "file(REMOVE ${OLD_CONFIG_FILES})")
endif()

install(EXPORT cifpp
	NAMESPACE cifpp::
	FILE "cifpp-targets.cmake"
	DESTINATION lib/cmake/cifpp)

if(CIFPP_DATA_DIR AND CIFPP_DOWNLOAD_CCD)
	install(
		FILES ${CMAKE_CURRENT_SOURCE_DIR}/rsrc/mmcif_ddl.dic
		${CMAKE_CURRENT_SOURCE_DIR}/rsrc/mmcif_pdbx.dic
		${CMAKE_CURRENT_SOURCE_DIR}/rsrc/mmcif_ma.dic ${COMPONENTS_CIF}
		DESTINATION ${CIFPP_DATA_DIR})
endif()

if(CIFPP_CACHE_DIR AND CIFPP_DOWNLOAD_CCD)
	install(
		FILES ${CMAKE_CURRENT_SOURCE_DIR}/rsrc/mmcif_ddl.dic
		${CMAKE_CURRENT_SOURCE_DIR}/rsrc/mmcif_pdbx.dic
		${CMAKE_CURRENT_SOURCE_DIR}/rsrc/mmcif_ma.dic ${COMPONENTS_CIF}
		DESTINATION ${CIFPP_CACHE_DIR})
endif()

set(CONFIG_TEMPLATE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cifpp-config.cmake.in)

configure_package_config_file(
	${CONFIG_TEMPLATE_FILE} ${CMAKE_CURRENT_BINARY_DIR}/cifpp/cifpp-config.cmake
	INSTALL_DESTINATION lib/cmake/cifpp
	PATH_VARS CIFPP_DATA_DIR)

install(
	FILES "${CMAKE_CURRENT_BINARY_DIR}/cifpp/cifpp-config.cmake"
	"${CMAKE_CURRENT_BINARY_DIR}/cifpp/cifpp-config-version.cmake"
	DESTINATION lib/cmake/cifpp)

set_target_properties(
	cifpp
	PROPERTIES VERSION ${PROJECT_VERSION}
	SOVERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
	INTERFACE_cifpp_MAJOR_VERSION ${PROJECT_VERSION_MAJOR})

set_property(
	TARGET cifpp
	APPEND
	PROPERTY COMPATIBLE_INTERFACE_STRING cifpp_MAJOR_VERSION)

write_basic_package_version_file(
	"${CMAKE_CURRENT_BINARY_DIR}/cifpp/cifpp-config-version.cmake"
	VERSION ${PROJECT_VERSION}
	COMPATIBILITY AnyNewerVersion)

# In case we're included as sub_directory:
if(NOT PROJECT_IS_TOP_LEVEL)
	set(CIFPP_SHARE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/rsrc PARENT_SCOPE)
endif()

if(BUILD_TESTING)
	add_subdirectory(test)
endif()

# Optionally install the update scripts for CCD and dictionary files
if(CIFPP_INSTALL_UPDATE_SCRIPT)
	if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux" OR ${CMAKE_SYSTEM_NAME} STREQUAL
		"GNU")
		if("${CMAKE_INSTALL_PREFIX}" STREQUAL "/usr/local")
			set(CIFPP_CRON_DIR
				"/etc/cron.weekly"
				CACHE PATH "The cron directory, for the update script")
		else()
			set(CIFPP_CRON_DIR
				"${CIFPP_ETC_DIR}/cron.weekly"
				CACHE PATH "The cron directory, for the update script")
		endif()
	elseif(${CMAKE_SYSTEM_NAME} STREQUAL "FreeBSD")
		set(CIFPP_CRON_DIR
			"${CIFPP_ETC_DIR}/periodic/weekly"
			CACHE PATH "The cron directory, for the update script")
	else()
		message(FATAL_ERROR "Don't know where to install the update script")
	endif()

	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tools/update-libcifpp-data.in
		update-libcifpp-data @ONLY)
	install(
		FILES ${CMAKE_CURRENT_BINARY_DIR}/update-libcifpp-data
		DESTINATION ${CIFPP_CRON_DIR}
		PERMISSIONS OWNER_EXECUTE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE
		WORLD_READ)

	install(DIRECTORY DESTINATION ${CIFPP_CACHE_DIR})

	# a config file, to make it complete
	if(NOT EXISTS "${CIFPP_ETC_DIR}/libcifpp.conf")
		file(
			WRITE ${CMAKE_CURRENT_BINARY_DIR}/libcifpp.conf
			[[# Uncomment the next line to enable automatic updates
# update=true
]])
		install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libcifpp.conf
			DESTINATION "${CIFPP_ETC_DIR}")
		install(
			CODE "message(\"A configuration file has been written to ${CIFPP_ETC_DIR}/libcifpp.conf, please edit this file to enable automatic updates\")"
		)

		install(DIRECTORY DESTINATION "${CIFPP_ETC_DIR}/libcifpp/cache-update.d")
	endif()

	target_compile_definitions(cifpp PUBLIC CACHE_DIR="${CIFPP_CACHE_DIR}")
endif()

if(BUILD_DOCUMENTATION)
	add_subdirectory(docs)
endif()

set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_SOURCE_TGZ ON)
set(CPACK_SOURCE_TBZ2 OFF)
set(CPACK_SOURCE_TXZ OFF)
set(CPACK_SOURCE_TZ OFF)
set(CPACK_SOURCE_IGNORE_FILES "/rsrc/components.cif;/build;/.vscode;/.git")
set(CPACK_PACKAGE_FILE_NAME "${PROJECT_NAME}-${PROJECT_VERSION}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME ${CPACK_PACKAGE_FILE_NAME})
include(CPack)
