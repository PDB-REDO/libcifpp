# SPDX-License-Identifier: BSD-2-Clause

# Copyright (c) 2021 NKI/AVL, Netherlands Cancer Institute

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.

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
project(cifpp VERSION 4.2.1 LANGUAGES CXX)

list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

enable_testing()

include(GNUInstallDirs)
include(CheckFunctionExists)
include(CheckIncludeFiles)
include(CheckLibraryExists)
include(CMakePackageConfigHelpers)
include(Dart)
include(GenerateExportHeader)

set(CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
if(BUILD_FOR_CCP4)
	if("$ENV{CCP4}" STREQUAL "" OR NOT EXISTS $ENV{CCP4})
		message(FATAL_ERROR "A CCP4 built was requested but CCP4 was not sourced")
	else()
		list(APPEND CMAKE_MODULE_PATH "$ENV{CCP4}")
		list(APPEND CMAKE_PREFIX_PATH "$ENV{CCP4}")
		set(CMAKE_INSTALL_PREFIX "$ENV{CCP4}")

		# This is the only option:
		if(WIN32)
			set(BUILD_SHARED_LIBS ON)
		endif()
	endif("$ENV{CCP4}" STREQUAL "" OR NOT EXISTS $ENV{CCP4})

	add_definitions(-DCCP4=1)
endif()

# Check if CCP4 is available
if(EXISTS "$ENV{CCP4}")
	set(CCP4 $ENV{CCP4})
	set(CLIBD ${CCP4}/lib/data)
endif()
if(CCP4 AND NOT CLIBD)
	set(CLIBD ${CCP4}/lib/data)
endif()

# When CCP4 is sourced in the environment, we can recreate the symmetry operations table
if(EXISTS "${CCP4}")
	if(CIFPP_RECREATE_SYMOP_DATA AND NOT EXISTS "${CLIBD}/syminfo.lib")
		message(WARNING "Symop data table recreation requested, but file syminfo.lib was not found in ${CLIBD}")
		set(CIFPP_RECREATE_SYMOP_DATA OFF)
	else()
		option(CIFPP_RECREATE_SYMOP_DATA "Recreate SymOp data table in case it is out of date" ON)
	endif()
else()
	set(CIFPP_RECREATE_SYMOP_DATA OFF)
	message("Not trying to recreate SymOpTable_data.hpp since CCP4 is not defined")
endif()

# set(CMAKE_DEBUG_POSTFIX d)

if(MSVC)
    # make msvc standards compliant...
    add_compile_options(/permissive-)

	macro(get_WIN32_WINNT version)
		if (WIN32 AND CMAKE_SYSTEM_VERSION)
			set(ver ${CMAKE_SYSTEM_VERSION})
			string(REPLACE "." "" ver ${ver})
			string(REGEX REPLACE "([0-9])" "0\\1" ver ${ver})

			set(${version} "0x${ver}")
		endif()
	endmacro()

	get_WIN32_WINNT(ver)
	add_definitions(-D_WIN32_WINNT=${ver})

	# On Windows, do not install in the system location
	if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT AND NOT BUILD_FOR_CCP4)
		message(STATUS "The library and auxiliary files will be installed in $ENV{LOCALAPPDATA}/${PROJECT_NAME}")
		set(CMAKE_INSTALL_PREFIX "$ENV{LOCALAPPDATA}/${PROJECT_NAME}" CACHE PATH "..." FORCE)
	endif()
endif()

if(UNIX AND NOT APPLE AND NOT BUILD_FOR_CCP4 AND CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	# On Linux, install in the $HOME/.local folder by default
	message(STATUS "The library and auxiliary files will be installed in $ENV{HOME}/.local")
	set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local" CACHE PATH "..." FORCE)
endif()

# Optionally use mrc to create resources

if(WIN32 AND BUILD_SHARED_LIBS)
	message("Not using resources when building shared libraries for Windows")
else()
	find_package(Mrc)

	if(MRC_FOUND)
		option(USE_RSRC "Use mrc to create resources" ON)
	else()
		message(WARNING "Not using resources since mrc was not found")
	endif()

	if(CIFPP_USE_RSRC STREQUAL "ON")
		set(CIFPP_USE_RSRC 1)

		message("Using resources compiled with ${MRC}")
		add_compile_definitions(CIFPP_USE_RSRC)
	endif()
endif()

# Libraries

set(CMAKE_THREAD_PREFER_PTHREAD)
set(THREADS_PREFER_PTHREAD_FLAG)
find_package(Threads)

find_package(Boost 1.70.0 REQUIRED COMPONENTS system iostreams regex program_options)

if(NOT MSVC AND Boost_USE_STATIC_LIBS)
	find_package(ZLIB REQUIRED)
	list(APPEND CIFPP_REQUIRED_LIBRARIES ZLIB::ZLIB)
endif()

include(FindFilesystem)
list(APPEND CIFPP_REQUIRED_LIBRARIES ${STDCPPFS_LIBRARY})

include(FindAtomic)
list(APPEND CIFPP_REQUIRED_LIBRARIES ${STDCPPATOMIC_LIBRARY})

# Create a revision file, containing the current git version info
include(VersionString)
write_version_header("LibCIFPP")

# SymOp data table
if(CIFPP_RECREATE_SYMOP_DATA)
	# The tool to create the table

	add_executable(symop-map-generator "${CMAKE_SOURCE_DIR}/tools/symop-map-generator.cpp")

	target_link_libraries(symop-map-generator Threads::Threads ${Boost_LIBRARIES} ${CIFPP_REQUIRED_LIBRARIES})
	if(Boost_INCLUDE_DIR)
		target_include_directories(symop-map-generator PUBLIC ${Boost_INCLUDE_DIR})
	endif()

	set($ENV{CLIBD} ${CLIBD})
	
	add_custom_command(
		OUTPUT ${CMAKE_SOURCE_DIR}/src/SymOpTable_data.hpp
		COMMAND $<TARGET_FILE:symop-map-generator> ${CLIBD}/syminfo.lib ${CMAKE_SOURCE_DIR}/src/SymOpTable_data.hpp
		)
	
	add_custom_target(
		OUTPUT ${CMAKE_SOURCE_DIR}/src/SymOpTable_data.hpp
		DEPENDS symop-map-generator "$ENV{CLIBD}/syminfo.lib"
	)
endif()

# Sources

set(project_sources 
	${PROJECT_SOURCE_DIR}/src/AtomType.cpp
	${PROJECT_SOURCE_DIR}/src/BondMap.cpp
	${PROJECT_SOURCE_DIR}/src/Cif++.cpp
	${PROJECT_SOURCE_DIR}/src/Cif2PDB.cpp
	${PROJECT_SOURCE_DIR}/src/CifParser.cpp
	${PROJECT_SOURCE_DIR}/src/CifUtils.cpp
	${PROJECT_SOURCE_DIR}/src/CifValidator.cpp
	${PROJECT_SOURCE_DIR}/src/Compound.cpp
	${PROJECT_SOURCE_DIR}/src/PDB2Cif.cpp
	${PROJECT_SOURCE_DIR}/src/PDB2CifRemark3.cpp
	${PROJECT_SOURCE_DIR}/src/Point.cpp
	${PROJECT_SOURCE_DIR}/src/Secondary.cpp
	${PROJECT_SOURCE_DIR}/src/Structure.cpp
	${PROJECT_SOURCE_DIR}/src/Symmetry.cpp
	${PROJECT_SOURCE_DIR}/src/TlsParser.cpp
)

set(project_headers 
	${PROJECT_SOURCE_DIR}/include/cif++/AtomType.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/BondMap.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/Cif++.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/Cif2PDB.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/CifParser.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/CifUtils.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/CifValidator.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/Compound.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/PDB2Cif.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/PDB2CifRemark3.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/Point.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/Secondary.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/Structure.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/Symmetry.hpp
	${PROJECT_SOURCE_DIR}/include/cif++/TlsParser.hpp
)

add_library(cifpp ${project_sources} ${project_headers} ${CMAKE_SOURCE_DIR}/src/SymOpTable_data.hpp)
set_target_properties(cifpp PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(cifpp
	PUBLIC
	"$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
	"$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
	${Boost_INCLUDE_DIR}
)

target_include_directories(cifpp
	PRIVATE
	${CMAKE_BINARY_DIR}
)

target_link_libraries(cifpp PUBLIC Threads::Threads Boost::regex Boost::iostreams ${CIFPP_REQUIRED_LIBRARIES})
# target_link_libraries(cifpp PRIVATE)

if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_link_options(cifpp PRIVATE -undefined dynamic_lookup)
endif (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")

option(CIFPP_DOWNLOAD_CCD "Download the CCD file components.cif during installation" OFF)
if(CIFPP_DOWNLOAD_CCD)
	# download the components.cif file from CCD
	set(COMPONENTS_CIF ${PROJECT_SOURCE_DIR}/data/components.cif)

	if (NOT EXISTS ${COMPONENTS_CIF})

		if (NOT EXISTS ${PROJECT_SOURCE_DIR}/data)
			file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/data/)
		endif()

		find_program(GUNZIP gunzip)

		if(GUNZIP)
			file(DOWNLOAD ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif.gz ${COMPONENTS_CIF}.gz
				SHOW_PROGRESS)
			add_custom_command(OUTPUT ${COMPONENTS_CIF}
				COMMAND ${GUNZIP} ${COMPONENTS_CIF}.gz
				WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/data/)
		else()
			file(DOWNLOAD ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif ${COMPONENTS_CIF}
				SHOW_PROGRESS)
		endif()
	endif()

	add_custom_target(COMPONENTS ALL DEPENDS ${COMPONENTS_CIF})
endif()

if(UNIX)
	option(CIFPP_INSTALL_UPDATE_SCRIPT "Install the script to update CCD and dictionary files" OFF)
	set(CIFPP_CACHE_DIR "/var/cache/libcifpp" CACHE STRING "The cache directory to use")
	target_compile_definitions(cifpp PUBLIC CACHE_DIR="${CIFPP_CACHE_DIR}")
endif()

generate_export_header(cifpp
	EXPORT_FILE_NAME cif++/Cif++Export.hpp)

set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR} )
set(LIBRARY_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR} )
set(SHARE_INSTALL_DIR ${CMAKE_INSTALL_DATADIR}/libcifpp)

set(CIFPP_DATA_DIR "${CMAKE_INSTALL_PREFIX}/${SHARE_INSTALL_DIR}" CACHE STRING "The directory containing the provided data files")

target_compile_definitions(cifpp PUBLIC DATA_DIR="${CIFPP_DATA_DIR}")

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
	FILES "${CMAKE_CURRENT_BINARY_DIR}/cif++/Cif++Export.hpp"
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cif++
	COMPONENT Devel
)

install(FILES
	${PROJECT_SOURCE_DIR}/rsrc/mmcif_ddl.dic
	${PROJECT_SOURCE_DIR}/rsrc/mmcif_pdbx_v50.dic
	${COMPONENTS_CIF}
	DESTINATION ${SHARE_INSTALL_DIR}
)

configure_package_config_file(Config.cmake.in
	${CMAKE_CURRENT_BINARY_DIR}/cifpp/cifppConfig.cmake
	INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cifpp
	PATH_VARS INCLUDE_INSTALL_DIR LIBRARY_INSTALL_DIR SHARE_INSTALL_DIR
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

set(prefix      ${CMAKE_INSTALL_PREFIX})
set(exec_prefix ${CMAKE_INSTALL_PREFIX})
set(libdir      ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(includedir  ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/libcifpp.pc.in
	${CMAKE_CURRENT_BINARY_DIR}/libcifpp.pc.in @ONLY)
file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/libcifpp.pc
	INPUT ${CMAKE_CURRENT_BINARY_DIR}/libcifpp.pc.in)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libcifpp.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

# Unit tests

option(CIFPP_BUILD_TESTS "Build test exectuables" OFF)

if(CIFPP_BUILD_TESTS)

	list(APPEND CIFPP_tests
		# pdb2cif
		rename-compound
		structure
		sugar
		unit)

	foreach(CIFPP_TEST IN LISTS CIFPP_tests)
		set(CIFPP_TEST "${CIFPP_TEST}-test")
		set(CIFPP_TEST_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/test/${CIFPP_TEST}.cpp")

		add_executable(${CIFPP_TEST} ${CIFPP_TEST_SOURCE})

		target_include_directories(${CIFPP_TEST} PRIVATE
			${CMAKE_CURRENT_SOURCE_DIR}/include
			${CMAKE_CURRENT_BINARY_DIR}  # for config.h
		)

		target_link_libraries(${CIFPP_TEST} PRIVATE Threads::Threads cifpp )

		if(CIFPP_USE_RSRC)
			mrc_target_resources(${CIFPP_TEST} ${CMAKE_SOURCE_DIR}/rsrc/mmcif_pdbx_v50.dic)
		endif()

		if(MSVC)
			# Specify unwind semantics so that MSVC knowns how to handle exceptions
			target_compile_options(${CIFPP_TEST} PRIVATE /EHsc)
		endif()

		add_custom_target("run-${CIFPP_TEST}" DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/Run${CIFPP_TEST}.touch ${CIFPP_TEST})

		add_custom_command(
			OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/Run${CIFPP_TEST}.touch
			COMMAND $<TARGET_FILE:${CIFPP_TEST}> -- ${CMAKE_SOURCE_DIR}/test)

		add_test(NAME ${CIFPP_TEST}
			COMMAND $<TARGET_FILE:${CIFPP_TEST}> -- ${CMAKE_SOURCE_DIR}/test)

	endforeach()
endif()

message("Will install in ${CMAKE_INSTALL_PREFIX}")

# Optionally install the update scripts for CCD and dictionary files

if(CIFPP_INSTALL_UPDATE_SCRIPT)
	set(CIFPP_CRON_DIR "$ENV{DESTDIR}/etc/cron.weekly")

	configure_file(${CMAKE_SOURCE_DIR}/tools/update-libcifpp-data.in update-libcifpp-data @ONLY)
	install(
		FILES ${CMAKE_CURRENT_BINARY_DIR}/update-libcifpp-data
		DESTINATION ${CIFPP_CRON_DIR}
		PERMISSIONS OWNER_EXECUTE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ
	)

	install(DIRECTORY DESTINATION ${CIFPP_CACHE_DIR})
	install(DIRECTORY DESTINATION "$ENV{DESTDIR}/etc/libcifpp/cache-update.d")

	# a config to, to make it complete
	if(NOT EXISTS "$ENV{DESTDIR}/etc/libcifpp.conf")
		file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/libcifpp.conf [[# Uncomment the next line to enable automatic updates
# update=true
]])
		install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libcifpp.conf DESTINATION "$ENV{DESTDIR}/etc")
		install(CODE "message(\"A configuration file has been written to $ENV{DESTDIR}/etc/libcifpp.conf, please edit this file to enable automatic updates\")")
	endif()

	target_compile_definitions(cifpp PUBLIC CACHE_DIR="${CIFPP_CACHE_DIR}")
endif()

