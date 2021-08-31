cmake_minimum_required(VERSION 3.16)

# set the project name
project(cifpp VERSION 1.1.0 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

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

if(CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter")
endif()
if(MSVC)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /W4")
endif()

# Build shared libraries by default (not my cup of tea, but hey)
option(BUILD_SHARED_LIBS "Build a shared library instead of a static one" OFF)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

# Only try to recreate SymOpTable_data.hpp if CCP4 is known
option(CCP4 "The location where ccp4 is installed" "")
if(DEFINED CCP4 AND EXISTS ${CCP4})
	set(CLIBD ${CCP4}/lib/data)
elseif(EXISTS $ENV{CCP4})
	set(CCP4 $ENV{CCP4})
	set(CLIBD ${CCP4}/lib/data)
endif()

if(EXISTS "${CCP4}")
	option(RECREATE_SYMOP_DATA "Recreate SymOp data table in case it is out of date" ON)

	if(RECREATE_SYMOP_DATA AND NOT EXISTS "${CLIBD}/syminfo.lib")
		message(FATAL_ERROR "Symop data table recreation requested, but file syminfo.lib was not found in ${CLIBD}")
	endif()
else()
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
	if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
		message(WARNING "The library and auxiliary files will be installed in $ENV{LOCALAPPDATA}/${PROJECT_NAME}")
		set(CMAKE_INSTALL_PREFIX "$ENV{LOCALAPPDATA}/${PROJECT_NAME}" CACHE PATH "..." FORCE)
	endif()

	# Find out the processor type for the target
	if(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "AMD64")
		set(COFF_TYPE "x64")
	elseif(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "i386")
		set(COFF_TYPE "x86")
	elseif(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "ARM64")
		set(COFF_TYPE "arm64")
	else()
		message(FATAL_ERROR "Unsupported or unknown processor type ${CMAKE_SYSTEM_PROCESSOR}")
	endif()	

	set(COFF_SPEC "--coff=${COFF_TYPE}")
endif()

if(UNIX AND NOT APPLE)
	# On Linux, install in the $HOME/.local folder by default
	if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
		message(WARNING "The library and auxiliary files will be installed in $ENV{HOME}/.local")
		set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local" CACHE PATH "..." FORCE)
	endif()
endif()

# Optionally use mrc to create resources
find_program(MRC mrc HINTS "$ENV{LOCALAPPDATA}/mrc" "$ENV{HOME}/.local/bin" "${CMAKE_INSTALL_PREFIX}/../mrc" "/usr/local/bin")

if(MRC)
	option(USE_RSRC "Use mrc to create resources" ON)
else()
	message(WARNING "Not using resources since mrc was not found")
endif()

if(USE_RSRC STREQUAL "ON")
	set(USE_RSRC 1)

	message("Using resources compiled with ${MRC}")
	add_compile_definitions(USE_RSRC)
endif()

# Libraries

set(CMAKE_THREAD_PREFER_PTHREAD)
set(THREADS_PREFER_PTHREAD_FLAG)
find_package(Threads)

set (Boost_DETAILED_FAILURE_MSG ON)
find_package(Boost 1.71.0 REQUIRED COMPONENTS system iostreams regex date_time program_options)

find_package(ZLIB)
find_package(BZip2)

include_directories(${Boost_INCLUDE_DIR})
link_libraries(${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

# Create a revision file, containing the current git version info

find_package(Git)
if(GIT_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
	include(GetGitRevisionDescription)
    get_git_head_revision(REFSPEC COMMITHASH)

	# Generate our own version string
	git_describe_working_tree(BUILD_VERSION_STRING --match=build --dirty)
else()
    message(WARNING "Git not found, cannot set version info")

    SET(BUILD_VERSION_STRING "unknown")
endif()

# generate version.h
include_directories(${CMAKE_BINARY_DIR} PRIVATE)
string(TIMESTAMP BUILD_DATE_TIME "%Y-%m-%dT%H:%M:%SZ" UTC)
configure_file("${CMAKE_SOURCE_DIR}/src/revision.hpp.in" "${CMAKE_BINARY_DIR}/revision.hpp" @ONLY)

# SymOp data table
if(RECREATE_SYMOP_DATA)
	# The tool to create the table

	add_executable(symop-map-generator "${CMAKE_SOURCE_DIR}/tools/symop-map-generator.cpp")

	target_link_libraries(symop-map-generator Threads::Threads ${Boost_LIBRARIES})

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

add_library(cifpp ${project_sources} ${project_headers} ${PROJECT_SOURCE_DIR}/src/SymOpTable_data.hpp)
set_target_properties(cifpp PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(cifpp
	PUBLIC
	"$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
	"$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_link_options(cifpp PRIVATE -undefined dynamic_lookup)
endif (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")

# download the components.cif file from CCD
set(COMPONENTS_CIF ${PROJECT_SOURCE_DIR}/data/components.cif)

if (NOT EXISTS ${COMPONENTS_CIF})

	if (NOT EXISTS ${PROJECT_SOURCE_DIR}/data)
		file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/data/)
	endif()

	if(${CMAKE_VERSION} VERSION_LESS "3.18.0")
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
	else()
		file(DOWNLOAD ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif.gz ${COMPONENTS_CIF}.gz
			SHOW_PROGRESS)
		file(ARCHIVE_EXTRACT INPUT ${COMPONENTS_CIF}.gz
			DESTINATION ${CMAKE_SOURCE_DIR}/data/
			VERBOSE)
	endif()
endif()

add_custom_target(COMPONENTS ALL DEPENDS ${COMPONENTS_CIF})

include_directories(${Boost_INCLUDE_DIRS})
link_libraries(${Boost_LIBRARIES})

set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR} )
set(LIBRARY_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR} )
set(SHARE_INSTALL_DIR ${CMAKE_INSTALL_DATADIR}/libcifpp)

# Fix cache dir
add_compile_definitions(
	# CACHE_DIR="${CMAKE_INSTALL_PREFIX}/${SHARE_INSTALL_DIR}"
	DATA_DIR="${CMAKE_INSTALL_PREFIX}/${SHARE_INSTALL_DIR}" )

generate_export_header(cifpp
	EXPORT_FILE_NAME cif++/Cif++Export.hpp)

# Install rules

install(TARGETS cifpp
	EXPORT cifppTargets
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

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
	${PROJECT_SOURCE_DIR}/data/components.cif
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

set(cifpp_MAJOR_VERSION 1)
set_target_properties(cifpp PROPERTIES
	VERSION ${PROJECT_VERSION}
	SOVERSION 1
	INTERFACE_cifpp_MAJOR_VERSION 1)

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

	if(USE_RSRC)
		add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/cifpp_test_rsrc.obj
			COMMAND ${MRC} -o ${CMAKE_CURRENT_BINARY_DIR}/cifpp_test_rsrc.obj ${CMAKE_SOURCE_DIR}/rsrc/mmcif_pdbx_v50.dic ${COFF_SPEC}
		)
		set(CIFPP_TEST_RESOURCE ${CMAKE_CURRENT_BINARY_DIR}/cifpp_test_rsrc.obj)
	endif()

	list(APPEND CIFPP_tests
		# pdb2cif
		rename-compound
		unit)

	foreach(CIFPP_TEST IN LISTS CIFPP_tests)
		set(CIFPP_TEST "${CIFPP_TEST}-test")
		set(CIFPP_TEST_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/test/${CIFPP_TEST}.cpp")

		add_executable(${CIFPP_TEST} ${CIFPP_TEST_SOURCE} ${CIFPP_TEST_RESOURCE})

		target_include_directories(${CIFPP_TEST} PRIVATE
			${CMAKE_CURRENT_SOURCE_DIR}/include
			${CMAKE_CURRENT_BINARY_DIR}  # for config.h
		)

		target_link_libraries(${CIFPP_TEST} Threads::Threads ${Boost_LIBRARIES} cifpp)

		if(${ZLIB_FOUND})
			target_link_libraries(${CIFPP_TEST} ZLIB::ZLIB)
		endif()

		if(${BZip2_FOUND})
			target_link_libraries(${CIFPP_TEST} BZip2::BZip2)
		endif()
	
		if(MSVC)
			# Specify unwind semantics so that MSVC knowns how to handle exceptions
			target_compile_options(${CIFPP_TEST} PRIVATE /EHsc)
		endif()

		add_custom_target("run-${CIFPP_TEST}" DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/Run${CIFPP_TEST}.touch ${CIFPP_TEST})

		add_custom_command(
			OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/Run${CIFPP_TEST}.touch
			COMMAND $<TARGET_FILE:${CIFPP_TEST}>
			WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test)

		add_test(NAME ${CIFPP_TEST}
			COMMAND $<TARGET_FILE:${CIFPP_TEST}>
			WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test)

	endforeach()
endif()