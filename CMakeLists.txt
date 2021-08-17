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

set(PACKAGE ${PROJECT_NAME})
set(PACKAGE_VERSION ${PROJECT_VERSION})
add_definitions(-DVERSION=${PACKAGE_VERSION})

set(CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter")
endif()
if(MSVC)
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /W4")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
endif()

set(CMAKE_DEBUG_POSTFIX d)

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

# Optionally use mrc to create resources
find_program(MRC mrc HINTS "$ENV{LOCALAPPDATA}/mrc" "${CMAKE_INSTALL_PREFIX}/../mrc" "/usr/local/bin")

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

include_directories(${Boost_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include)
link_libraries(${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

# Create a revision file, containing the current git version info

find_package(Git)
if(GIT_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
    include(GetGitRevisionDescription)

	execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD OUTPUT_VARIABLE SHORT_SHA OUTPUT_STRIP_TRAILING_WHITESPACE)

    set(REVISION ${SHORT_SHA} CACHE STRING "git short sha" FORCE)

    # only use the plugin to tie the configure state to the sha to force rebuilds
    # of files that depend on version.h

	include(GetGitRevisionDescription)
    get_git_head_revision(REFSPEC COMMITHASH)

	# Generate our own version string
	git_describe_working_tree(BUILD_VERSION_STRING --match=build --dirty)
else()
    message(WARNING "Git not found, cannot set version info")

    SET(BUILD_VERSION_STRING "unknown")
endif()

# generate version.h
include_directories(${CMAKE_BINARY_DIR})
configure_file("${CMAKE_SOURCE_DIR}/src/revision.hpp.in" "${CMAKE_BINARY_DIR}/revision.hpp" @ONLY)
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

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
)

add_library(cifpp_obj OBJECT ${project_sources} ${project_headers})
set_target_properties(cifpp_obj PROPERTIES POSITION_INDEPENDENT_CODE 1)

add_library(cifpp SHARED $<TARGET_OBJECTS:cifpp_obj>)
set_target_properties(cifpp PROPERTIES
	SOVERSION ${PROJECT_VERSION}
	INTERFACE_cifpp_MAJOR_VERSION cifpp_MAJOR_VERSION
	COMPATIBLE_INTERFACE_STRING cifpp_MAJOR_VERSION
)

add_library(cifpp_static STATIC $<TARGET_OBJECTS:cifpp_obj>)

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
add_compile_definitions(CACHE_DIR="${SHARE_INSTALL_DIR}/libcifpp" DATA_DIR="${SHARE_INSTALL_DIR}/libcifpp" )

configure_file(${PROJECT_SOURCE_DIR}/src/Config-cmake.hpp.in Config.hpp)

generate_export_header(cifpp)

# Install rules

install(TARGETS cifpp cifpp_static
	EXPORT cifppTargets
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})


install(
	DIRECTORY include/cif++
	DESTINATION include
	COMPONENT Devel
)

install(
	FILES "${CMAKE_CURRENT_BINARY_DIR}/cifpp_export.h"
	DESTINATION include/cif++
	COMPONENT Devel
)

install(FILES
	${PROJECT_SOURCE_DIR}/rsrc/mmcif_ddl.dic
	${PROJECT_SOURCE_DIR}/rsrc/mmcif_pdbx_v50.dic
	${PROJECT_SOURCE_DIR}/data/components.cif
	DESTINATION ${CMAKE_INSTALL_DATADIR}
)

set(INCLUDE_INSTALL_DIR include )
set(LIBRARY_INSTALL_DIR lib )
set(SHARE_INSTALL_DIR share/libcifpp)

set(ConfigPackageLocation lib/cmake/cifpp)

configure_package_config_file(cmake/cifppConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cifpp/cifppConfig.cmake
  INSTALL_DESTINATION ${ConfigPackageLocation}
  PATH_VARS INCLUDE_INSTALL_DIR LIBRARY_INSTALL_DIR SHARE_INSTALL_DIR)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/cifpp/cifppConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)

export(EXPORT cifppTargets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/cifpp/cifppTargets.cmake"
  NAMESPACE cifpp::
)

install(EXPORT cifppTargets
  FILE cifppTargets.cmake
  NAMESPACE cifpp::
  DESTINATION ${ConfigPackageLocation}
  COMPONENT Devel
)

install(
  FILES
	"${CMAKE_CURRENT_BINARY_DIR}/cifpp/cifppConfig.cmake"
	"${CMAKE_CURRENT_BINARY_DIR}/cifpp/cifppConfigVersion.cmake"
  DESTINATION ${ConfigPackageLocation}
  COMPONENT Devel
)

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

		target_link_libraries(${CIFPP_TEST} cifpp_static Threads::Threads ${Boost_LIBRARIES})

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
			COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${CIFPP_TEST}
		)

		add_test(NAME ${CIFPP_TEST}
			COMMAND $<TARGET_FILE:${CIFPP_TEST}>
			WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test)

	endforeach()
endif()

