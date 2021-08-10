cmake_minimum_required(VERSION 3.10)

# set the project name
project(cifpp VERSION 1.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set sharedir, if not use -DSHARE_INSTALL_DIR
if(NOT SHARE_INSTALL_DIR)
  if(CMAKE_INSTALL_DATADIR)
    set(SHARE_INSTALL_DIR "${CMAKE_INSTALL_DATADIR}")
  else(CMAKE_INSTALL_DATADIR)
    set(SHARE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/share")
  endif(CMAKE_INSTALL_DATADIR)
endif(NOT SHARE_INSTALL_DIR)

set (Boost_DETAILED_FAILURE_MSG ON)
# set (BOOST_ROOT ${PROJECT_SOURCE_DIR}/../boost_1_75_0)
# set (Boost_COMPILER "-vc")
# set (Boost_USE_STATIC_RUNTIME ON)  

find_package(Boost 1.73.0 REQUIRED COMPONENTS system iostreams regex date_time)

set(CMAKE_THREAD_PREFER_PTHREAD)
set(THREADS_PREFER_PTHREAD_FLAG)
find_package(Threads)

find_package(ZLIB)
find_package(BZip2)

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
	${PROJECT_SOURCE_DIR}/src/SymOpTable_data.cpp
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

# download the components.cif file from CCD
if (NOT EXISTS ${PROJECT_SOURCE_DIR}/data)
	file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/data/)
endif()

# unfortunately, extract_archive does not recognize .gz files on Windows?
if (NOT EXISTS ${PROJECT_SOURCE_DIR}/data/components.cif)
	file(
		DOWNLOAD ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif ${PROJECT_SOURCE_DIR}/data/components.cif
		SHOW_PROGRESS)
endif()

add_library(${PROJECT_NAME} ${project_sources} ${project_headers})
add_executable(unit-test ${PROJECT_SOURCE_DIR}/test/unit-test.cpp)

configure_file(${PROJECT_SOURCE_DIR}/src/Config.hpp.in Config.hpp)
target_include_directories(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${Boost_INCLUDE_DIRS}>)

include_directories(${PROJECT_NAME} PUBLIC ${Boost_INCLUDE_DIRS})
link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})

target_compile_definitions(${PROJECT_NAME}
	PUBLIC "" CACHE_DIR="${SHARE_INSTALL_DIR}/libcifpp" DATA_DIR="${SHARE_INSTALL_DIR}/libcifpp" )

target_compile_definitions(unit-test
	PUBLIC "" CACHE_DIR="${SHARE_INSTALL_DIR}/libcifpp" DATA_DIR="${SHARE_INSTALL_DIR}/libcifpp" )

target_link_libraries(unit-test ${PROJECT_NAME} ${Boost_LIBRARIES} Threads::Threads)

if(${ZLIB_FOUND})
	target_link_libraries(unit-test ZLIB::ZLIB)
endif()

if(${BZip2_FOUND})
	target_link_libraries(unit-test BZip2::BZip2)
endif()

if(MSVC)
    # make msvc standards compliant...
    target_compile_options(${PROJECT_NAME} PRIVATE /permissive-)
    target_compile_options(unit-test PRIVATE /permissive-)
endif()

include(GenerateExportHeader)

generate_export_header(cifpp)
set_property(TARGET cifpp PROPERTY VERSION ${PROJECT_VERSION})
set_property(TARGET cifpp PROPERTY SOVERSION 1)
set_property(TARGET cifpp PROPERTY INTERFACE_cifpp_MAJOR_VERSION 1)
set_property(TARGET cifpp APPEND PROPERTY
  COMPATIBLE_INTERFACE_STRING cifpp_MAJOR_VERSION
)

install(TARGETS cifpp EXPORT cifppTargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
)

install(
	DIRECTORY include/cif++
	DESTINATION include
	COMPONENT Devel
)

install(
	FILES "${CMAKE_CURRENT_BINARY_DIR}/cifpp_export.h"
	DESTINATION include
	COMPONENT Devel
)

install(FILES
	${PROJECT_SOURCE_DIR}/rsrc/mmcif_ddl.dic
	${PROJECT_SOURCE_DIR}/rsrc/mmcif_pdbx_v50.dic
	${PROJECT_SOURCE_DIR}/data/components.cif
	DESTINATION share/libcifpp/
)

include(CMakePackageConfigHelpers)

set(INCLUDE_INSTALL_DIR include )
set(LIBRARY_INSTALL_DIR lib )

configure_package_config_file(cmake/cifppConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cifpp/cifppConfig.cmake
  INSTALL_DESTINATION lib/cmake/cifpp
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

set(ConfigPackageLocation lib/cmake/cifpp)
install(EXPORT cifppTargets
  FILE
    cifppTargets.cmake
  NAMESPACE
    cifpp::
  DESTINATION
    ${ConfigPackageLocation}
)

install(
  FILES
	"${CMAKE_CURRENT_BINARY_DIR}/cifpp/cifppConfig.cmake"
	"${CMAKE_CURRENT_BINARY_DIR}/cifpp/cifppConfigVersion.cmake"
  DESTINATION
    ${ConfigPackageLocation}
  COMPONENT
    Devel
)
