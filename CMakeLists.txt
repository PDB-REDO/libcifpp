cmake_minimum_required(VERSION 3.10)

# set the project name
project(cifpp VERSION 1.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set bindir, if not use -DBIN_INSTALL_DIR
if(NOT BIN_INSTALL_DIR)
  if(CMAKE_INSTALL_BINDIR)
    set(BIN_INSTALL_DIR ${CMAKE_INSTALL_BINDIR})
  else(CMAKE_INSTALL_BINDIR)
    set(BIN_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/bin")
  endif(CMAKE_INSTALL_BINDIR)
endif(NOT BIN_INSTALL_DIR)

# Set libdir, if not use -DLIB_INSTALL_DIR
if(NOT LIB_INSTALL_DIR)
  if(CMAKE_INSTALL_LIBDIR)
    set(LIB_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR})
  else(CMAKE_INSTALL_LIBDIR)
    set(LIB_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}")
  endif(CMAKE_INSTALL_LIBDIR)
endif(NOT LIB_INSTALL_DIR)

# Set includedir, if not use -DINCLUDE_INSTALL_DIR
if(NOT INCLUDE_INSTALL_DIR)
  if(CMAKE_INSTALL_INCLUDEDIR)
    set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR})
  else(CMAKE_INSTALL_INCLUDEDIR)
    set(INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include")
  endif(CMAKE_INSTALL_INCLUDEDIR)
endif(NOT INCLUDE_INSTALL_DIR)

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

add_library(${PROJECT_NAME} ${project_sources} ${project_headers})
add_executable(unit-test ${PROJECT_SOURCE_DIR}/test/unit-test.cpp)

configure_file(libcifpp.pc.in ${CMAKE_BINARY_DIR}/libcifpp.pc IMMEDIATE @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/libcifpp.pc DESTINATION ${LIB_INSTALL_DIR}/pkgconfig)

configure_file(${PROJECT_SOURCE_DIR}/src/Config.hpp.in Config.hpp)
target_include_directories(${PROJECT_NAME} PUBLIC "${PROJECT_BINARY_DIR}")
# install(FILES ${CMAKE_BINARY_DIR}/Config.hpp DESTINATION ${INCLUDE_INSTALL_DIR}/${PROJECT_NAME})

if(Boost_FOUND)
    include_directories(${PROJECT_NAME} PUBLIC ${Boost_INCLUDE_DIRS})
    link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
endif()

install(DIRECTORY include/cif++ DESTINATION ${INCLUDE_INSTALL_DIR})

install(FILES
	${PROJECT_SOURCE_DIR}/rsrc/mmcif_ddl.dic
	${PROJECT_SOURCE_DIR}/rsrc/mmcif_pdbx_v50.dic
	DESTINATION ${SHARE_INSTALL_DIR}/libcifpp/
)

install(TARGETS ${PROJECT_NAME}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}
    RUNTIME DESTINATION ${BIN_INSTALL_DIR}
)

target_compile_definitions(${PROJECT_NAME}
	PUBLIC "" CACHE_DIR="${SHARE_INSTALL_DIR}/libcifpp" DATA_DIR="${SHARE_INSTALL_DIR}/libcifpp" )

target_compile_definitions(unit-test
	PUBLIC "" CACHE_DIR="${SHARE_INSTALL_DIR}/libcifpp" DATA_DIR="${SHARE_INSTALL_DIR}/libcifpp" )

target_link_libraries(unit-test ${PROJECT_NAME} ${Boost_LIBRARIES})

if(MSVC)
    # make msvc standards compliant...
    target_compile_options(${PROJECT_NAME} PRIVATE /permissive-)
    target_compile_options(unit-test PRIVATE /permissive-)
endif()

