# Simple check to see if we need a library for std::atomic

if(TARGET std::atomic)
	return()
endif()

cmake_minimum_required(VERSION 3.10)

include(CMakePushCheckState)
include(CheckIncludeFileCXX)
include(CheckCXXSourceRuns)

cmake_push_check_state()

check_include_file_cxx("atomic" _CXX_ATOMIC_HAVE_HEADER)
mark_as_advanced(_CXX_ATOMIC_HAVE_HEADER)

set(code [[
#include <atomic>
int main(int argc, char** argv) {
  std::atomic<long long> s;
  ++s;
  return 0;
}
]])

check_cxx_source_runs("${code}" _CXX_ATOMIC_BUILTIN)

if(_CXX_ATOMIC_BUILTIN)
	set(_found 1)
else()
  list(APPEND CMAKE_REQUIRED_LIBRARIES atomic)
  list(APPEND FOLLY_LINK_LIBRARIES atomic)

  check_cxx_source_runs("${code}" _CXX_ATOMIC_LIB_NEEDED)
  if (NOT _CXX_ATOMIC_LIB_NEEDED)
    message(FATAL_ERROR "unable to link C++ std::atomic code: you may need \
      to install GNU libatomic")
  else()
	set(_found 1)
  endif()
endif()

if(_found)
	add_library(std::atomic INTERFACE IMPORTED)
	set_property(TARGET std::atomic APPEND PROPERTY INTERFACE_COMPILE_FEATURES cxx_std_14)

	if(_CXX_ATOMIC_BUILTIN)
		# Nothing to add...
	elseif(_CXX_ATOMIC_LIB_NEEDED)
		set_target_properties(std::atomic PROPERTIES IMPORTED_LIBNAME atomic)
		set(STDCPPATOMIC_LIBRARY atomic)
	endif()
endif()

cmake_pop_check_state()

set(Atomic_FOUND ${_found} CACHE BOOL "TRUE if we can run a program using std::atomic" FORCE)
mark_as_advanced(Atomic_FOUND)

if(Atomic_FIND_REQUIRED AND NOT Atomic_FOUND)
    message(FATAL_ERROR "Cannot run simple program using std::atomic")
endif()
