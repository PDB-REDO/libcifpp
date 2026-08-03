# SPDX-License-Identifier: BSD-2-Clause
# 
# Copyright (c) 2025  NKI/AVL, Netherlands Cancer Institute
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

if(NOT (Catch2_FOUND OR TARGET Catch2))
	find_package(Catch2 3 QUIET)

	if(NOT Catch2_FOUND)
		FetchContent_Declare(
			Catch2
			GIT_REPOSITORY https://github.com/catchorg/Catch2.git
			GIT_TAG v3.4.0)

		FetchContent_MakeAvailable(Catch2)

		target_compile_features(Catch2 PRIVATE cxx_std_20)
	endif()
endif()

list(
	APPEND
	CIFPP_tests
	unit-v2
	unit-3d
	model
	# query
	rename-compound
	sugar
	spinner
	reconstruction
	validate-pdbx
	cql
	matrix
)

add_library(test-main OBJECT "${CMAKE_CURRENT_SOURCE_DIR}/test-main.cpp")

target_link_libraries(test-main cifpp::cifpp Catch2::Catch2)

foreach(CIFPP_TEST IN LISTS CIFPP_tests)
	set(CIFPP_TEST "${CIFPP_TEST}-test")
	set(CIFPP_TEST_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/${CIFPP_TEST}.cpp")

	add_executable(
		${CIFPP_TEST} ${CIFPP_TEST_SOURCE} $<TARGET_OBJECTS:test-main>)

	target_link_libraries(${CIFPP_TEST} PRIVATE cifpp::cifpp Catch2::Catch2)
	target_include_directories(${CIFPP_TEST} PRIVATE "${EIGEN_INCLUDE_DIR}")
	target_compile_features(${CIFPP_TEST} PUBLIC cxx_std_23)

	if(MSVC)
		# Specify unwind semantics so that MSVC knowns how to handle exceptions
		target_compile_options(${CIFPP_TEST} PRIVATE /EHsc)
	endif()

	if(NOT (CIFPP_TEST STREQUAL "spinner-test"))
		add_test(NAME ${CIFPP_TEST}
			COMMAND $<TARGET_FILE:${CIFPP_TEST}> --data-dir ${CMAKE_CURRENT_SOURCE_DIR})
	endif()
endforeach()
