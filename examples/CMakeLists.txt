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

cmake_minimum_required(VERSION 3.15)
project(cifpp_example LANGUAGES CXX)

if(PROJECT_IS_TOP_LEVEL AND NOT (TARGET cifpp OR cifpp_FOUND))
    find_package(cifpp REQUIRED)
endif()

add_executable(example example.cpp)
target_link_libraries(example cifpp::cifpp)

add_executable(example-cql example-cql.cpp)
target_link_libraries(example-cql cifpp::cifpp)

# Used in the smoke test in CI:

include(CTest)
add_test(NAME example_test COMMAND $<TARGET_FILE:example> ${CMAKE_CURRENT_SOURCE_DIR}/1cbs.cif.gz)
add_test(NAME example_test_cql COMMAND $<TARGET_FILE:example-cql> ${CMAKE_CURRENT_SOURCE_DIR}/1cbs.cif.gz)
