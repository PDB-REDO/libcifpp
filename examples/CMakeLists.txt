cmake_minimum_required(VERSION 3.15)
project(cifpp_example LANGUAGES CXX)

find_package(cifpp REQUIRED)

add_executable(example example.cpp)
target_link_libraries(example cifpp::cifpp)