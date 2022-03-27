cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

project(1_finite_field LANGUAGES CXX)

include(CTest)

set (CMAKE_CXX_STANDARD 17)

add_subdirectory(lib)
add_subdirectory(app)
