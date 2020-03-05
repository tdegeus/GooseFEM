cmake_minimum_required(VERSION 3.1)

project(example)

set(HIGHFIVE_USE_XTENSOR 1)
set(HIGHFIVE_USE_BOOST 0)
find_package(HighFive REQUIRED)
find_package(GooseFEM REQUIRED)
find_package(GMatElastic REQUIRED)

add_executable(example example.cpp)
target_link_libraries(example PRIVATE GMatElastic GooseFEM HighFive)
target_link_libraries(example INTERFACE xtensor::optimize xtensor::use_xsimd)
