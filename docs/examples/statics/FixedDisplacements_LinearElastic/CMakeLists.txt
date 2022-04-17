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

# Add more executables used as illustration

add_executable(example_realloc example_realloc.cpp)
target_link_libraries(example_realloc PRIVATE GMatElastic GooseFEM HighFive)
target_link_libraries(example_realloc INTERFACE xtensor::optimize xtensor::use_xsimd)

add_executable(manual_partition manual_partition.cpp)
target_link_libraries(manual_partition PRIVATE GMatElastic GooseFEM HighFive)
target_link_libraries(manual_partition INTERFACE xtensor::optimize xtensor::use_xsimd)
