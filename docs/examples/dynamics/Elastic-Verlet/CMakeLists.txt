cmake_minimum_required(VERSION 3.1)

project(main)

set(HIGHFIVE_USE_XTENSOR 1)
set(HIGHFIVE_USE_BOOST 0)
find_package(HighFive REQUIRED)
find_package(GooseFEM REQUIRED)
find_package(GMatElastoPlasticQPot REQUIRED)

add_executable(main main.cpp)
target_link_libraries(main PRIVATE GMatElastoPlasticQPot GooseFEM HighFive)
target_link_libraries(main INTERFACE xtensor::optimize xtensor::use_xsimd)
