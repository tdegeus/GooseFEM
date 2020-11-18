
cmake_minimum_required(VERSION 3.0)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    project(GooseFEM-test-gmat)
    find_package(GooseFEM REQUIRED CONFIG)
endif()

set(ASSERT ON)
set(DEBUG ON)
option(SIMD "Enable xsimd" ON)
option(WARNINGS "Enable warnings" ON)

set(test_name "test-gmat")

find_package(Catch2 REQUIRED)

if(NOT WIN32)
    find_package(GMatElastic REQUIRED)
    find_package(GMatElastoPlasticQPot REQUIRED)

    add_executable(${test_name}
        main.cpp
        hybrid-elastic.cpp
        hybrid-elastoplasticqpot.cpp)

    target_link_libraries(${test_name} PRIVATE
        Catch2::Catch2
        GooseFEM
        GMatElastic
        GMatElastoPlasticQPot)
else()
    find_package(GMatElastoPlasticQPot REQUIRED)

    add_executable(${test_name}
        main.cpp
        hybrid-elastoplasticqpot.cpp)

    target_link_libraries(${test_name} PRIVATE
        Catch2::Catch2
        GooseFEM
        GMatElastoPlasticQPot)
endif()

if(SIMD)
    target_link_libraries(${test_name} PRIVATE xtensor::optimize xtensor::use_xsimd)
endif()

if(WARNINGS)
    target_link_libraries(${test_name} PRIVATE GooseFEM::compiler_warnings)
endif()

if(ASSERT)
    target_link_libraries(${test_name} PRIVATE GooseFEM::assert)
endif()

if(DEBUG)
    target_link_libraries(${test_name} PRIVATE GooseFEM::debug)
endif()

add_test(NAME ${test_name} COMMAND ${test_name})
