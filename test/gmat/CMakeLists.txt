
cmake_minimum_required(VERSION 3.0)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    project(GooseFEM-test-gmat)
    find_package(GooseFEM REQUIRED CONFIG)
endif()

option(WARNINGS "Show build warnings" ON)
option(SIMD "Enable xsimd" ON)
option(ASSERT "Enable assertions" ON)
option(DEBUG "Enable all assertions" ON)

find_package(Catch2 REQUIRED)

if(NOT WIN32)
    find_package(GMatElastic REQUIRED)

    add_executable(test-gmat
        main.cpp
        hybrid-elastic.cpp)

    target_link_libraries(test-gmat PRIVATE Catch2::Catch2 GooseFEM GMatElastic)
else()
    add_executable(test-gmat
        main.cpp)

    target_link_libraries(test-gmat PRIVATE Catch2::Catch2 GooseFEM)
endif()

if (SIMD)
    target_link_libraries(test-gmat PRIVATE xtensor::optimize xtensor::use_xsimd)
endif()

if (WARNINGS)
    target_link_libraries(test-gmat PRIVATE GooseFEM::compiler_warnings)
endif()

if (ASSERT)
    target_link_libraries(test-gmat PRIVATE GooseFEM::assert)
endif()

if (DEBUG)
    target_link_libraries(test-gmat PRIVATE GooseFEM::debug)
endif()
