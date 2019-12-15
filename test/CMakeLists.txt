
cmake_minimum_required(VERSION 3.0)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    project(GooseFEM-test)
    find_package(GooseFEM REQUIRED CONFIG)
endif()

option(WARNINGS "Show build warnings" ON)
option(SIMD "Enable xsimd" ON)
option(DEBUG "Enable all assertions" ON)

find_package(Catch2 REQUIRED)

add_executable(main main.cpp)

target_link_libraries(main PRIVATE Catch2::Catch2 GooseFEM)

if (SIMD)
    target_link_libraries(main PRIVATE xtensor::optimize xtensor::use_xsimd)
endif()

if (WARNINGS)
    target_link_libraries(main PRIVATE GooseFEM::compiler_warnings)
endif()

if (DEBUG)
    target_link_libraries(main PRIVATE GooseFEM::debug)
endif()
