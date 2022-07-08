cmake_minimum_required(VERSION 3.19..3.21)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    project(GooseFEM)
    find_package(GooseFEM REQUIRED CONFIG)
    option(USE_ASSERT "${PROJECT_NAME}: Build with assertions" ON)
    option(USE_DEBUG "${PROJECT_NAME}: Build in debug mode" OFF)
    option(USE_SIMD "${PROJECT_NAME}: Build with hardware optimization" OFF)
endif()

set(MYPROJECT "${PROJECT_NAME}-test-basic")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Catch2 REQUIRED)
find_package(xtensor REQUIRED)

add_library(mytarget INTERFACE IMPORTED)

target_link_libraries(mytarget INTERFACE
    ${PROJECT_NAME}
    ${PROJECT_NAME}::compiler_warnings
    Catch2::Catch2WithMain)

if(USE_ASSERT)
    target_link_libraries(mytarget INTERFACE ${PROJECT_NAME}::assert)
    message(STATUS "Compiling ${MYPROJECT} with assertions")
endif()

if(USE_DEBUG)
    target_link_libraries(mytarget INTERFACE ${PROJECT_NAME}::debug)
    message(STATUS "Compiling ${MYPROJECT} in debug mode")
endif()

if(USE_SIMD)
    find_package(xtensor REQUIRED)
    find_package(xsimd REQUIRED)
    target_link_libraries(mytarget INTERFACE xtensor::use_xsimd xtensor::optimize)
    message(STATUS "Compiling ${MYPROJECT} with hardware optimization")
endif()

file(GLOB APP_SOURCES *.cpp)

foreach(mysource ${APP_SOURCES})
    string(REPLACE ".cpp" "" myexec ${mysource})
    get_filename_component(myexec ${myexec} NAME)
    add_executable(${myexec} ${mysource})
    target_link_libraries(${myexec} PRIVATE mytarget)
    add_test(NAME ${myexec} COMMAND ${myexec})
endforeach()
