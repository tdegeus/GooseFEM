
cmake_minimum_required(VERSION 3.0)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    project(GooseFEM-test-basic)
    find_package(GooseFEM REQUIRED CONFIG)
endif()

set(ASSERT ON)
set(DEBUG ON)
option(XSIMD "Use xsimd optimisations" ON)

set(test_name "test-basic")

find_package(Catch2 REQUIRED)
find_package(xtensor REQUIRED)

add_executable(${test_name}
    main.cpp
    Allocate.cpp
    Element.cpp
    ElementHex8.cpp
    ElementQuad4.cpp
    Iterate.cpp
    Matrix.cpp
    MatrixDiagonal.cpp
    Mesh.cpp
    MeshQuad4.cpp
    Vector.cpp
    VectorPartitioned.cpp
    version.cpp)

target_link_libraries(${test_name} PRIVATE
    Catch2::Catch2
    GooseFEM
    GooseFEM::compiler_warnings
    xtensor::optimize)

if(ASSERT)
    target_link_libraries(${test_name} PRIVATE GooseFEM::assert)
endif()

if(DEBUG)
    target_link_libraries(${test_name} PRIVATE GooseFEM::debug)
endif()

if(SIMD)
    find_package(xsimd REQUIRED)
    target_link_libraries(${test_name} PRIVATE xtensor::use_xsimd)
endif()

add_test(NAME ${test_name} COMMAND ${test_name})
