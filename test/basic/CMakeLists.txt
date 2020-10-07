
cmake_minimum_required(VERSION 3.0)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    project(GooseFEM-test-basic)
    find_package(GooseFEM REQUIRED CONFIG)
endif()

option(WARNINGS "Show build warnings" ON)
option(SIMD "Enable xsimd" ON)
option(ASSERT "Enable assertions" ON)
option(DEBUG "Enable all assertions" ON)

find_package(Catch2 REQUIRED)

add_executable(test-basic
    main.cpp
    Allocate.cpp
    ElementHex8.cpp
    ElementQuad4.cpp
    Iterate.cpp
    Matrix.cpp
    MatrixDiagonal.cpp
    Mesh.cpp
    MeshQuad4.cpp
    Vector.cpp
    VectorPartitioned.cpp)

target_link_libraries(test-basic PRIVATE Catch2::Catch2 GooseFEM)

if (SIMD)
    target_link_libraries(test-basic PRIVATE xtensor::optimize xtensor::use_xsimd)
endif()

if (WARNINGS)
    target_link_libraries(test-basic PRIVATE GooseFEM::compiler_warnings)
endif()

if (ASSERT)
    target_link_libraries(test-basic PRIVATE GooseFEM::assert)
endif()

if (DEBUG)
    target_link_libraries(test-basic PRIVATE GooseFEM::debug)
endif()
