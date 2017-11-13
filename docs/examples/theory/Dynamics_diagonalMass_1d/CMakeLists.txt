cmake_minimum_required(VERSION 3.1)

# define a project name
project(main)

# set optimization level
set(CMAKE_BUILD_TYPE Release)
add_definitions(-DNDEBUG)

# set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# load pkg-config
find_package(PkgConfig)
# load header-only modules using pkg-config
pkg_check_modules(EIGEN3 REQUIRED eigen3)
include_directories(${EIGEN3_INCLUDE_DIRS})

# find pybind11
find_package(pybind11 REQUIRED)

# create the module
pybind11_add_module(main main.cpp)
