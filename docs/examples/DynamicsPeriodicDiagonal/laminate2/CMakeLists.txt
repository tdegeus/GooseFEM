cmake_minimum_required(VERSION 3.1)

# define a project name
project(example)

# set optimization level
set(CMAKE_BUILD_TYPE Release)
add_definitions(-DNDEBUG)

# set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# load pkg-config
find_package(PkgConfig)

# find cppmat
pkg_check_modules(CPPMAT REQUIRED cppmat)
include_directories(${CPPMAT_INCLUDE_DIRS})

# find GooseFEM
pkg_check_modules(GOOSEFEM REQUIRED GooseFEM)
include_directories(${GOOSEFEM_INCLUDE_DIRS})

# find GooseMaterial
pkg_check_modules(GOOSEMATERIAL REQUIRED GooseMaterial)
include_directories(${GOOSEMATERIAL_INCLUDE_DIRS})

# find eigen3
pkg_check_modules(EIGEN3 REQUIRED eigen3)
include_directories(${EIGEN3_INCLUDE_DIRS})

# find HDF5
find_package(HDF5 COMPONENTS CXX REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
set(LINK_LIBS ${HDF5_LIBS} ${HDF5_LIBRARIES})

# find HDF5pp
pkg_check_modules(HDF5PP REQUIRED HDF5pp)
include_directories(${HDF5PP_INCLUDE_DIRS})

# add executable
add_executable(${PROJECT_NAME} main.cpp)

# link libraries
target_link_libraries(${PROJECT_NAME} ${LINK_LIBS})
