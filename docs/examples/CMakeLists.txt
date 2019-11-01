cmake_minimum_required(VERSION 3.1)

# define a project name
project(examples)

# set optimization level
set(CMAKE_BUILD_TYPE Release)

# set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set warnings on
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -g")

# enable optimisation
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DXTENSOR_USE_XSIMD=ON")

# assertions
if(NOT ASSERT)
  add_definitions(-DNDEBUG)
else()
  add_definitions(-DXTENSOR_ENABLE_ASSERT=ON)
endif()

# add current library
string(CONCAT GooseFEM "${CMAKE_CURRENT_SOURCE_DIR}" "/../../include/")
include_directories("${GooseFEM}")

# load xtensor-stack libraries
find_package(xtl REQUIRED)
find_package(xtensor REQUIRED)

# define empty list of libraries to link
set(PROJECT_LIBS "")

# load HDF5
find_package(HDF5 COMPONENTS CXX REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
set(PROJECT_LIBS ${PROJECT_LIBS} ${HDF5_LIBS} ${HDF5_LIBRARIES})

# load eigen
find_package(PkgConfig)
pkg_check_modules(EIGEN3 REQUIRED eigen3)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIRS})

# add custom include paths
if(NOT "$ENV{INCLUDE_PATH}" STREQUAL "")
  string(REPLACE ":" ";" INCLUDE_LIST "$ENV{INCLUDE_PATH}")
  include_directories(${INCLUDE_LIST})
endif()

# create executables

add_executable(
  "statics_FixedDisplacements_LinearElastic_example"
  "statics/FixedDisplacements_LinearElastic/example.cpp")

target_link_libraries(
  "statics_FixedDisplacements_LinearElastic_example"
  ${PROJECT_LIBS})

add_executable(
  "statics_FixedDisplacements_LinearElastic_manual_partition"
  "statics/FixedDisplacements_LinearElastic/manual_partition.cpp")

target_link_libraries(
  "statics_FixedDisplacements_LinearElastic_manual_partition"
  ${PROJECT_LIBS})

# copy Python scripts

configure_file(
  "statics/FixedDisplacements_LinearElastic/example.py"
  "statics_FixedDisplacements_LinearElastic_example.py" COPYONLY)

configure_file(
  "statics/FixedDisplacements_LinearElastic/manual_partition.py"
  "statics_FixedDisplacements_LinearElastic_manual_partition.py" COPYONLY)

configure_file(
  "statics/FixedDisplacements_LinearElastic/plot.py"
  "statics_FixedDisplacements_LinearElastic_plot.py" COPYONLY)
