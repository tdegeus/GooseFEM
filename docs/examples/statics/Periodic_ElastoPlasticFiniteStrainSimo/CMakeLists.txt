cmake_minimum_required(VERSION 3.1)

# basic info
# ----------

# define a project name
project(main)

# define empty list of libraries to link
set(PROJECT_LIBS "")

# copy files to build folder
configure_file("plot.py" "plot.py" COPYONLY)

# basic compiler options
# ----------------------

# set optimization level
set(CMAKE_BUILD_TYPE Release)

# set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set warnings on
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -g")

# enable optimisation
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DXTENSOR_USE_XSIMD=ON")

# assertions
if(NOT ASSERT)
  add_definitions(-DNDEBUG)
else()
  add_definitions(-DXTENSOR_ENABLE_ASSERT=ON)
  add_definitions(-DGOOSEFEM_ENABLE_ASSERT=ON)
  add_definitions(-DGMATELASTOPLASTICFINITESTRAINSIMO_ENABLE_ASSERT=ON)
endif()

# load packages
# -------------

# add custom include paths
if(NOT "$ENV{INCLUDE_PATH}" STREQUAL "")
  string(REPLACE ":" ";" INCLUDE_LIST "$ENV{INCLUDE_PATH}")
  include_directories(${INCLUDE_LIST})
endif()

# load HDF5
find_package(HDF5 COMPONENTS CXX REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
set(PROJECT_LIBS ${PROJECT_LIBS} ${HDF5_LIBS} ${HDF5_LIBRARIES})

# load eigen
find_package(PkgConfig)
pkg_check_modules(EIGEN3 REQUIRED eigen3)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIRS})

# create executable
# -----------------

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} ${PROJECT_LIBS})
