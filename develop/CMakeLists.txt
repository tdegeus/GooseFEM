
cmake_minimum_required(VERSION 2.8.12)

project(test)

option(WARNINGS "Show build warnings" ON)

# basic compiler options
# ----------------------

# set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# switch on warnings
if(WARNINGS)
  if(MSVC)
    if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
      string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    endif()
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
  endif()
endif()

# load packages
# -------------

# load pkg-config
find_package(PkgConfig)

# load header-only modules using pkg-config
# - cppmat
pkg_check_modules(CPPMAT REQUIRED cppmat)
include_directories(${CPPMAT_INCLUDE_DIRS})
# - eigen3
pkg_check_modules(EIGEN3 REQUIRED eigen3)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIRS})
# - GooseFEM
pkg_check_modules(GOOSEFEM REQUIRED GooseFEM)
include_directories(${GOOSEFEM_INCLUDE_DIRS})

# compile
# -------

add_executable(${PROJECT_NAME}
  main.cpp
  Vector.cpp
  MatrixDiagonal.cpp
  ElementQuad4.cpp
)
