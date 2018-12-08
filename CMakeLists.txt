#==================================================================================================#
#                                                                                                  #
# (c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM       #
#                                                                                                  #
#==================================================================================================#

# initialization
# --------------

# required to specify the c++ standard
cmake_minimum_required(VERSION 3.0)

# required for install
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

# options
option(PKGCONFIG "Build pkg-config" ON)
option(BUILD_TESTS "${PROJECT_NAME} test suite" OFF)

# project settings
# ----------------

# name
project(GooseFEM)

# file that contains the version information
set(parse_version include/GooseFEM/config.h)

# header files
set(headers
  include/GooseFEM/config.h
  include/GooseFEM/GooseFEM.h
  include/GooseFEM/Mesh.hpp
  include/GooseFEM/Mesh.h
  include/GooseFEM/MeshTri3.hpp
  include/GooseFEM/MeshTri3.h
  include/GooseFEM/MeshQuad4.hpp
  include/GooseFEM/MeshQuad4.h
  include/GooseFEM/MeshHex8.hpp
  include/GooseFEM/MeshHex8.h
  include/GooseFEM/Element.hpp
  include/GooseFEM/Element.h
  include/GooseFEM/ElementQuad4.hpp
  include/GooseFEM/ElementQuad4.h
  include/GooseFEM/ElementQuad4Planar.hpp
  include/GooseFEM/ElementQuad4Planar.h
  include/GooseFEM/ElementQuad4Axisymmetric.hpp
  include/GooseFEM/ElementQuad4Axisymmetric.h
  include/GooseFEM/ElementHex8.hpp
  include/GooseFEM/ElementHex8.h
  include/GooseFEM/Vector.hpp
  include/GooseFEM/Vector.h
  include/GooseFEM/VectorPartitioned.hpp
  include/GooseFEM/VectorPartitioned.h
  include/GooseFEM/Matrix.hpp
  include/GooseFEM/Matrix.h
  include/GooseFEM/MatrixPartitioned.hpp
  include/GooseFEM/MatrixPartitioned.h
  include/GooseFEM/MatrixDiagonal.hpp
  include/GooseFEM/MatrixDiagonal.h
  include/GooseFEM/MatrixDiagonalPartitioned.hpp
  include/GooseFEM/MatrixDiagonalPartitioned.h
  include/GooseFEM/Iterate.hpp
  include/GooseFEM/Iterate.h
)

# automatically parse the version number
file(READ "${parse_version}" version)
string(REGEX MATCH "define[ \t]+GOOSEFEM_WORLD_VERSION[ \t]+([0-9]+)" _ "${version}")
set(world_version "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+GOOSEFEM_MAJOR_VERSION[ \t]+([0-9]+)" _ "${version}")
set(major_version "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+GOOSEFEM_MINOR_VERSION[ \t]+([0-9]+)" _ "${version}")
set(minor_version "${CMAKE_MATCH_1}")
set(GOOSEFEM_VERSION ${world_version}.${major_version}.${minor_version})

# print information to screen
message(STATUS "Building ${PROJECT_NAME} v${GOOSEFEM_VERSION}")

# paths
# -----

set(GOOSEFEM_ROOT_DIR        "${CMAKE_INSTALL_PREFIX}")
set(GOOSEFEM_INCLUDE_DIR     "${CMAKE_INSTALL_PREFIX}/${INCLUDE_INSTALL_DIR}")
set(INCLUDE_INSTALL_DIR      "${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}")
set(CMAKEPACKAGE_INSTALL_DIR "${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}")
set(PKGCONFIG_INSTALL_DIR    "${CMAKE_INSTALL_DATADIR}/pkgconfig")
set(fcmake                   "GooseFEMConfig.cmake")
set(fpkg                     "GooseFEM.pc")

# configure CMake
# ---------------

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/GooseFEMConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/GooseFEMConfig.cmake
  PATH_VARS GOOSEFEM_INCLUDE_DIR GOOSEFEM_ROOT_DIR
  INSTALL_DESTINATION ${CMAKEPACKAGE_INSTALL_DIR}
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

# install/build
# -------------

# build test cases
if(BUILD_TESTS)
  add_subdirectory(test)
endif()

# disable pkg-config for native Windows builds
if(WIN32 OR CMAKE_HOST_SYSTEM_NAME MATCHES Windows)
  option(PKGCONFIG "Build pkg-config ${fpkg} file" OFF)
endif()

# pkg-config
if(PKGCONFIG)
  configure_file(${fpkg}.in ${fpkg} @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${fpkg} DESTINATION ${PKGCONFIG_INSTALL_DIR})
endif()

# CMake
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${fcmake} DESTINATION ${CMAKEPACKAGE_INSTALL_DIR})

# header files
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/${headers} DESTINATION ${INCLUDE_INSTALL_DIR})

# print information to screen
# ---------------------------

message(STATUS "")
message(STATUS "+---------------------------------------------------------------------------------")
message(STATUS "|")
message(STATUS "| Use 'make install' to install in")
message(STATUS "|   ${CMAKE_INSTALL_PREFIX}")
if(BUILD_TESTS)
  message(STATUS "|")
  message(STATUS "| Use 'make' and './test/test' to run the tests")
endif()
message(STATUS "|")
message(STATUS "| To specify a custom directory call")
message(STATUS "|   cmake /path/to/${PROJECT_NAME} -DCMAKE_INSTALL_PREFIX=yourprefix")
message(STATUS "|")
message(STATUS "| For custom paths, add the following line to your '~/.bashrc'")
message(STATUS "|   export PKG_CONFIG_PATH=${CMAKE_INSTALL_PREFIX}/share/pkgconfig:$PKG_CONFIG_PATH")
message(STATUS "|")
message(STATUS "+---------------------------------------------------------------------------------")
message(STATUS "")
