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
project(xGooseFEM)

# file that contains the version information
set(parse_version include/xGooseFEM/GooseFEM.h)

# header files
set(headers
  include/xGooseFEM/GooseFEM.h
  include/xGooseFEM/Mesh.hpp
  include/xGooseFEM/Mesh.h
  include/xGooseFEM/MeshTri3.hpp
  include/xGooseFEM/MeshTri3.h
  include/xGooseFEM/MeshQuad4.hpp
  include/xGooseFEM/MeshQuad4.h
  include/xGooseFEM/MeshHex8.hpp
  include/xGooseFEM/MeshHex8.h
  include/xGooseFEM/Element.hpp
  include/xGooseFEM/Element.h
  include/xGooseFEM/ElementQuad4.hpp
  include/xGooseFEM/ElementQuad4.h
  include/xGooseFEM/ElementHex8.hpp
  include/xGooseFEM/ElementHex8.h
  include/xGooseFEM/Vector.hpp
  include/xGooseFEM/Vector.h
  include/xGooseFEM/MatrixDiagonal.hpp
  include/xGooseFEM/MatrixDiagonal.h
  include/xGooseFEM/Iterate.hpp
  include/xGooseFEM/Iterate.h
  include/xGooseFEM/Dynamics.hpp
  include/xGooseFEM/Dynamics.h
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
  message(STATUS "| Use 'make' and './test' to run the tests")
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
