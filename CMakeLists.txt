# required to specify the c++ standard
cmake_minimum_required(VERSION 3.0)

# required for install
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

# project settings
# ----------------

# name
project(GooseFEM)

# file that contains the version information
set(parse_version src/${PROJECT_NAME}/${PROJECT_NAME}.h)

# header files
set(headers
  src/GooseFEM/GooseFEM.h
  src/GooseFEM/Mesh.hpp
  src/GooseFEM/Mesh.h
  src/GooseFEM/MeshTri3.hpp
  src/GooseFEM/MeshTri3.h
  src/GooseFEM/MeshQuad4.hpp
  src/GooseFEM/MeshQuad4.h
  src/GooseFEM/MeshHex8.hpp
  src/GooseFEM/MeshHex8.h
  src/GooseFEM/Element.hpp
  src/GooseFEM/Element.h
  src/GooseFEM/ElementQuad4.hpp
  src/GooseFEM/ElementQuad4.h
  src/GooseFEM/ElementHex8.hpp
  src/GooseFEM/ElementHex8.h
  src/GooseFEM/Vector.hpp
  src/GooseFEM/Vector.h
  src/GooseFEM/MatrixDiagonal.hpp
  src/GooseFEM/MatrixDiagonal.h
  src/GooseFEM/Iterate.hpp
  src/GooseFEM/Iterate.h
  src/GooseFEM/Dynamics.hpp
  src/GooseFEM/Dynamics.h
)

# automatically parse the version number
file(READ "${parse_version}" version)
string(REGEX MATCH "define[ \t]+GOOSEFEM_WORLD_VERSION[ \t]+([0-9]+)" _ "${version}")
set(world_version "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+GOOSEFEM_MAJOR_VERSION[ \t]+([0-9]+)" _ "${version}")
set(major_version "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+GOOSEFEM_MINOR_VERSION[ \t]+([0-9]+)" _ "${version}")
set(minor_version "${CMAKE_MATCH_1}")
set(GOOSEFEM_VERSION_NUMBER ${world_version}.${major_version}.${minor_version})

# paths
# -----

set(GOOSEFEM_ROOT_DIR        "${CMAKE_INSTALL_PREFIX}")
set(GOOSEFEM_INCLUDE_DIR     "${CMAKE_INSTALL_PREFIX}/${INCLUDE_INSTALL_DIR}")
set(INCLUDE_INSTALL_DIR      "${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}")
set(CMAKEPACKAGE_INSTALL_DIR "${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}")
set(PKGCONFIG_INSTALL_DIR    "${CMAKE_INSTALL_DATADIR}/pkgconfig")
set(fcmake                   "${PROJECT_NAME}Config.cmake")
set(fpkg                     "${PROJECT_NAME}.pc")

# options
# -------

# configure pkg-config (default: on)
option(PKGCONFIG "Build pkg-config ${fpkg} file" ON)

# disable pkg-config for native Windows builds
if(WIN32 OR CMAKE_HOST_SYSTEM_NAME MATCHES Windows)
  option(PKGCONFIG "Build pkg-config ${fpkg} file" OFF)
endif()

# C++ standard
# ------------

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_CXX_STANDARD OR CMAKE_CXX_STANDARD LESS 14)
  set(CMAKE_CXX_STANDARD 14)
endif()

# configure CMake
# ---------------

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  PATH_VARS GOOSEFEM_INCLUDE_DIR GOOSEFEM_ROOT_DIR
  INSTALL_DESTINATION ${CMAKEPACKAGE_INSTALL_DIR}
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

# install
# -------

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
message(STATUS "|")
message(STATUS "| To specify a custom directory call")
message(STATUS "|   cmake /path/to/${PROJECT_NAME} -DCMAKE_INSTALL_PREFIX=yourprefix")
message(STATUS "|")
message(STATUS "| For custom paths, add the following line to your '~/.bashrc'")
message(STATUS "|   export PKG_CONFIG_PATH=${CMAKE_INSTALL_PREFIX}/share/pkgconfig:$PKG_CONFIG_PATH")
message(STATUS "|")
message(STATUS "+---------------------------------------------------------------------------------")
message(STATUS "")
