# required to specify the c++ standard
cmake_minimum_required(VERSION 3.0)

# required for install
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

# project settings
# ----------------

# define project information
project(GooseFEM)

# list of header files
set(headers
  GooseFEM/GooseFEM.h
  GooseFEM/Mesh.cpp
  GooseFEM/Mesh.h
  GooseFEM/MeshTri3.cpp
  GooseFEM/MeshTri3.h
  GooseFEM/MeshQuad4.cpp
  GooseFEM/MeshQuad4.h
  GooseFEM/MeshHex8.cpp
  GooseFEM/MeshHex8.h
  GooseFEM/Element.cpp
  GooseFEM/Element.h
  GooseFEM/ElementQuad4.cpp
  GooseFEM/ElementQuad4.h
  GooseFEM/Vector.cpp
  GooseFEM/Vector.h
  GooseFEM/MatrixDiagonal.cpp
  GooseFEM/MatrixDiagonal.h
  GooseFEM/Iterate.cpp
  GooseFEM/Iterate.h
  GooseFEM/Dynamics.cpp
  GooseFEM/Dynamics.h
)

# options
# -------

# set options
option(PKGCONFIG "Build pkg-config ${PROJECT_NAME}.pc file" ON)

# disable pkg-config for native Windows builds
if(WIN32 OR CMAKE_HOST_SYSTEM_NAME MATCHES Windows)
  option(PKGCONFIG "Build pkg-config ${PROJECT_NAME}.pc file" OFF)
endif()

# paths
# -----

set(GOOSEFEM_ROOT_DIR        "${CMAKE_INSTALL_PREFIX}")
set(GOOSEFEM_INCLUDE_DIR     "${CMAKE_INSTALL_PREFIX}/${INCLUDE_INSTALL_DIR}")
set(INCLUDE_INSTALL_DIR      "${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}")
set(CMAKEPACKAGE_INSTALL_DIR "${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}")
set(PKGCONFIG_INSTALL_DIR    "${CMAKE_INSTALL_DATADIR}/pkgconfig")

# version
# -------

# automatically parse the version number
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}/${PROJECT_NAME}.h" version)
string(REGEX MATCH "define[ \t]+GOOSEFEM_WORLD_VERSION[ \t]+([0-9]+)" _ "${version}")
set(world_version "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+GOOSEFEM_MAJOR_VERSION[ \t]+([0-9]+)" _ "${version}")
set(major_version "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+GOOSEFEM_MINOR_VERSION[ \t]+([0-9]+)" _ "${version}")
set(minor_version "${CMAKE_MATCH_1}")
set(GOOSEFEM_VERSION_NUMBER ${world_version}.${major_version}.${minor_version})

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
  configure_file(${PROJECT_NAME}.pc.in ${PROJECT_NAME}.pc @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc DESTINATION ${PKGCONFIG_INSTALL_DIR})
endif()

# CMake
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  DESTINATION
  ${CMAKEPACKAGE_INSTALL_DIR}
)

# headers
install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/${headers}
  DESTINATION
  ${INCLUDE_INSTALL_DIR}/${PROJECT_NAME}
)

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
