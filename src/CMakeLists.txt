# required to specify the c++ standard
cmake_minimum_required(VERSION 3.0)

# define project information
project(GooseFEM)

# automatically parse the version number
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/GooseFEM/GooseFEM.h" _goosefem_version_header)
string(REGEX MATCH "define[ \t]+GOOSEFEM_WORLD_VERSION[ \t]+([0-9]+)" _goosefem_world_version_match "${_goosefem_version_header}")
set(GOOSEFEM_WORLD_VERSION "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+GOOSEFEM_MAJOR_VERSION[ \t]+([0-9]+)" _goosefem_major_version_match "${_goosefem_version_header}")
set(GOOSEFEM_MAJOR_VERSION "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+GOOSEFEM_MINOR_VERSION[ \t]+([0-9]+)" _goosefem_minor_version_match "${_goosefem_version_header}")
set(GOOSEFEM_MINOR_VERSION "${CMAKE_MATCH_1}")
set(GOOSEFEM_VERSION_NUMBER ${GOOSEFEM_WORLD_VERSION}.${GOOSEFEM_MAJOR_VERSION}.${GOOSEFEM_MINOR_VERSION})

# required for install
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

# set c++ standard -> c++14
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_CXX_STANDARD OR CMAKE_CXX_STANDARD LESS 14)
  set(CMAKE_CXX_STANDARD 14)
endif()

# set paths where 'GooseFEM' will be installed
# - root
set(GOOSEFEM_ROOT_DIR ${CMAKE_INSTALL_PREFIX})
# - headers
set(GOOSEFEM_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/${INCLUDE_INSTALL_DIR}")
# - headers
set(INCLUDE_INSTALL_DIR
  "${CMAKE_INSTALL_INCLUDEDIR}/GooseFEM"
  CACHE PATH "Path, relative to CMAKE_PREFIX_PATH, where the GooseFEM header files are installed"
)
# - CMake configuration
set(CMAKEPACKAGE_INSTALL_DIR
  "${CMAKE_INSTALL_DATADIR}/GooseFEM"
  CACHE PATH "Path where GooseFEMConfig.cmake is installed"
)
# - pkg-config
set(PKGCONFIG_INSTALL_DIR
  "${CMAKE_INSTALL_DATADIR}/pkgconfig"
  CACHE PATH "Path where GooseFEM.pc is installed"
)

# list all header files
set(GOOSEFEM_HEADERS
  GooseFEM/GooseFEM.h
  GooseFEM/Mesh.h
  GooseFEM/Mesh.cpp
  GooseFEM/MeshTri3.h
  GooseFEM/MeshTri3.cpp
  GooseFEM/MeshQuad4.h
  GooseFEM/MeshQuad4.cpp
  GooseFEM/MeshHex8.h
  GooseFEM/MeshHex8.cpp
  GooseFEM/DynamicsDiagonal.h
  GooseFEM/DynamicsDiagonal.cpp
  GooseFEM/OverdampedDynamicsDiagonal.h
  GooseFEM/OverdampedDynamicsDiagonal.cpp
)

# Disable pkg-config for native Windows builds
if(NOT WIN32 OR NOT CMAKE_HOST_SYSTEM_NAME MATCHES Windows)
  option(GOOSEFEM_BUILD_PKGCONFIG "Build pkg-config .pc file for GooseFEM" ON)
endif()

# build pkg-config file -> fill in relevant values
if(GOOSEFEM_BUILD_PKGCONFIG)
  configure_file(GooseFEM.pc.in GooseFEM.pc @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/GooseFEM.pc DESTINATION ${PKGCONFIG_INSTALL_DIR})
endif()

# configure CMake file
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/GooseFEMConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/GooseFEMConfig.cmake
  PATH_VARS GOOSEFEM_INCLUDE_DIR GOOSEFEM_ROOT_DIR
  INSTALL_DESTINATION ${CMAKEPACKAGE_INSTALL_DIR}
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

# install
# - CMake-file
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/GooseFEMConfig.cmake
  DESTINATION
  ${CMAKEPACKAGE_INSTALL_DIR}
)
# - headers
install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/${GOOSEFEM_HEADERS}
  DESTINATION
  ${INCLUDE_INSTALL_DIR}/GooseFEM
)

# print information to screen
message(STATUS "")
message(STATUS "+-------------------------------------------------------------------------------")
message(STATUS "|")
message(STATUS "| Use 'make install' to install in '${CMAKE_INSTALL_PREFIX}'")
message(STATUS "| To specify a custom directory call")
message(STATUS "|   cmake .. -DCMAKE_INSTALL_PREFIX=yourprefix")
message(STATUS "| For custom paths, add the following line to your '~/.bashrc'")
message(STATUS "|   export PKG_CONFIG_PATH=${CMAKE_INSTALL_PREFIX}/share/pkgconfig:$PKG_CONFIG_PATH")
message(STATUS "|")
message(STATUS "+-------------------------------------------------------------------------------")
message(STATUS "")
