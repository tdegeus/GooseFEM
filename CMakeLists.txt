#
# (c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM
#

cmake_minimum_required(VERSION 3.0)

# Basic settings
# ==============

project(GooseFEM)

option(BUILD_TESTS "Build tests" OFF)
option(BUILD_EXAMPLES "Build examples" OFF)

# Version
# =======

IF(DEFINED ENV{PKG_VERSION})
    set(GOOSEFEM_VERSION $ENV{PKG_VERSION})
else()
    execute_process(
        COMMAND python -c "from setuptools_scm import get_version; print(get_version())"
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GOOSEFEM_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

message(STATUS "Building GooseFEM ${GOOSEFEM_VERSION}")

# Set target
# ==========

find_package(xtensor REQUIRED)

add_library(GooseFEM INTERFACE)

target_include_directories(GooseFEM INTERFACE
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

target_link_libraries(GooseFEM INTERFACE xtensor)

# Installation
# ============

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/" DESTINATION include)

configure_file("include/${PROJECT_NAME}/version.h"
               "${CMAKE_CURRENT_BINARY_DIR}/version.h"
                @ONLY)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/version.h"
        DESTINATION "include/${PROJECT_NAME}/")

install(TARGETS GooseFEM EXPORT GooseFEM-targets)

install(
    EXPORT GooseFEM-targets
    FILE GooseFEMTargets.cmake
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/GooseFEM")

set(_GOOSEFEM ${CMAKE_SIZEOF_VOID_P})
unset(CMAKE_SIZEOF_VOID_P)

write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/GooseFEMConfigVersion.cmake"
    VERSION ${GOOSEFEM_VERSION}
    COMPATIBILITY AnyNewerVersion)

set(CMAKE_SIZEOF_VOID_P ${_GOOSEFEM})

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/GooseFEMConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/GooseFEMConfigVersion.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/GooseFEM")

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GooseFEM.pc.in"
               "${CMAKE_CURRENT_BINARY_DIR}/GooseFEM.pc" @ONLY)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/GooseFEM.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig/")

# Add builds
# ==========

include("GooseFEMConfig.cmake")

include(CTest)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(test/basic)

    enable_testing()
    add_subdirectory(test/gmat)
endif()

if(BUILD_EXAMPLES)
    enable_testing()
    add_subdirectory(docs/examples)
endif()
