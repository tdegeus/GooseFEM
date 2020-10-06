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

file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/include/GooseFEM/config.h" GooseFEM_version_defines
     REGEX "#define GOOSEFEM_VERSION_(MAJOR|MINOR|PATCH)")

foreach(ver ${GooseFEM_version_defines})
    if(ver MATCHES "#define GOOSEFEM_VERSION_(MAJOR|MINOR|PATCH) +([^ ]+)$")
        set(GOOSEFEM_VERSION_${CMAKE_MATCH_1}
            "${CMAKE_MATCH_2}"
            CACHE INTERNAL "")
    endif()
endforeach()

set(GOOSEFEM_VERSION ${GOOSEFEM_VERSION_MAJOR}.${GOOSEFEM_VERSION_MINOR}.${GOOSEFEM_VERSION_PATCH})

message(STATUS "Building GooseFEM v${GOOSEFEM_VERSION}")

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

if(BUILD_TESTS)
    add_subdirectory(test/basic)
    if (NOT WIN32)
        add_subdirectory(test/gmat)
    endif()
endif()

if(BUILD_EXAMPLES)
    add_subdirectory(docs/examples)
endif()
