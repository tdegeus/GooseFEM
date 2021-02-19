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
    if(ver MATCHES "#define GOOSEFEM_VERSION_(MAJOR|MINOR|PATCH) ([0-9]*)(.*)")
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

# Adding git-hash to separate header
# ==================================

# Get the current working branch
execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE)

# Get the latest commit hash
execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE)

# Write a "git.h" overwrite with the current branch/hash
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/git.h"
    "#define GOOSEFEM_GIT_HASH \"${GIT_COMMIT_HASH}\"\n"
    "#define GOOSEFEM_GIT_BRANCH \"${GIT_BRANCH}\"\n")

# Overwrite the installed "include/GooseFEM/git.h" dummy
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/git.h" DESTINATION "include/GooseFEM/")

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
