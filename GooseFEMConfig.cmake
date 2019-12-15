# GooseFEM cmake module
#
# This module sets the target:
#
#   GooseFEM
#
# In addition, it sets the following variables:
#
#   GooseFEM_FOUND - true if GooseFEM found
#   GooseFEM_VERSION - GooseFEM's version
#   GooseFEM_INCLUDE_DIRS - the directory containing GooseFEM headers
#
# The following support targets are defined to simplify things:
#
#   GooseFEM::compiler_warnings - enable compiler warnings
#   GooseFEM::assert - enable GooseFEM assertions
#   GooseFEM::debug - enable all assertions (slow)

include(CMakeFindDependencyMacro)

# Define target "GooseFEM"
# ========================

if(NOT TARGET GooseFEM)
    include("${CMAKE_CURRENT_LIST_DIR}/GooseFEMTargets.cmake")
    get_target_property(GooseFEM_INCLUDE_DIRS GooseFEM INTERFACE_INCLUDE_DIRECTORIES)
endif()

# Find dependencies
# =================

find_dependency(xtensor)

find_package(Eigen3 QUIET)

if(NOT Eigen3_FOUND)
    find_package(PkgConfig)
    pkg_check_modules(EIGEN3 QUIET eigen3)
endif()

if(Eigen3_FOUND)
    target_include_directories(GooseFEM INTERFACE ${EIGEN3_INCLUDE_DIRS})
endif()

# Define support target "GooseFEM::compiler_warnings"
# ===================================================

if(NOT TARGET GooseFEM::compiler_warnings)
    add_library(GooseFEM::compiler_warnings INTERFACE IMPORTED)
    if(MSVC)
        target_compile_options(GooseFEM::compiler_warnings INTERFACE
            /W4)
    else()
        target_compile_options(GooseFEM::compiler_warnings INTERFACE
            -Wall
            -Wextra
            -pedantic
            -Wno-unknown-pragmas)
    endif()
endif()

# Define support target "GooseFEM::assert"
# ========================================

if(NOT TARGET GooseFEM::assert)
    add_library(GooseFEM::assert INTERFACE IMPORTED)
    target_compile_definitions(GooseFEM::assert INTERFACE GOOSEFEM_ENABLE_ASSERT)
endif()

# Define support target "GooseEYE::debug"
# =======================================

if(NOT TARGET GooseFEM::debug)
    add_library(GooseFEM::debug INTERFACE IMPORTED)
    target_compile_definitions(GooseFEM::debug INTERFACE GOOSEFEM_DEBUG)
endif()
