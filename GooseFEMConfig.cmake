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

if(NOT TARGET GooseFEM)
    include("${CMAKE_CURRENT_LIST_DIR}/GooseFEMTargets.cmake")
    get_target_property(GooseFEM_INCLUDE_DIRS GooseFEM INTERFACE_INCLUDE_DIRECTORIES)
endif()

# Find dependencies

find_dependency(xtensor)

find_package(Eigen3)

if (TARGET Eigen3::Eigen)
    target_link_libraries(GooseFEM INTERFACE Eigen3::Eigen)
endif()

# Define support target "GooseFEM::compiler_warnings"

if(NOT TARGET GooseFEM::compiler_warnings)
    add_library(GooseFEM::compiler_warnings INTERFACE IMPORTED)
    if(MSVC)
        set_property(
            TARGET GooseFEM::compiler_warnings
            PROPERTY INTERFACE_COMPILE_OPTIONS
            /W4)
    else()
        set_property(
            TARGET GooseFEM::compiler_warnings
            PROPERTY INTERFACE_COMPILE_OPTIONS
            -Wall -Wextra -pedantic -Wno-unknown-pragmas)
    endif()
endif()

# Define support target "GooseFEM::warnings"

if(NOT TARGET GooseFEM::warnings)
    add_library(GooseFEM::warnings INTERFACE IMPORTED)
    set_property(
        TARGET GooseFEM::warnings
        PROPERTY INTERFACE_COMPILE_DEFINITIONS
        GOOSEFEM_ENABLE_WARNING_PYTHON)
endif()

# Define support target "GooseFEM::assert"

if(NOT TARGET GooseFEM::assert)
    add_library(GooseFEM::assert INTERFACE IMPORTED)
    set_property(
        TARGET GooseFEM::assert
        PROPERTY INTERFACE_COMPILE_DEFINITIONS
        GOOSEFEM_ENABLE_ASSERT)
endif()

# Define support target "GooseFEM::debug"

if(NOT TARGET GooseFEM::debug)
    add_library(GooseFEM::debug INTERFACE IMPORTED)
    set_property(
        TARGET GooseFEM::debug
        PROPERTY INTERFACE_COMPILE_DEFINITIONS
        XTENSOR_ENABLE_ASSERT GOOSEFEM_ENABLE_ASSERT)
endif()
