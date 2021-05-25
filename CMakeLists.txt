cmake_minimum_required(VERSION 3.18..3.20)

# Basic settings
# ==============

project(GooseFEM)

option(BUILD_ALL "${PROJECT_NAME}: Build tests, Python API & docs" OFF)
option(BUILD_TESTS "${PROJECT_NAME}: Build tests" OFF)
option(BUILD_PYTHON "${PROJECT_NAME}: Build Python API" OFF)
option(BUILD_DOCS "${PROJECT_NAME}: Build docs" OFF)

execute_process(
    COMMAND python -c "from setuptools_scm import get_version; print(get_version())"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE PROJECT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE)

message(STATUS "Building ${PROJECT_NAME} ${PROJECT_VERSION}")

string(TOUPPER "${PROJECT_NAME}" PROJECT_NAME_UPPER)

# Set target
# ==========

find_package(xtensor REQUIRED)

add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(${PROJECT_NAME} INTERFACE
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

target_link_libraries(${PROJECT_NAME} INTERFACE xtensor)

target_compile_definitions(${PROJECT_NAME} INTERFACE
    ${PROJECT_NAME_UPPER}_VERSION="${PROJECT_VERSION}")

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

install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}-targets)

install(
    EXPORT ${PROJECT_NAME}-targets
    FILE "${PROJECT_NAME}Targets.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

set(${PROJECT_NAME}_TMP ${CMAKE_SIZEOF_VOID_P})
unset(CMAKE_SIZEOF_VOID_P)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY AnyNewerVersion)

set(CMAKE_SIZEOF_VOID_P ${${PROJECT_NAME}_TMP})

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}Config.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}.pc.in"
               "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc" @ONLY)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig/")

# Add builds
# ==========

include("${PROJECT_NAME}Config.cmake")
include(CTest)

# Tests

if(BUILD_TESTS OR BUILD_ALL)
    enable_testing()
    add_subdirectory(test/basic)

    enable_testing()
    add_subdirectory(test/gmat)
endif()

# Examples

if(BUILD_EXAMPLES)
    enable_testing()
    add_subdirectory(docs/examples)
endif()

# Python API

if(BUILD_PYTHON OR BUILD_ALL)
    find_package(Python REQUIRED COMPONENTS Interpreter Development NumPy)
    find_package(pybind11 REQUIRED CONFIG)
    find_package(xtensor-python REQUIRED)

    pybind11_add_module(mypymod python/main.cpp)
    set_target_properties(mypymod PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
    target_link_libraries(mypymod PUBLIC ${PROJECT_NAME}
        pybind11::module xtensor-python Python::NumPy)

    target_compile_definitions(mypymod PRIVATE VERSION_INFO=${PROJECT_VERSION})
endif()

# Documentation

if(BUILD_DOCS OR BUILD_ALL)
    find_package(Doxygen REQUIRED)

    set(DOXYGEN_EXCLUDE_SYMBOLS detail)
    set(DOXYGEN_CASE_SENSE_NAMES YES)
    set(DOXYGEN_USE_MATHJAX YES)
    set(DOXYGEN_GENERATE_TREEVIEW YES)
    set(DOXYGEN_JAVADOC_AUTOBRIEF YES)
    set(DOXYGEN_MACRO_EXPANSION YES)
    set(DOXYGEN_SOURCE_BROWSER YES)
    set(DOXYGEN_GENERATE_XML YES)
    set(DOXYGEN_QUIET YES)
    set(DOXYGEN_WARN_IF_UNDOCUMENTED YES)
    set(DOXYGEN_WARN_AS_ERROR YES)
    set(DOXYGEN_ALIASES "license=@par License:")
    set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "README.md")

    set(DOXYGEN_HTML_EXTRA_STYLESHEET
        "${CMAKE_CURRENT_SOURCE_DIR}/docs/doxystyle/doxygen-awesome.css"
        "${CMAKE_CURRENT_SOURCE_DIR}/docs/doxystyle/custom.css")

    set(DOXYGEN_STRIP_FROM_INC_PATH "${CMAKE_CURRENT_SOURCE_DIR}/include")
    set(DOXYGEN_STRIP_FROM_PATH "${CMAKE_CURRENT_SOURCE_DIR}/include")
    doxygen_add_docs(docs "${CMAKE_CURRENT_SOURCE_DIR}/include" "README.md")
endif()
