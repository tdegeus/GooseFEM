
cmake_minimum_required(VERSION 3.0)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    project(GooseFEM-examples)
    find_package(GooseFEM REQUIRED CONFIG)
endif()

option(WARNINGS "Show build warnings" ON)
option(SIMD "Enable xsimd" ON)
option(DEBUG "Enable all assertions" ON)

find_package(HDF5 REQUIRED)

add_library(libraries INTERFACE IMPORTED)

target_link_libraries(libraries INTERFACE GooseFEM ${HDF5_LIBRARIES})

if (SIMD)
    target_link_libraries(libraries INTERFACE xtensor::optimize xtensor::use_xsimd)
endif()

if (WARNINGS)
    target_link_libraries(libraries INTERFACE GooseFEM::compiler_warnings)
endif()

if (DEBUG)
    target_link_libraries(libraries INTERFACE GooseFEM::debug)
endif()


# create executable

set(exec "statics_FixedDisplacements_LinearElastic_example")
set(source "statics/FixedDisplacements_LinearElastic/example.cpp")

add_executable(${exec} ${source})

target_link_libraries(${exec} PRIVATE libraries)

configure_file(
     "statics/FixedDisplacements_LinearElastic/example.py"
     "statics_FixedDisplacements_LinearElastic_example.py" COPYONLY)

configure_file(
     "statics/FixedDisplacements_LinearElastic/plot.py"
     "statics_FixedDisplacements_LinearElastic_plot.py" COPYONLY)

# create executable

set(exec "statics_FixedDisplacements_LinearElastic_manual_partition")
set(source "statics/FixedDisplacements_LinearElastic/manual_partition.cpp")

add_executable(${exec} ${source})

target_link_libraries(${exec} PRIVATE libraries)

configure_file(
     "statics/FixedDisplacements_LinearElastic/manual_partition.py"
     "statics_FixedDisplacements_LinearElastic_manual_partition.py" COPYONLY)
