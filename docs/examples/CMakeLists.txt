
cmake_minimum_required(VERSION 3.16)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    project(GooseFEM-examples)
    find_package(GooseFEM REQUIRED CONFIG)
endif()

option(WARNINGS "Show build warnings" ON)
option(SIMD "Enable xsimd" ON)
option(DEBUG "Enable all assertions" OFF)

set(HIGHFIVE_USE_BOOST 0)
set(HIGHFIVE_USE_XTENSOR 1)
find_package(HighFive REQUIRED)
find_package(XDMFWrite_HighFive REQUIRED)

add_library(libraries INTERFACE IMPORTED)

target_link_libraries(libraries INTERFACE GooseFEM HighFive XDMFWrite_HighFive)

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

# create executable

set(exec "statics_FixedDisplacements_LinearElastic_manual_partition")
set(source "statics/FixedDisplacements_LinearElastic/manual_partition.cpp")
add_executable(${exec} ${source})
target_link_libraries(${exec} PRIVATE libraries)

# create executable

set(exec "statics_MixedPeriodic_LinearElastic_example")
set(source "statics/MixedPeriodic_LinearElastic/example.cpp")
add_executable(${exec} ${source})
target_link_libraries(${exec} PRIVATE libraries)

# create executable

set(exec "statics_Periodic_ElastoPlastic_main")
set(source "statics/Periodic_ElastoPlastic/main.cpp")
add_executable(${exec} ${source})
target_link_libraries(${exec} PRIVATE libraries)

# create executable

set(exec "statics_Periodic_ElastoPlasticFiniteStrainSimo_main")
set(source "statics/Periodic_ElastoPlasticFiniteStrainSimo/main.cpp")
add_executable(${exec} ${source})
target_link_libraries(${exec} PRIVATE libraries)

# create executable

set(exec "statics_Periodic_LinearElastic_main")
set(source "statics/Periodic_LinearElastic/main.cpp")
add_executable(${exec} ${source})
target_link_libraries(${exec} PRIVATE libraries)

# create executable

set(exec "statics_Periodic_NonLinearElastic_main")
set(source "statics/Periodic_NonLinearElastic/main.cpp")
add_executable(${exec} ${source})
target_link_libraries(${exec} PRIVATE libraries)


