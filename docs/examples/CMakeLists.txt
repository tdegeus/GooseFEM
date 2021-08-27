cmake_minimum_required(VERSION 3.1..3.21)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    project(QPot)
    find_package(QPot REQUIRED CONFIG)
    option(USE_ASSERT "${PROJECT_NAME}: Build with assertions" ON)
    option(USE_DEBUG "${PROJECT_NAME}: Build in debug mode" OFF)
    option(USE_SIMD "${PROJECT_NAME}: Build with hardware optimization" OFF)
endif()

set(MYPROJECT "${PROJECT_NAME}-examples")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(HIGHFIVE_USE_BOOST 0)
set(HIGHFIVE_USE_XTENSOR 1)
find_package(HighFive REQUIRED)
find_package(XDMFWrite_HighFive REQUIRED)
find_package(GMatElastic REQUIRED)
find_package(GMatElastoPlastic REQUIRED)
find_package(GMatElastoPlasticFiniteStrainSimo REQUIRED)
find_package(GMatNonLinearElastic REQUIRED)

add_library(mytarget INTERFACE IMPORTED)

target_link_libraries(mytarget INTERFACE
    ${PROJECT_NAME}
    ${PROJECT_NAME}::compiler_warnings
    HighFive
    XDMFWrite_HighFive
    GMatElastic
    GMatElastoPlastic
    GMatElastoPlasticFiniteStrainSimo
    GMatNonLinearElastic)

if(USE_ASSERT)
    target_link_libraries(mytarget INTERFACE ${PROJECT_NAME}::assert)
    message(STATUS "Compiling ${MYPROJECT} with assertions")
endif()

if(USE_DEBUG)
    target_link_libraries(mytarget INTERFACE ${PROJECT_NAME}::debug)
    message(STATUS "Compiling ${MYPROJECT} in debug mode")
endif()

if(USE_SIMD)
    find_package(xtensor REQUIRED)
    find_package(xsimd REQUIRED)
    target_link_libraries(mytarget INTERFACE xtensor::use_xsimd xtensor::optimize)
    message(STATUS "Compiling ${MYPROJECT} with hardware optimization")
endif()

set(myexec "statics_FixedDisplacements_LinearElastic_example")
set(mysource "statics/FixedDisplacements_LinearElastic/example.cpp")
add_executable(${myexec} ${mysource})
target_link_libraries(${myexec} PRIVATE mytarget)
add_test(NAME ${myexec} COMMAND ${myexec})


set(myexec "statics_FixedDisplacements_LinearElastic_manual_partition")
set(mysource "statics/FixedDisplacements_LinearElastic/manual_partition.cpp")
add_executable(${myexec} ${mysource})
target_link_libraries(${myexec} PRIVATE mytarget)
add_test(NAME ${myexec} COMMAND ${myexec})


set(myexec "statics_MixedPeriodic_LinearElastic_example")
set(mysource "statics/MixedPeriodic_LinearElastic/example.cpp")
add_executable(${myexec} ${mysource})
target_link_libraries(${myexec} PRIVATE mytarget)
add_test(NAME ${myexec} COMMAND ${myexec})


set(myexec "statics_Periodic_ElastoPlastic_main")
set(mysource "statics/Periodic_ElastoPlastic/main.cpp")
add_executable(${myexec} ${mysource})
target_link_libraries(${myexec} PRIVATE mytarget)
add_test(NAME ${myexec} COMMAND ${myexec})


set(myexec "statics_Periodic_ElastoPlasticFiniteStrainSimo_main")
set(mysource "statics/Periodic_ElastoPlasticFiniteStrainSimo/main.cpp")
add_executable(${myexec} ${mysource})
target_link_libraries(${myexec} PRIVATE mytarget)
add_test(NAME ${myexec} COMMAND ${myexec})


set(myexec "statics_Periodic_LinearElastic_main")
set(mysource "statics/Periodic_LinearElastic/main.cpp")
add_executable(${myexec} ${mysource})
target_link_libraries(${myexec} PRIVATE mytarget)
add_test(NAME ${myexec} COMMAND ${myexec})


set(myexec "statics_Periodic_NonLinearElastic_main")
set(mysource "statics/Periodic_NonLinearElastic/main.cpp")
add_executable(${myexec} ${mysource})
target_link_libraries(${myexec} PRIVATE mytarget)
add_test(NAME ${myexec} COMMAND ${myexec})


set(myexec "dynamics_Elastic-Verlet_main")
set(mysource "dynamics/Elastic-Verlet/main.cpp")
add_executable(${myexec} ${mysource})
target_link_libraries(${myexec} PRIVATE mytarget)
# add_test(NAME ${myexec} COMMAND ${myexec}) # takes too long for testing
