cmake_minimum_required(VERSION 3.13)
set(CMAKE_USER_MAKE_RULES_OVERRIDE "set_defaults.cmake")

################################################################################
# project
project(tuv-x)
enable_language(Fortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake-modules)

################################################################################
# options

include(CMakeDependentOption)
option(ENABLE_MPI "Enable MPI parallel support" OFF)
cmake_dependent_option(ENABLE_OPENMP "Enable OpenMP support" OFF "ENABLE_MPI" OFF)
option(ENABLE_TESTS "Build tests" ON)
option(ENABLE_COVERAGE "Enable code coverage output" ON)
option(ENABLE_MEMCHECK "Enable memory checking in tests" ON)
option(ENABLE_NC_CONFIG "Use nc-config to determine NetCDF libraries" OFF)

################################################################################
# Set up include and lib directories

set(TUVX_MOD_DIR "${PROJECT_BINARY_DIR}/include")
set(TUVX_LIB_DIR "${PROJECT_BINARY_DIR}/lib")

################################################################################
# Compiler specific flags

if(${CMAKE_Fortran_COMPILER_ID} MATCHES "NAG")
  add_definitions(-DMUSICA_IS_NAG_COMPILER)
endif()

################################################################################
# MPI

if(ENABLE_MPI)
  add_definitions(-DMUSICA_USE_MPI)
endif()

################################################################################
# Memory check

if(ENABLE_MEMCHECK)
  find_file(MEMCHECK_SUPPRESS_FILE
    DOC "Suppression file for memory checking"
    NAMES openmpi-valgrind.supp
    PATHS
      /usr/share/openmpi
      /usr/lib64/openmpi/share
      /usr/lib64/openmpi/share/openmpi
      /usr/share)
  if(MEMCHECK_SUPPRESS_FILE)
    set(MEMCHECK_SUPPRESS "--suppressions=${PROJECT_SOURCE_DIR}/test/valgrind.supp --suppressions=${MEMCHECK_SUPPRESS_FILE}")
  else()
    set(MEMCHECK_SUPPRESS "--suppressions=${PROJECT_SOURCE_DIR}/test/valgrind.supp")
  endif()
endif()

################################################################################
# OpenMP

if(ENABLE_OPENMP)
  find_package(OpenMP)
  if(OpenMP_Fortran_FOUND)
    message(STATUS "Compiling with OpenMP support")
    add_definitions(-DMUSICA_USE_OPENMP)
  else()
    message(FATAL_ERROR "OpenMP package not found")
  endif()
endif()

################################################################################
# Test code coverage

if(ENABLE_COVERAGE)
  include(CodeCoverage)
  if(CMAKE_COMPILER_IS_GNUCXX)
    setup_target_for_coverage_lcov(
        NAME coverage
        EXECUTABLE "ctest"
        EXCLUDE "${PROJECT_SOURCE_DIR}/test/*"
        BASE_DIRECTORY "${PROJECT_SOURCE_DIR}/src")
  endif()
endif()

################################################################################
# json-fortran library

find_path(JSON_INCLUDE_DIR json_module.mod
  DOC "json-fortran include directory (must include json_*.mod files)"
  PATHS
    $ENV{JSON_FORTRAN_HOME}/lib
    /opt/local/lib
    /usr/local/lib
    /usr/local/lib64)
find_library(JSON_LIB jsonfortran
  DOC "json-fortran library"
  PATHS
    $ENV{JSON_FORTRAN_HOME}/lib
    /opt/local/lib
    /usr/local/lib
    /usr/local/lib64)
include_directories(${JSON_INCLUDE_DIR})

################################################################################
# NetCDF library

find_path(NETCDF_INCLUDE_DIR netcdf.mod NETCDF.mod
  DOC "NetCDF include directory (must contain netcdf.mod)"
  PATHS
    $ENV{NETCDF_HOME}/include
    /usr/lib/gfortran/modules
    /usr/lib64/gfortran/modules
    /opt/local/include)
find_library(NETCDF_C_LIB netcdf
  DOC "NetCDF C library"
  PATHS
    $ENV{NETCDF_HOME}/lib
    $ENV{NETCDF_HOME}/lib64
    opt/local/lib)
find_library(NETCDF_FORTRAN_LIB netcdff
  DOC "NetCDF Fortran library"
  PATHS
    $ENV{NETCDF_HOME}/lib
    $ENV{NETCDF_HOME}/lib64
    /opt/local/lib)
if(ENABLE_NC_CONFIG)
  execute_process(COMMAND nc-config --flibs OUTPUT_VARIABLE NC_CONFIG_FLIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(NETCDF_LIBS ${NC_CONFIG_FLIBS} )
else()
  set(NETCDF_LIBS ${NETCDF_LIB} ${NETCDF_FORTRAN_LIB})
endif()
include_directories(${NETCDF_INCLUDE_DIR})

################################################################################
# musica-core library

find_path(MUSICA_CORE_INCLUDE_DIR musica_math.mod
  DOC "musica-core include directory (must include musica_*.mod files)"
  PATHS
    $ENV{MUSICA_CORE_HOME}/include
    /opt/local/lib
    /usr/local/lib
    /usr/local/lib64)
find_library(MUSICA_CORE_LIB musica_core
  DOC "musica-core library"
  PATHS
    $ENV{MUSICA_CORE_HOME}/lib
    /opt/local/lib
    /usr/local/lib
    /usr/local/lib64)
include_directories(${MUSICA_CORE_INCLUDE_DIR})

################################################################################
# copy data

add_custom_target(copy-data ALL COMMAND ${CMAKE_COMMAND}
  -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_BINARY_DIR}/data)
add_custom_target(copy-tool ALL COMMAND ${CMAKE_COMMAND}
  -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/tool ${CMAKE_BINARY_DIR}/tool)
add_custom_target(copy-examples ALL COMMAND ${CMAKE_COMMAND}
  -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/examples ${CMAKE_BINARY_DIR}/examples)

################################################################################
# TUV-x source

add_subdirectory(src)
include_directories(${CMAKE_BINARY_DIR}/src)

################################################################################
# TUV-x tests

enable_testing()
if(ENABLE_TESTS)
  add_subdirectory(test)
endif()

################################################################################
