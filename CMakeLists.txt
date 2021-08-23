cmake_minimum_required(VERSION 3.12)
project(photo-decomp)
enable_language(Fortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

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
set(NETCDF_LIBS ${NETCDF_LIB} ${NETCDF_FORTRAN_LIB})
include_directories(${NETCDF_INCLUDE_DIR})

################################################################################
# nc4fortran library

find_path(NC4F_INCLUDE_DIR nc4fortran.mod
  DOC "nc4fortran include directory (must contain nc4fortran.mod and related *.smod files)"
  PATHS
    /usr/lib/gfortran/modules
    /usr/lib64/gfortran/modules
    /usr/local/include
    /opt/local/include)
find_library(NC4F_LIB nc4fortran
  DOC "nc4fortran library"
  PATHS
    /usr/local/lib
    /opt/local/lib)
include_directories(${NC4F_INCLUDE_DIR})

################################################################################
# musica-core library

add_subdirectory(libs/musica-core/src)
include_directories(${CMAKE_BINARY_DIR}/libs/musica-core/src)

################################################################################
# Photo-decomp source

add_custom_target(copy-DATAE1 ALL COMMAND ${CMAKE_COMMAND}
  -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data/DATAE1 ${CMAKE_BINARY_DIR}/DATAE1)
add_custom_target(copy-DATAJ1 ALL COMMAND ${CMAKE_COMMAND}
  -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data/DATAJ1 ${CMAKE_BINARY_DIR}/DATAJ1)
add_custom_target(copy-DATAS1 ALL COMMAND ${CMAKE_COMMAND}
  -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data/DATAS1 ${CMAKE_BINARY_DIR}/DATAS1)
add_custom_target(copy-INPUTS ALL COMMAND ${CMAKE_COMMAND}
  -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/test/data/INPUTS ${CMAKE_BINARY_DIR}/INPUTS)
add_custom_target(copy-OUTPUTS ALL COMMAND ${CMAKE_COMMAND}
  -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/test/data/OUTPUTS ${CMAKE_BINARY_DIR}/OUTPUTS)
add_subdirectory(src)
include_directories(${CMAKE_BINARY_DIR}/src)

################################################################################
# Photo-decomp tests

enable_testing()
add_subdirectory(test)

################################################################################
