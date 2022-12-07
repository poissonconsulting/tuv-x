################################################################################
# Preamble

cmake_minimum_required(VERSION 3.21)
set(CMAKE_USER_MAKE_RULES_OVERRIDE "cmake/set_defaults.cmake")

project(
  tuv-x
  VERSION 0.2.0
  LANGUAGES Fortran
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/cmake/cmake-modules")

message ( STATUS "CMake build configuration for tuv-x(${CMAKE_BUILD_TYPE}) ${PROJECT_VERSION}" )

################################################################################
# Projet wide setup options

include(CMakeDependentOption)
option(ENABLE_MPI "Enable MPI parallel support" OFF)
cmake_dependent_option(ENABLE_OPENMP "Enable OpenMP support" OFF "ENABLE_MPI" OFF)
option(ENABLE_TESTS "Build tests" ON)
option(ENABLE_COVERAGE "Enable code coverage output" ON)
option(ENABLE_MEMCHECK "Enable memory checking in tests" ON)
option(ENABLE_NC_CONFIG "Use nc-config to determine NetCDF libraries" OFF)

# Set up include and lib directories
set(TUVX_MOD_DIR "${PROJECT_BINARY_DIR}/include")
set(TUVX_LIB_DIR "${PROJECT_BINARY_DIR}/lib")
include(GNUInstallDirs)
set(
  INSTALL_PREFIX 
  "${CMAKE_INSTALL_PREFIX}/ncar/tuvx-${PROJECT_VERSION}-${CMAKE_Fortran_COMPILER_ID}"
)
set(INSTALL_MOD_DIR "${INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}")

# Compiler specific flags
if(${CMAKE_Fortran_COMPILER_ID} MATCHES "NAG")
  add_definitions(-DMUSICA_IS_NAG_COMPILER)
endif()

# MPI
if(ENABLE_MPI)
  add_definitions(-DMUSICA_USE_MPI)
endif()

# Test code coverage
if(ENABLE_COVERAGE)
  include(CodeCoverage)
  if(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
    setup_target_for_coverage_lcov(
        NAME coverage
        EXECUTABLE "ctest"
        EXCLUDE "${PROJECT_SOURCE_DIR}/test/*"
        BASE_DIRECTORY "${PROJECT_SOURCE_DIR}/src")
  endif()
endif()

# copy data
add_custom_target(copy-data ALL COMMAND ${CMAKE_COMMAND}
  -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_BINARY_DIR}/data)
add_custom_target(copy-tool ALL COMMAND ${CMAKE_COMMAND}
  -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/tool ${CMAKE_BINARY_DIR}/tool)
add_custom_target(copy-examples ALL COMMAND ${CMAKE_COMMAND}
  -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/examples ${CMAKE_BINARY_DIR}/examples)

################################################################################
# Dependencies

include(cmake/dependencies.cmake)

################################################################################
# TUV-x targets

# add alias so that this project can be included with add_subdirectory or 
# be used with find_pacakge
add_library(tuvxlib)
add_library(tuvx::tuvx ALIAS tuvxlib)

set_target_properties(tuvxlib PROPERTIES 
  OUTPUT_NAME tuvx
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
  ARCHIVE_OUTPUT_DIRECTORY ${TUVX_LIB_DIR}
  Fortran_MODULE_DIRECTORY ${TUVX_MOD_DIR}
)

add_executable(tuv-x src/tuvx.F90)

add_subdirectory(src)

target_link_libraries(tuvxlib 
  PUBLIC
    musicacore::musicacore
    ${NETCDF_LIBS} 
    ${JSON_LIB}
)

target_include_directories(tuvxlib 
  PUBLIC  
    $<BUILD_INTERFACE:${TUVX_MOD_DIR}>
    $<INSTALL_INTERFACE:${INSTALL_MOD_DIR}>
)

target_link_libraries(tuv-x PUBLIC tuvxlib)
if(ENABLE_OPENMP)
  target_link_libraries(tuv-x PUBLIC OpenMP::OpenMP_Fortran)
endif()


################################################################################
# TUV-x tests

enable_testing()
if(ENABLE_TESTS)
  add_subdirectory(test)
endif()

################################################################################
# Packaging

# only include packaging if we are the top level project being built
if(PROJECT_IS_TOP_LEVEL)
  add_subdirectory(packaging)
endif()

################################################################################
