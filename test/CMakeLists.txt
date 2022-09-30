################################################################################
# Test utilities

include(test_util)

################################################################################
# Copy test data

add_custom_target(copy_test_data ALL ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR}/test)

################################################################################
# Add subdirectories

add_subdirectory(unit)
if(NOT ${CMAKE_Fortran_COMPILER_ID} MATCHES "NAG")
  # oldtuv doesn't build with NAG, so bypass the regression tests
  add_subdirectory(oldtuv)
  add_subdirectory(regression)
endif()

################################################################################
# Run examples as tests

add_test(NAME full_example COMMAND tuv-x examples/full_config.json
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

################################################################################
