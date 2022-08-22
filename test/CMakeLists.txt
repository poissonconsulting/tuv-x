################################################################################
# Test utilities

include(test_util)

################################################################################
# Copy test data

add_custom_target(copy_test_data ALL ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR}/test)

################################################################################
# Add subdirectories

add_subdirectory(oldtuv)
add_subdirectory(unit)
add_subdirectory(regression)

################################################################################
