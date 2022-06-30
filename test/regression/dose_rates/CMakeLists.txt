################################################################################
# Test utilities

include(test_util)

################################################################################
# Dose rate regression tests

add_test(NAME regression_dose_rates COMMAND test/regression/dose_rates/dose_rates.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_test(NAME memcheck_regression_dose_rates COMMAND test/regression/dose_rates/dose_rates_memcheck.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

################################################################################
