################################################################################
# Copy test data

add_custom_target(copy_test_data ALL ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR}/test)

################################################################################
# Add subdirectories

add_subdirectory(oldtuv)

################################################################################
# Link photolib to a test and add it to the suite

macro(add_std_test test_name)
  target_include_directories(${test_name} PUBLIC ${CMAKE_BINARY_DIR}/src)
  target_link_libraries(${test_name} tuvxlib)
  add_test(NAME ${test_name} COMMAND ${test_name})
endmacro(add_std_test)

################################################################################
# Photo-decomp tests

add_executable(grid_test unit/grid/test.F90)
add_std_test(grid_test)

add_executable(profile_test Profiles/test.F90)
add_std_test(profile_test)

add_executable(radXfer_test radXfer/radXfer_core.F90 radXfer/radXfer_test.F90)
add_std_test(radXfer_test)

add_executable(radiator_test radiator/radiator_core.F90 radiator/radiator_test.F90)
add_std_test(radiator_test)

add_test(NAME regression_photolysis_rates COMMAND test/regression/photolysis_rate/photolysis_rates.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_test(NAME regression_aerosols_radiator COMMAND test/regression/radiation.aerosols.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME regression_o2_radiator       COMMAND test/regression/radiation.o2.sh       WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME regression_o3_radiator       COMMAND test/regression/radiation.o3.sh       WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME regression_rayleigh_radiator COMMAND test/regression/radiation.rayleigh.sh WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME regression_all_radiators     COMMAND test/regression/radiation.all.sh      WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(config_test config_tst/config.F90 config_tst/test.F90)
add_std_test(config_test)

################################################################################