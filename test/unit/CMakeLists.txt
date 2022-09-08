################################################################################
# Test utilities

include(test_util)

################################################################################
# Add subdirectories

add_subdirectory(cross_section)
add_subdirectory(grid)
add_subdirectory(profile)
add_subdirectory(radiative_transfer)
add_subdirectory(radiator)
add_subdirectory(spectral_weight)

################################################################################
# TUV-x tests

create_standard_test(NAME la_sr_bands SOURCES test_ls_sr_bands.F90 test_utils.F90)

################################################################################