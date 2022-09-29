################################################################################
# Photo-decomp tool source

set(TUVX_SRC constants.F90
             core.F90
             cross_section.F90
             cross_section_factory.F90
             cross_section_warehouse.F90
             diagnostic_util.F90
             dose_rates.F90
             grid.F90
             grid_factory.F90
             grid_warehouse.F90
             interpolate.F90
             la_sr_bands.F90
             linear_algebra.F90
             netcdf.F90
             photolysis_rates.F90
             profile.F90
             profile_factory.F90
             profile_warehouse.F90
             quantum_yield.F90
             quantum_yield_factory.F90
             spectral_weight.F90
             spectral_weight_factory.F90
             spherical_geometry.F90
             tuvx.F90
             util.F90
             )

add_subdirectory(cross_sections)
add_subdirectory(grids)
add_subdirectory(linear_algebras)
add_subdirectory(profiles)
add_subdirectory(quantum_yields)
add_subdirectory(radiative_transfer)
add_subdirectory(spectral_weights)

add_library(tuvxlib ${TUVX_SRC} ${GRID_SRC} ${PROFILE_SRC}
                    ${CROSS_SECTION_SRC} ${QUANTUM_YIELD_SRC} ${SPECTRAL_WGHT_SRC}
                    ${RADIATIVE_TRANSFER_SRC} ${LINEAR_ALGEBRA_SRC})
target_link_libraries(tuvxlib musicalib ${NETCDF_LIBS})
set_target_properties(tuvxlib PROPERTIES OUTPUT_NAME tuvx)

add_executable(tuv-x tuvx.F90)
target_link_libraries(tuv-x tuvxlib)

################################################################################
