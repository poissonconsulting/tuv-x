################################################################################
# Photo-decomp tool source

set(PHOTO_SRC abs_linalgebra.f90
              abstract_radXfer.f90
              delta_eddington.f90
              disord.f90
              disord_subs.f
              functs.f
              grids.f
              la_srb.f
              lepoly.f
              linalgebra.f90
              newlst.f
              numer.f
              odo3.f
              odrl.f
              orbit.f
              psndo.f
              qys.f
              rdetfl.f
              rdinp.f
              rdxs.f
              rtlink.f
              rtrans.f
              rxn.f
              savout.f
              setaer.f
              setalb.f
              setcld.f
              seth2o.f
              setmtx.f
              setno2.f
              seto2.f
              setsnw.f
              setso2.f
              sphers.f
              solvec.f
              swbiol.f
              swchem.f
              swdom.f
              swphys.f
              terint.f
              tuv_params.f
              vpair.f
              vpo3.f
              vptmp.f
              waters.f
              wrflut.f
              wshift.f
              zeroal.f
              )

add_subdirectory(cross_section)
add_subdirectory(quantum_yield)
add_subdirectory(radXfer_cross_section)
add_subdirectory(spectral_wght)

add_library(photolib ${PHOTO_SRC} ${CROSS_SRC} ${QY_SRC} ${RAD_SRC} ${SW_SRC})
target_link_libraries(photolib musicalib ${NETCDF_LIBS})
set_target_properties(photolib PROPERTIES OUTPUT_NAME photo)

add_executable(photo tuv.f)
target_link_libraries(photo photolib)

################################################################################
