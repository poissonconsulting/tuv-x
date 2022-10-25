################################################################################
# Radiative transfer source

add_subdirectory(radiators)
add_subdirectory(solvers)

set(SRC radiative_transfer.F90
        radiator.F90
        radiator_factory.F90
        radiator_warehouse.F90
        solver.F90
        solver_factory.F90
        ${RADIATOR_SRC}
        ${SOLVER_SRC}
        )

list(TRANSFORM SRC PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")
set(RADIATIVE_TRANSFER_SRC ${SRC} PARENT_SCOPE)

################################################################################
