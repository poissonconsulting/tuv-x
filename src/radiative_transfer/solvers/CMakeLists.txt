################################################################################
# Radiative transfer solver source

set(SRC delta_eddington.F90
        discrete_ordinant.F90
        discrete_ordinant_util.F90
        )

list(TRANSFORM SRC PREPEND "solvers/")
set(SOLVER_SRC ${SRC} PARENT_SCOPE)

################################################################################
