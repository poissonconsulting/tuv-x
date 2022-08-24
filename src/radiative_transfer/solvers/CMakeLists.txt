################################################################################
# Radiative transfer solver source

set(SRC delta_eddington.F90
        )

list(TRANSFORM SRC PREPEND "solvers/")
set(SOLVER_SRC ${SRC} PARENT_SCOPE)

################################################################################
