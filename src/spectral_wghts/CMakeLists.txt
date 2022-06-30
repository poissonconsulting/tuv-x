################################################################################
# Quantum yield source

set(SRC notch_filter.F90
        gaussian_filter.F90
        eppley.F90
        )

list(TRANSFORM SRC PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")
set(SPECTRAL_WGHT_SRC ${SRC} PARENT_SCOPE)

################################################################################
