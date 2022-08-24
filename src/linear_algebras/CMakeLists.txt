################################################################################
# Linear algebra source

set(SRC linpack.F90
        )

list(TRANSFORM SRC PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")
set(LINEAR_ALGEBRA_SRC ${SRC} PARENT_SCOPE)

################################################################################
