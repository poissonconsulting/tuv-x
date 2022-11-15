################################################################################
# Radiator transfer source

set(SRC aerosol.F90
        from_host.F90
        )

list(TRANSFORM SRC PREPEND "radiators/")
set(RADIATOR_SRC ${SRC} PARENT_SCOPE)

################################################################################
