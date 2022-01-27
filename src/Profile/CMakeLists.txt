################################################################################
# One dimensional grid source

set(PROFILE_SRC profile/abstract.Profile.type.F90
                profile/from_csv_file.type.F90
                profile/from_config.type.F90
                profile/srfAlbedo_from_config.type.F90
                profile/air.from_csv_file.type.F90
                profile/o2.from_csv_file.type.F90
                profile/o3.from_csv_file.type.F90
                profile/sza_from_time.type.F90
                profile/earth_sun_distance.type.F90
                profile/Profile_factory.F90
                profile/Profile_warehouse.F90
                PARENT_SCOPE)

################################################################################
