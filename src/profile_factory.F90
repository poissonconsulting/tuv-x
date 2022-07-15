! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The tuvx_profile_factory module

!> Build profile objects
module tuvx_profile_factory


  implicit none

  private
  public :: profile_builder

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function profile_builder( config, grid_warehouse ) result( new_profile_t )

    use musica_assert,                   only : die_msg
    use musica_config,                   only : config_t
    use musica_string,                   only : string_t
    use tuvx_profile_from_csv_file,      only : from_csv_file_t
    use tuvx_profile_air,                only : profile_air_t
    use tuvx_profile_o2,                 only : o2_from_csv_file_t
    use tuvx_profile_o3,                 only : o3_from_csv_file_t
    use tuvx_profile_solar_zenith_angle, only : sza_from_time_t
    use tuvx_profile_earth_sun_distance, only : earth_sun_distance_t
    use tuvx_profile,                    only : profile_t
    use tuvx_grid_warehouse,             only : grid_warehouse_t
    use tuvx_profile_from_config,        only : from_config_t
    use tuvx_profile_surface_albedo,     only : surface_albedo_t
    use tuvx_profile_extraterrestrial_flux, only : extraterrestrial_flux_t

    !> Grid configuration data
    type(config_t), intent(inout)         :: config
    type(grid_warehouse_t), intent(inout) :: grid_warehouse

    !> New profile object
    class(profile_t), pointer :: new_profile_t

    ! Local variables
    character(len=*), parameter :: Iam = 'profile builder: '
    type(string_t) :: profile_type

    new_profile_t => null()
    call config%get( 'type', profile_type, Iam )

    select case( profile_type%to_char() )
      case( 'from csv file' )
        new_profile_t => from_csv_file_t( config, grid_warehouse )
      case( 'extraterrestrial flux' )
        new_profile_t => extraterrestrial_flux_t( config, grid_warehouse )
      case( 'from config file' )
        new_profile_t => from_config_t( config, grid_warehouse )
      case( 'surface albedo' )
        new_profile_t => surface_albedo_t( config, grid_warehouse )
      case( 'air' )
        new_profile_t => profile_air_t( config, grid_warehouse )
      case( 'O2' )
        new_profile_t => o2_from_csv_file_t( config, grid_warehouse )
      case( 'O3' )
        new_profile_t => o3_from_csv_file_t( config, grid_warehouse )
      case( 'solar zenith angle' )
        new_profile_t => sza_from_time_t( config, grid_warehouse )
      case( 'Earth-Sun distance' )
        new_profile_t => earth_sun_distance_t( config, grid_warehouse )
      case default
        call die_msg( 884374015, "Invalid profile type: '" // &
          profile_type%to_char()//"'" )
    end select

  end function profile_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_profile_factory
