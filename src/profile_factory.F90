! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_profile_factory
  ! Build :f:type:`~tuvx_profile/profile_t` objects

  implicit none

  private
  public :: profile_builder

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function profile_builder( config, grid_warehouse ) result( new_profile_t )
    ! Build an instance of a :f:type:`~tuvx_profile/profile_t`

    use musica_assert,                   only : die_msg
    use musica_config,                   only : config_t
    use musica_string,                   only : string_t
    use tuvx_profile_from_csv_file,      only : profile_from_csv_file_t
    use tuvx_profile_air,                only : profile_air_t
    use tuvx_profile_o2,                 only : profile_o2_t
    use tuvx_profile_o3,                 only : profile_o3_t
    use tuvx_profile_solar_zenith_angle, only : profile_solar_zenith_angle_t
    use tuvx_profile_earth_sun_distance, only : profile_earth_sun_distance_t
    use tuvx_profile,                    only : profile_t
    use tuvx_grid_warehouse,             only : grid_warehouse_t
    use tuvx_profile_from_config,        only : profile_from_config_t
    use tuvx_profile_surface_albedo,     only : profile_surface_albedo_t
    use tuvx_profile_extraterrestrial_flux, only : profile_extraterrestrial_flux_t

    type(config_t), intent(inout)         :: config ! Grid configuration data
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`

    class(profile_t), pointer :: new_profile_t ! New :f:type:`~tuvx_profile/profile_t` object

    ! Local variables
    character(len=*), parameter :: Iam = 'profile builder: '
    type(string_t) :: profile_type

    new_profile_t => null()
    call config%get( 'type', profile_type, Iam )

    select case( profile_type%to_char() )
      case( 'from csv file' )
        new_profile_t => profile_from_csv_file_t( config, grid_warehouse )
      case( 'extraterrestrial flux' )
        new_profile_t => profile_extraterrestrial_flux_t( config, grid_warehouse )
      case( 'from config file' )
        new_profile_t => profile_from_config_t( config, grid_warehouse )
      case( 'surface albedo' )
        new_profile_t => profile_surface_albedo_t( config, grid_warehouse )
      case( 'air' )
        new_profile_t => profile_air_t( config, grid_warehouse )
      case( 'O2' )
        new_profile_t => profile_o2_t( config, grid_warehouse )
      case( 'O3' )
        new_profile_t => profile_o3_t( config, grid_warehouse )
      case( 'solar zenith angle' )
        new_profile_t => profile_solar_zenith_angle_t( config, grid_warehouse )
      case( 'Earth-Sun distance' )
        new_profile_t => profile_earth_sun_distance_t( config, grid_warehouse )
      case default
        call die_msg( 884374015, "Invalid profile type: '" // &
          profile_type%to_char()//"'" )
    end select

  end function profile_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_profile_factory
