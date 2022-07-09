! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The tuvx_spectral_wght_factory module

!> Build spectral weight calculators
module tuvx_spectral_wght_factory

  use tuvx_spectral_wght,                      only : spectral_wght_t
  use tuvx_spectral_wght_notch_filter,         only : spectral_wght_notch_filter_t
  use tuvx_spectral_wght_gaussian,             only : spectral_wght_gaussian_t
  use tuvx_spectral_wght_eppley,                           only : spectral_wght_eppley_t
  use tuvx_spectral_wght_par,                              only : spectral_wght_par_t
  use tuvx_spectral_wght_exp_decay,                        only : spectral_wght_exp_decay_t
  use tuvx_spectral_wght_scup_mice,                        only : spectral_wght_scup_mice_t
  use tuvx_spectral_wght_standard_human_erythema,          only : spectral_wght_standard_human_erythema_t
  use tuvx_spectral_wght_uv_index,                         only : spectral_wght_uv_index_t
  use tuvx_spectral_wght_phytoplankton_boucher,            only : spectral_wght_phytoplankton_boucher_t
  use tuvx_spectral_wght_plant_damage,                     only : spectral_wght_plant_damage_t
  use tuvx_spectral_wght_plant_damage_flint_caldwell,      only : spectral_wght_plant_damage_flint_caldwell_t
  use tuvx_spectral_wght_plant_damage_flint_caldwell_ext,  only : spectral_wght_plant_damage_flint_caldwell_ext_t

  implicit none

  private
  public :: spectral_wght_builder

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function spectral_wght_builder( config, grid_warehouse, profile_warehouse ) result( new_spectral_wght_t )

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use musica_constants,              only : dk => musica_dk
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    !> Arguments
    !> New rate constant calculator
    class(spectral_wght_t), pointer :: new_spectral_wght_t
    !> Spectral wght configuration
    type(config_t), intent(inout)   :: config
    !> Warehouses
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse
    type(profile_warehouse_t), intent(inout) :: profile_warehouse

    type(string_t) :: spectral_wght_type
    character(len=*), parameter :: Iam = 'spectral wght builder: '

    new_spectral_wght_t => null()
    call config%get( 'type', spectral_wght_type, Iam )

    select case( spectral_wght_type%to_char() )
      case( 'base' )
        new_spectral_wght_t => spectral_wght_t( config, grid_warehouse, profile_warehouse )
      case( 'Notch Filter' )
        new_spectral_wght_t => spectral_wght_notch_filter_t( config, grid_warehouse, profile_warehouse )
      case( 'Gaussian' )
        new_spectral_wght_t => spectral_wght_gaussian_t( config, grid_warehouse, profile_warehouse )
      case( 'Eppley UV Photometer' )
        new_spectral_wght_t => spectral_wght_eppley_t( config, grid_warehouse, profile_warehouse )
      case( 'PAR, 400-700 nm, umol m-2 s-1' )
        new_spectral_wght_t => spectral_wght_par_t( config, grid_warehouse, profile_warehouse )
      case( 'Exponential decay, 14 nm/10' )
        new_spectral_wght_t => spectral_wght_exp_decay_t( config, grid_warehouse, profile_warehouse )
      case( 'SCUP-mice (de Gruijl et al., 1993)' )
        new_spectral_wght_t => spectral_wght_scup_mice_t( config, grid_warehouse, profile_warehouse )
      case( 'Standard human erythema (Webb et al., 2011)' )
        new_spectral_wght_t => spectral_wght_standard_human_erythema_t( config, grid_warehouse, profile_warehouse )
      case( 'UV index (WMO, 1994; Webb et al., 2011)' )
        new_spectral_wght_t => spectral_wght_uv_index_t( config, grid_warehouse, profile_warehouse )
      case( 'Phytoplankton (Boucher et al., 1994)' )
        new_spectral_wght_t => spectral_wght_phytoplankton_boucher_t( config, grid_warehouse, profile_warehouse )
      case( 'Plant damage (Caldwell, 1971)' )
        new_spectral_wght_t => spectral_wght_plant_damage_t( config, grid_warehouse, profile_warehouse )
      case( 'Plant damage,Flint&Caldwell,2003,orig.' )
        new_spectral_wght_t => spectral_wght_plant_damage_flint_caldwell_t( config, grid_warehouse, profile_warehouse )
      case( 'Plant damage,Flint&Caldwell,2003,ext390' )
        new_spectral_wght_t => spectral_wght_plant_damage_flint_caldwell_ext_t( config, grid_warehouse, profile_warehouse )
      case default
        call die_msg( 450768215, "Invalid spectral wght type: '"//              &
                                 spectral_wght_type%to_char()//"'" )
    end select

  end function spectral_wght_builder

end module tuvx_spectral_wght_factory
