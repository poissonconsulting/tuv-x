! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The tuvx_spectral_wght_factory module

!> Build spectral weight calculators
module tuvx_spectral_wght_factory

  use tuvx_spectral_wght,                      only : spectral_wght_t
! use tuvx_uv_b_280_315_nm_spectral_wght_type, only : uv_b_280_315_nm_spectral_wght_t
! use tuvx_uv_b_280_320_nm_spectral_wght_type, only : uv_b_280_320_nm_spectral_wght_t
! use tuvx_spectral_wght_uv_a_315_400_nm,      only : spectral_wght_uv_a_315_400_nm_t
  use tuvx_spectral_wght_notch_filter,         only : spectral_wght_notch_filter_t
  use tuvx_spectral_wght_gaussian,             only : spectral_wght_gaussian_t
  use tuvx_spectral_wght_eppley,                           only : spectral_wght_eppley_t
! use tuvx_par_400_700nm_spectral_wght_type,               only : par_400_700nm_spectral_wght_t
! use tuvx_exponential_decay_spectral_wght_type,           only : exponential_decay_spectral_wght_t
! use tuvx_scup_mice_spectral_wght_type,                   only : scup_mice_spectral_wght_t
! use tuvx_standard_human_erythema_spectral_wght_type,     only : standard_human_erythema_spectral_wght_t
! use tuvx_uv_index_spectral_wght_type,                    only : uv_index_spectral_wght_t
! use tuvx_phytoplankton_boucher_spectral_wght_type,       only : phytoplankton_boucher_spectral_wght_t
! use tuvx_plant_damage_spectral_wght_type,                only : plant_damage_spectral_wght_t
! use tuvx_plant_damage_flint_caldwell_spectral_wght_type, only : plant_damage_flint_caldwell_spectral_wght_t
! use tuvx_plant_damage_flint_caldwell_ext_spectral_wght_type, only : plant_damage_flint_caldwell_ext_spectral_wght_t

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
      case( 'UV-B,280-315nm' )
!       new_spectral_wght_t => uv_b_280_315_nm_spectral_wght_t( config, grid_warehouse, profile_warehouse )
      case( 'UV-B*,280-320nm' )
!       new_spectral_wght_t => uv_b_280_320_nm_spectral_wght_t( config, grid_warehouse, profile_warehouse )
      case( 'UV-A,315-400nm' )
!       new_spectral_wght_t => spectral_wght_uv_a_315_400_nm_t( config, grid_warehouse, profile_warehouse )
      case( 'Notch Filter' )
        new_spectral_wght_t => spectral_wght_notch_filter_t( config, grid_warehouse, profile_warehouse )
      case( 'vis+,> 400nm' )
!       new_spectral_wght_t => visplus_spectral_wght_t( config, grid_warehouse, profile_warehouse )
      case( 'Gaussian' )
        new_spectral_wght_t => spectral_wght_gaussian_t( config, grid_warehouse, profile_warehouse )
      case( 'Gaussian, 320 nm, 10 nm FWHM' )
!       new_spectral_wght_t => gaussian_320_nm_10_nm_FWHM_spectral_wght_t( config, grid_warehouse, profile_warehouse )
      case( 'Gaussian, 340 nm, 10 nm FWHM' )
!       new_spectral_wght_t => gaussian_340_nm_10_nm_FWHM_spectral_wght_t( config, grid_warehouse, profile_warehouse )
      case( 'Gaussian, 380 nm, 10 nm FWHM' )
!       new_spectral_wght_t => gaussian_380_nm_10_nm_FWHM_spectral_wght_t( config, grid_warehouse, profile_warehouse )
      case( 'Eppley UV Photometer' )
        new_spectral_wght_t => spectral_wght_eppley_t( config, grid_warehouse, profile_warehouse )
      case( 'PAR, 400-700nm, umol m-2 s-1' )
!       new_spectral_wght_t => par_400_700nm_spectral_wght_t( config, grid_warehouse, profile_warehouse )
      case( 'Exponential decay' )
!       new_spectral_wght_t => exponential_decay_spectral_wght_t( config, grid_warehouse, profile_warehouse )
      case( 'SCUP-mice (de Gruijl et al., 1993)' )
!       new_spectral_wght_t => scup_mice_spectral_wght_t( config, grid_warehouse, profile_warehouse )
      case( 'Standard human erythema(Webb et al., 2011)' )
!       new_spectral_wght_t => standard_human_erythema_spectral_wght_t( config, grid_warehouse, profile_warehouse )
      case( 'UV Index(WMO, 1994; Webb et al., 2011)' )
!       new_spectral_wght_t => uv_index_spectral_wght_t( config, grid_warehouse, profile_warehouse )
      case( 'Phytoplankton (Boucher et al., 1994)' )
!       new_spectral_wght_t => phytoplankton_boucher_spectral_wght_t( config, grid_warehouse, profile_warehouse )
      case( 'Plant damage (Caldwell, 1971)' )
!       new_spectral_wght_t => plant_damage_spectral_wght_t( config, grid_warehouse, profile_warehouse )
      case( 'Plant damage,Flint&Caldwell,2003,orig' )
!       new_spectral_wght_t => plant_damage_flint_caldwell_spectral_wght_t( config, grid_warehouse, profile_warehouse )
      case( 'Plant damage,Flint&Caldwell,2003,ext390' )
!       new_spectral_wght_t => plant_damage_flint_caldwell_ext_spectral_wght_t( config, grid_warehouse, profile_warehouse )
      case default
        call die_msg( 450768215, "Invalid spectral wght type: '"//              &
                                 spectral_wght_type%to_char()//"'" )
    end select

  end function spectral_wght_builder

end module tuvx_spectral_wght_factory
