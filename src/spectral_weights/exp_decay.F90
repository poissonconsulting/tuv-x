! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_spectral_weight_exp_decay
  ! The exponential decay spectral weight type and related functions

  use tuvx_spectral_weight,            only : spectral_weight_t
  use musica_constants,                only : dk => musica_dk

  implicit none

  private
  public :: spectral_weight_exp_decay_t

  type, extends(spectral_weight_t) :: spectral_weight_exp_decay_t
    ! Calculator for exponential decay spectral weights
  contains
    procedure :: calculate => run
  end type spectral_weight_exp_decay_t

  !> Constructor
  interface spectral_weight_exp_decay_t
    module procedure constructor
  end interface spectral_weight_exp_decay_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result ( this )
    ! Initialize the spectral weight

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(spectral_weight_t),  pointer       :: this
    type(config_t),            intent(inout) :: config
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse
    type(profile_warehouse_t), intent(inout) :: profile_warehouse

    ! Local variables
    type(string_t) :: required_keys(1), optional_keys(1)

    required_keys(1) = "type"
    optional_keys(1) = "name"
    call assert_msg( 399197164,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "exp decay spectral wght." )

    allocate( spectral_weight_exp_decay_t :: this )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function run( this, grid_warehouse, profile_warehouse )                     &
      result( spectral_weight )
    ! Calculate the spectral weight

    use musica_string,                 only : string_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t
    use tuvx_grid,                     only : grid_t

    class(spectral_weight_exp_decay_t), intent(in)    :: this
    type(grid_warehouse_t),             intent(inout) :: grid_warehouse
    type(profile_warehouse_t),          intent(inout) :: profile_warehouse
    real(kind=dk), allocatable                        :: spectral_weight(:)

    ! Local variables

    class(grid_t), pointer      :: lambdaGrid => null( )

    lambdaGrid => grid_warehouse%get_grid( "wavelength", "nm" )

    allocate( spectral_weight( lambdaGrid%ncells_ ) )

    spectral_weight = 10._dk**( ( 300._dk - lambdaGrid%mid_ ) / 14._dk )

    if( associated( lambdaGrid ) ) deallocate( lambdaGrid )

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_spectral_weight_exp_decay
