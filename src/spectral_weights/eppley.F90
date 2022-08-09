! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_spectral_weight_eppley
  ! The eppley type and related functions

  use tuvx_spectral_weight,            only : spectral_weight_t
  use musica_constants,                only : dk => musica_dk

  implicit none

  private
  public :: spectral_weight_eppley_t

  type, extends(spectral_weight_t) :: spectral_weight_eppley_t
    ! Calculator for Eppley spectral weight
  contains
    procedure :: calculate => run
  end type spectral_weight_eppley_t

  !> Constructor
  interface spectral_weight_eppley_t
    module procedure constructor
  end interface spectral_weight_eppley_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result ( this )
    ! Initialize the Epply spectral weight

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_spectral_weight,          only : base_constructor
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(spectral_weight_t),  pointer       :: this
    type(config_t),            intent(inout) :: config
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse
    type(profile_warehouse_t), intent(inout) :: profile_warehouse

    allocate( spectral_weight_eppley_t :: this )
    call base_constructor( this, config, grid_warehouse, profile_warehouse )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function run( this, grid_warehouse, profile_warehouse )                     &
      result( spectral_weight )
    ! Calculate the spectral weight

    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(spectral_weight_eppley_t), intent(in)    :: this
    type(grid_warehouse_t),          intent(inout) :: grid_warehouse
    type(profile_warehouse_t),       intent(inout) :: profile_warehouse
    real(kind=dk), allocatable                     :: spectral_weight(:)

    ! Local variables
    real(dk), parameter         :: NINETY = 90.0_dk
    real(dk)                    :: accum

    class(grid_t), pointer      :: lambdaGrid => null()

    lambdaGrid => grid_warehouse%get_grid( "wavelength", "nm" )

    spectral_weight = this%spectral_weight_parms(1)%array( :, 1 )
    accum = dot_product( spectral_weight, lambdaGrid%delta_ )
    spectral_weight = NINETY * spectral_weight / accum

    if( associated( lambdaGrid ) ) deallocate( lambdaGrid )

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_spectral_weight_eppley
