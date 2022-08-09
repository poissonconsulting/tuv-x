! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_spectral_weight_gaussian
  ! The gaussian spectral weight type and related functions

  use tuvx_spectral_weight,            only : spectral_weight_t
  use musica_constants,                only : dk => musica_dk

  implicit none

  private
  public :: spectral_weight_gaussian_t

  type, extends(spectral_weight_t) :: spectral_weight_gaussian_t
    ! Calculator for Gaussian spectral weight
    real(dk) :: centroid ! The gaussian centroid
  contains
    procedure :: calculate => run
  end type spectral_weight_gaussian_t

  interface spectral_weight_gaussian_t
    module procedure constructor
  end interface spectral_weight_gaussian_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result ( this )
    ! Initialize the Gaussian spectral weight

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

    ! Local variables
    character(len=*), parameter :: Iam = 'Gaussian spectral weight constructor'
    type(string_t) :: required_keys(2), optional_keys(1)

    required_keys(1) = "type"
    required_keys(2) = "centroid"
    optional_keys(1) = "name"
    call assert_msg( 399197574,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "gaussian spectral wght." )

    allocate( spectral_weight_gaussian_t :: this )
    select type( this )
      class is( spectral_weight_gaussian_t )
        call config%get( 'centroid', this%centroid, Iam )
    end select

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function run( this, grid_warehouse, profile_warehouse )                     &
      result( spectral_weight )
    ! Calculate the Gaussian spectral weight

    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(spectral_weight_gaussian_t), intent(in)    :: this
    type(grid_warehouse_t),            intent(inout) :: grid_warehouse
    type(profile_warehouse_t),         intent(inout) :: profile_warehouse
    real(kind=dk), allocatable                       :: spectral_weight(:)

    ! Local variables
    real(dk), parameter         :: rTWO = 2.0_dk
    real(dk)                    :: accum

    class(grid_t), pointer      :: lambdaGrid => null()

    lambdaGrid => grid_warehouse%get_grid( "wavelength", "nm" )

    spectral_weight = exp( -( log( rTWO ) * .04_dk                            &
                      * ( lambdaGrid%mid_ - this%centroid )**2 ) )
    accum = sum( spectral_weight )
    spectral_weight = spectral_weight / accum

    if( associated( lambdaGrid ) ) deallocate( lambdaGrid )

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_spectral_weight_gaussian
