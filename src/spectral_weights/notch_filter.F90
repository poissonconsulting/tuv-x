! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_spectral_weight_notch_filter
  ! The notch filter spectral weight type and related functions

  use tuvx_spectral_weight,            only : spectral_weight_t
  use musica_constants,                only : dk => musica_dk

  implicit none

  private
  public :: spectral_weight_notch_filter_t

  type, extends(spectral_weight_t) :: spectral_weight_notch_filter_t
    ! Calculator for notch filter spectral weight
    real(dk) :: notch_filter_begin  ! Notch filter lower end point
    real(dk) :: notch_filter_end    ! Notch filter upper end point
  contains
    procedure :: calculate => run
  end type spectral_weight_notch_filter_t

  interface spectral_weight_notch_filter_t
    module procedure constructor
  end interface spectral_weight_notch_filter_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result ( this )
    ! Initialize the notch fileter spectral weight

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
    character(len=*), parameter :: Iam =                                      &
        'Notch filter spectral weight constrcutor'
    type(string_t) :: required_keys(3), optional_keys(1)

    required_keys(1) = "type"
    required_keys(2) = "notch filter begin"
    required_keys(3) = "notch filter end"
    optional_keys(1) = "name"
    call assert_msg( 399197574,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "notch_filter spectral wght." )

    allocate( spectral_weight_notch_filter_t :: this )
    select type( this )
      class is( spectral_weight_notch_filter_t )
        call config%get( 'notch filter begin', this%notch_filter_begin, Iam )
        call config%get( 'notch filter end', this%notch_filter_end, Iam )
    end select

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function run( this, grid_warehouse, profile_warehouse )                     &
      result( spectral_weight )
    ! Calculate the Notch filter spectral weight

    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(spectral_weight_notch_filter_t), intent(in)    :: this
    type(grid_warehouse_t),                intent(inout) :: grid_warehouse
    type(profile_warehouse_t),             intent(inout) :: profile_warehouse
    real(kind=dk), allocatable                           :: spectral_weight(:)

    ! Local variables
    real(dk), parameter         :: rZERO = 0.0_dk
    real(dk), parameter         :: rONE  = 1.0_dk

    class(grid_t), pointer      :: lambdaGrid => null()

    lambdaGrid => grid_warehouse%get_grid( "wavelength", "nm" )

    allocate( spectral_weight( lambdaGrid%ncells_ ) )

    where( this%notch_filter_begin < lambdaGrid%mid_ .and. &
           lambdaGrid%mid_ < this%notch_filter_end )
      spectral_weight = rONE
    elsewhere
      spectral_weight = rZERO
    endwhere

    if( associated( lambdaGrid ) ) deallocate( lambdaGrid )

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_spectral_weight_notch_filter
