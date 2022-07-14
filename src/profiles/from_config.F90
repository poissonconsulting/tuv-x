! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! profile specified in json config file
module tuvx_profile_from_config

  use musica_constants, only : &
    dk => musica_dk, ik => musica_ik, lk => musica_lk
  use tuvx_profile,     only : profile_t

  implicit none

  private
  public :: from_config_t

  type, extends(profile_t) :: from_config_t
  contains
  end type from_config_t

  !> Constructor
  interface from_config_t
    module procedure constructor
  end interface from_config_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize grid
  function constructor( profile_config, grid_warehouse ) result( this )

    use musica_config, only : config_t
    use musica_string, only : string_t
    use musica_assert, only : assert_msg
    use tuvx_grid,     only : grid_t
    use tuvx_grid_warehouse,  only : grid_warehouse_t

    type(from_config_t), pointer          :: this
    type(config_t), intent(inout)         :: profile_config
    type(grid_warehouse_t), intent(inout) :: grid_warehouse

    ! Local variables
    character(len=*), parameter :: Iam = 'From config profile initialize: '
    integer(ik)                 :: ndx
    real(dk)                    :: uniformValue
    logical(lk)                 :: found
    type(config_t)              :: grid_config
    type(string_t)              :: grid_name, grid_units
    class(grid_t), pointer :: theGrid

    allocate( this )

    ! Get the handle
    call profile_config%get( 'name', this%handle_, Iam, default = 'none' )

    ! Get values from config file
    call profile_config%get( "values", this%edge_val_, Iam, found=found )
    if( .not. found ) then
      call profile_config%get( "uniform value", uniformValue, Iam,            &
                               found = found )
      call assert_msg( 715232062, found,                                      &
                      "Neither 'Values' or 'Uniform value' keyword specified" )

      call profile_config%get( "grid", grid_config, Iam, found = found )
      call assert_msg( 376823788, found,                                      &
                       "Grid missing from profile configuration" )
      call grid_config%get( "name", grid_name, Iam )
      call grid_config%get( "units", grid_units, Iam )

      theGrid => grid_warehouse%get_grid( grid_name, grid_units )
      this%edge_val_ = (/ (uniformValue,ndx=1,theGrid%ncells_+1_ik) /)
    endif

    this%ncells_ = size(this%edge_val_) - 1_ik

    this%mid_val_ = .5_dk * ( &
      this%edge_val_(1_ik:this%ncells_) + &
      this%edge_val_(2_ik:this%ncells_+1_ik) &
    )

    this%delta_val_ = (this%edge_val_(2_ik:this%ncells_+1_ik) - &
      this%edge_val_(1_ik:this%ncells_))

    deallocate( theGrid )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_profile_from_config
