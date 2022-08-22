! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_profile_from_config
  ! Profile specified in json config file. See 
  ! :ref:`configuration-profiles-from-config` for more information.

  use musica_constants,                only : dk => musica_dk
  use tuvx_profile,                    only : profile_t

  implicit none

  private
  public :: profile_from_config_t

  type, extends(profile_t) :: profile_from_config_t
  contains
  end type profile_from_config_t

  interface profile_from_config_t
    module procedure constructor
  end interface profile_from_config_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse ) result( this )
    ! Initialize this profile

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t

    type(profile_from_config_t), pointer   :: this ! This f:type:`~tuvx_profile_from_config/profile_from_config_t`
    type(config_t), intent(inout)         :: config ! A profile config
    type(grid_warehouse_t), intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`

    ! Local variables
    character(len=*), parameter :: Iam = 'From config profile initialize: '
    integer                :: ndx
    real(dk)               :: uniformValue
    logical                :: found
    type(config_t)         :: grid_config
    type(string_t)         :: grid_name, grid_units
    class(grid_t), pointer :: theGrid
    type(string_t) :: required_keys(3), optional_keys(3)

    required_keys(1) = "type"
    required_keys(2) = "units"
    required_keys(3) = "grid"
    optional_keys(1) = "name"
    optional_keys(2) = "values"
    optional_keys(3) = "uniform value"

    call assert_msg( 530095878,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "profile from configuration file." )

    allocate( this )

    ! Get the handle
    call config%get( 'name', this%handle_, Iam, default = 'none' )
    call config%get( 'units', this%units_, Iam )

    call config%get( "grid", grid_config, Iam, found = found )
    call assert_msg( 376823788, found,                                      &
                      "Grid missing from profile configuration" )
    call grid_config%get( "name", grid_name, Iam )
    call grid_config%get( "units", grid_units, Iam )

    theGrid => grid_warehouse%get_grid( grid_name, grid_units )

    ! Get values from config file
    call config%get( "values", this%edge_val_, Iam, found = found )

    if( .not. found ) then
      call config%get( "uniform value", uniformValue, Iam, found = found )
      call assert_msg( 715232062, found,                                      &
                      "Neither 'values' or 'Uniform value' keyword specified" )

      this%edge_val_ = (/ ( uniformValue, ndx = 1, theGrid%ncells_ + 1 ) /)
    else
      ! Grid only required if we are using a constant value, but we should
      ! still perform a check if we are reading "values" from the config file
      ! to make sure that the "values" are the same length as the grid.
      ! each profile has to provide values at every element of 
      ! whatever grid they're on (height, wavelength, etc.)
      call assert_msg( 778098283,                                             &
        size( this%edge_val_, dim = 1 ) == size( theGrid%edge_, dim = 1 ),    &
        "The length of `values` must match the length of the specified grid")
    endif

    this%ncells_ = size( this%edge_val_ ) - 1

    this%mid_val_ = .5_dk * ( this%edge_val_( 1 : this%ncells_ ) +            &
                              this%edge_val_( 2 : this%ncells_ + 1 ) )

    this%delta_val_ = ( this%edge_val_( 2 : this%ncells_ + 1 ) -              &
                        this%edge_val_( 1 : this%ncells_ ) )


      deallocate( theGrid )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_profile_from_config
