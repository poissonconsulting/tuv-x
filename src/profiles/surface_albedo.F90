! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_profile_surface_albedo
  ! Surface albedo profile type

  use musica_constants,                only : dk => musica_dk
  use tuvx_profile,                    only : profile_t

  implicit none

  private
  public :: profile_surface_albedo_t

  type, extends(profile_t) :: profile_surface_albedo_t
  contains
  end type profile_surface_albedo_t

  interface profile_surface_albedo_t
    module procedure constructor
  end interface profile_surface_albedo_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse ) result( this )
    ! Initialize grid

    use musica_assert,                 only : assert_msg, die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t

    type(profile_surface_albedo_t), pointer   :: this ! This f:type:`~tuvx_profile_surface_albedo/profile_surface_albedo_t`
    type(config_t), intent(inout)         :: config ! A profile config
    type(grid_warehouse_t), intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`

    ! Local variables
    real(dk)                    :: uniformValue
    integer                     :: ndx
    class(grid_t), pointer      :: lambdaGrid
    character(len=*), parameter :: Iam = 'Surface albedo profile initialize: '
    type(string_t) :: required_keys(3), optional_keys(1)

    required_keys(1) = "type"
    required_keys(2) = "units"
    required_keys(3) = "uniform value"
    optional_keys(1) = "name"

    call assert_msg( 423806855,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "surface albedo profile." )

    allocate( this )

    ! Get the handle
    call config%get( 'name', this%handle_, Iam, default = 'none' )
    call config%get( 'units', this%units_, Iam )

    lambdaGrid => grid_warehouse%get_grid( "wavelength", "nm" )

    ! Get values from config file
    call config%get( "uniform Value", uniformValue, Iam )

    this%ncells_ = lambdaGrid%ncells_

    this%edge_val_ = (/ ( uniformValue, ndx = 1, this%ncells_ + 1 ) /)

    this%mid_val_ = .5_dk * ( this%edge_val_( 1 : this%ncells_ ) +            &
                              this%edge_val_( 2 : this%ncells_ + 1 ) )

    this%delta_val_ = this%edge_val_( 2 : this%ncells_ + 1 ) -                &
                      this%edge_val_( 1 : this%ncells_ )

    deallocate( lambdaGrid )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_profile_surface_albedo
