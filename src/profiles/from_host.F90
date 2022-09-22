! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_profile_from_host
  ! Profile whose values will be provided by the host application at runtime
  ! See :ref:`configuration-profiles-from-host` for more information.

  use musica_constants,                only : dk => musica_dk
  use tuvx_profile,                    only : profile_t

  implicit none

  public :: profile_from_host_t, profile_updater_t

  type, extends(profile_t) :: profile_from_host_t
    ! profile that can be updated from a host application
  contains
  end type profile_from_host_t

  ! profile constructor
  interface profile_from_host_t
    module procedure :: constructor_char
    module procedure :: constructor_string
  end interface profile_from_host_t

  type :: profile_updater_t
    private
    ! updater for `profile_from_host_t` profiles
    type(profile_from_host_t), pointer :: profile_ => null( )
  contains
    procedure :: update
  end type profile_updater_t

  ! profile updater constructor
  interface profile_updater_t
    module procedure :: updater_constructor
  end interface profile_updater_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor_char( name, units, number_of_cells ) result( this )
    ! Initializes the profile

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : to_char

    character(len=*),          intent(in) :: name            ! name of the profile
    character(len=*),          intent(in) :: units           ! units for the profile
    integer,                   intent(in) :: number_of_cells ! number of discreet profile cells
    type(profile_from_host_t), pointer    :: this            ! constructor profile

    allocate( this )

    this%handle_ = name
    this%units_  = units
    this%ncells_ = number_of_cells

    call assert_msg( 253742246, this%ncells_ >= 0,                            &
                     "Invalid profile size for profile from host: '"//        &
                     trim( to_char( number_of_cells ) ) )
    allocate( this%mid_val_(    this%ncells_     ) )
    allocate( this%edge_val_(   this%ncells_ + 1 ) )
    allocate( this%delta_val_(  this%ncells_     ) )
    allocate( this%layer_dens_( this%ncells_     ) )
    this%mid_val_(:)    = 0.0_dk
    this%edge_val_(:)   = 0.0_dk
    this%delta_val_(:)  = 0.0_dk
    this%layer_dens_(:) = 0.0_dk

  end function constructor_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor_string( name, units, number_of_cells ) result( this )
    ! Initializes the profile

    use musica_string,                 only : string_t

    type(string_t),            intent(in) :: name            ! name of the profile
    type(string_t),            intent(in) :: units           ! units for the profile
    integer,                   intent(in) :: number_of_cells ! number of discreet profile cells
    type(profile_from_host_t), pointer    :: this            ! constructed profile

    this => constructor_char( name%to_char( ), units%to_char( ),              &
                              number_of_cells )

  end function constructor_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function updater_constructor( profile ) result( this )
    ! Constructs an updater for a `profile_from_host_t` profile

    class(profile_from_host_t), target :: profile ! profile to be updated
    type(profile_updater_t)            :: this    ! new updater

    this%profile_ => profile

  end function updater_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update( this, mid_point_values, edge_values, layer_densities )
    ! Updates the target profile

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : to_char

    class(profile_updater_t), intent(inout) :: this                ! profile updater
    real(kind=dk), optional,  intent(in)    :: mid_point_values(:) ! new mid-point values
    real(kind=dk), optional,  intent(in)    :: edge_values(:)      ! new edge values
    real(kind=dk), optional,  intent(in)    :: layer_densities(:)  ! new layer densities

    integer :: size_profile, size_host

    call assert_msg( 314048995, associated( this%profile_ ),                  &
                     "Cannot update an unspecified profile" )
    if( present( mid_point_values ) ) then
      size_profile = size( this%profile_%mid_val_ )
      size_host    = size( mid_point_values )
      call assert_msg( 307331449, size_profile == size_host,                  &
                       "Size mismatch for profile mid-point values for "//    &
                       "profile '"//this%profile_%handle_//"'. Expected "//   &
                       trim( to_char( size_profile ) )//", got "//            &
                       trim( to_char( size_host ) ) )
      this%profile_%mid_val_(:) = mid_point_values(:)
    end if
    if( present( edge_values ) ) then
      size_profile = size( this%profile_%edge_val_ )
      size_host    = size( edge_values )
      call assert_msg( 573012833, size_profile == size_host,                  &
                       "Size mismatch for profile edge values for "//         &
                       "profile '"//this%profile_%handle_//"'. Expected "//   &
                       trim( to_char( size_profile ) )//", got "//            &
                       trim( to_char( size_host ) ) )
      this%profile_%edge_val_(:) = edge_values(:)
      this%profile_%delta_val_(:) =                                           &
          this%profile_%edge_val_( 2 : this%profile_%ncells_ + 1 ) -          &
          this%profile_%edge_val_( 1 : this%profile_%ncells_ )
    end if
    if( present( layer_densities ) ) then
      size_profile = size( this%profile_%layer_dens_ )
      size_host    = size( layer_densities )
      call assert_msg( 229340252, size_profile == size_host,                  &
                       "Size mismatch for profile layer densities for "//     &
                       "profile '"//this%profile_%handle_//"'. Expected "//   &
                       trim( to_char( size_profile ) )//", got "//            &
                       trim( to_char( size_host ) ) )
      this%profile_%layer_dens_(:) = layer_densities(:)
    end if

  end subroutine update

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_profile_from_host
