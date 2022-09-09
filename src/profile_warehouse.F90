! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_profile_warehouse
  ! The profile warehouse type and related functions. Holds one or more
  ! :f:type:`~tuvx_profile/profile_t` s created by the
  ! :f:mod:`tuvx_profile_factory`

  use tuvx_profile, only : profile_ptr

  implicit none

  private
  public :: profile_warehouse_t

  type profile_warehouse_t
    private
    type(profile_ptr), allocatable :: profiles_(:)
  contains
    procedure, private :: get_profile_char, get_profile_string
    generic :: get_profile => get_profile_char, get_profile_string
    procedure :: exists_char, exists_string
    generic :: exists => exists_char, exists_string
    ! returns the number of bytes required to pack the warehouse onto a buffer
    procedure :: pack_size
    ! packs the warehouse onto a character buffer
    procedure :: mpi_pack
    ! unpacks a warehouse from a character buffer into the object
    procedure :: mpi_unpack
    final :: finalize
  end type profile_warehouse_t

  interface profile_warehouse_t
    ! profile warehouse_t constructor
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse ) result( profile_warehouse )
    ! profile warehouse constructor

    use musica_config,        only : config_t
    use musica_iterator,      only : iterator_t
    use musica_string,        only : string_t
    use tuvx_grid_warehouse,  only : grid_warehouse_t
    use tuvx_profile_factory, only : profile_builder

    type(config_t),             intent(inout) :: config ! profile configuration data
    type(grid_warehouse_t),     intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    class(profile_warehouse_t), pointer       :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`

    ! local variables
    type(config_t)              :: profile_config
    class(iterator_t), pointer  :: iter
    type(profile_ptr)           :: profile_obj
    character(len=32)           :: keychar
    type(string_t)              :: aswkey
    character(len=*), parameter :: Iam = "profile warehouse constructor: "
    class(profile_warehouse_t), pointer :: profile_warehouse_ptr

    allocate( profile_warehouse )
    allocate( profile_warehouse%profiles_(0) )

    ! iterate over profiles
    iter => config%get_iterator( )
    do while( iter%next( ) )
      keychar = config%key( iter )
      aswkey  = keychar
      call config%get( iter, profile_config, Iam )
      call profile_config%add( 'name', aswkey, Iam )

      ! Build profile objects
      profile_obj%val_ => profile_builder( profile_config, grid_warehouse )
      profile_warehouse%profiles_ =                                       &
          [ profile_warehouse%profiles_, profile_obj ]
    end do

    deallocate( iter )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_profile_char( this, name, units ) result( a_profile_ptr )
    ! Get copy of a profile object

    use musica_assert,     only : assert_msg
    use tuvx_profile,      only : profile_t

    class(profile_warehouse_t), intent(inout) :: this
    character(len=*),           intent(in)    :: name
    character(len=*),           intent(in)    :: units
    class(profile_t),           pointer       :: a_profile_ptr

    ! Local variables
    character(len=*), parameter :: Iam = 'profile warehouse get_profile: '
    integer :: ndx
    logical :: found

    found = .false.
    do ndx = 1, size( this%profiles_ )
      if( name .eq. this%profiles_( ndx )%val_%handle_ ) then
        found = .true.
        exit
      endif
    end do

    call assert_msg( 460768214, found, "Invalid profile handle: '"//name//"'" )
    call assert_msg( 465479557,                                               &
                     units .eq. this%profiles_( ndx )%val_%units( ),      &
                     "Profile '"//name//"' has units of '"//                  &
                     this%profiles_( ndx )%val_%units( )//"' not '"//     &
                     units//"' as requested." )
    allocate( a_profile_ptr, source = this%profiles_( ndx )%val_ )

  end function get_profile_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_profile_string( this, name, units ) result( a_profile_ptr )
    ! Get a copy of a profile object

    use musica_string,                 only : string_t
    use tuvx_profile,                  only : profile_t

    class(profile_warehouse_t), intent(inout) :: this
    type(string_t),             intent(in)    :: name
    type(string_t),             intent(in)    :: units
    class(profile_t), pointer                 :: a_profile_ptr

    a_profile_ptr => this%get_profile_char( name%to_char( ), units%to_char( ) )

  end function get_profile_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function exists_char( this, name, units ) result( exists )
    ! checks if a profile exists in the warehouse

    use musica_string,                 only : string_t
    use tuvx_profile,                  only : profile_t

    class(profile_warehouse_t), intent(inout) :: this
    character(len=*),           intent(in)    :: name
    character(len=*),           intent(in)    :: units

    integer :: ndx

    exists = .false.
    do ndx = 1, size( this%profiles_ )
      if( name .eq. this%profiles_( ndx )%val_%handle_ ) then
        exists  = .true.
        exit
      endif
    end do

  end function exists_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function exists_string( this, name, units ) result( exists )
    ! checks if a profile exists in the warehouse

    use musica_string,                 only : string_t
    use tuvx_profile,                  only : profile_t

    class(profile_warehouse_t), intent(inout) :: this
    type(string_t),             intent(in)    :: name
    type(string_t),             intent(in)    :: units

    exists = this%exists_char( name%to_char( ), units%to_char( ) )

  end function exists_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! returns the number of bytes required to pack the warehouse onto a buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack_size
    use musica_string,                 only : string_t
    use tuvx_profile_factory,          only : profile_type_name

    class(profile_warehouse_t), intent(in) :: this ! warehouse to be packed
    integer, optional,          intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: i_profile
    type(string_t) :: type_name

    call assert( 821308000, allocated( this%profiles_ ) )
    pack_size = musica_mpi_pack_size( size( this%profiles_ ) )
    do i_profile = 1, size( this%profiles_ )
    associate( profile => this%profiles_( i_profile )%val_ )
      type_name = profile_type_name( profile )
      pack_size = pack_size                                                   &
                  + type_name%pack_size( comm )                               &
                  + profile%pack_size(   comm )
    end associate
    end do
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, position, comm )
    ! packs the warehouse onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack
    use musica_string,                 only : string_t
    use tuvx_profile_factory,          only : profile_type_name

    class(profile_warehouse_t), intent(in)    :: this      ! warehouse to be packed
    character,                  intent(inout) :: buffer(:) ! memory buffer
    integer,                    intent(inout) :: position  ! current buffer position
    integer, optional,          intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_profile
    type(string_t) :: type_name

    prev_pos = position
    call assert( 333182583, allocated( this%profiles_ ) )
    call musica_mpi_pack( buffer, position, size( this%profiles_ ), comm )
    do i_profile = 1, size( this%profiles_ )
    associate( profile => this%profiles_( i_profile )%val_ )
      type_name = profile_type_name( profile )
      call type_name%mpi_pack( buffer, position, comm )
      call profile%mpi_pack(   buffer, position, comm )
    end associate
    end do
    call assert( 777643951, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! unpacks a warehouse from a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack
    use musica_string,                 only : string_t
    use tuvx_profile_factory,          only : profile_allocate

    class(profile_warehouse_t), intent(out)   :: this      ! warehouse to be unpacked
    character,                  intent(inout) :: buffer(:) ! memory buffer
    integer,                    intent(inout) :: position  ! current buffer position
    integer, optional,          intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_profile, n_profiles
    type(string_t) :: type_name

    prev_pos = position
    call musica_mpi_unpack( buffer, position, n_profiles, comm )
    if( allocated( this%profiles_ ) ) deallocate( this%profiles_ )
    allocate( this%profiles_( n_profiles ) )
    do i_profile = 1, n_profiles
    associate( profile => this%profiles_( i_profile )%val_ )
      call type_name%mpi_unpack( buffer, position, comm )
      profile => profile_allocate( type_name )
      call profile%mpi_unpack( buffer, position, comm )
    end associate
    end do
    call assert( 294782841, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Finalize profile warehouse

    use musica_constants, only : ik => musica_ik

    type(profile_warehouse_t), intent(inout) :: this

    ! Local variables
    integer(kind=ik) :: ndx
    character(len=*), parameter :: Iam = 'profile warehouse finalize: '

    if( allocated( this%profiles_ ) ) then
      do ndx = 1, size( this%profiles_ )
        if( associated( this%profiles_( ndx )%val_ ) ) then
          deallocate( this%profiles_( ndx )%val_ )
        end if
      end do
      deallocate( this%profiles_ )
    endif

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_profile_warehouse
