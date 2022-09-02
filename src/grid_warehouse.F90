! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_grid_warehouse
! A warehouse to hold and distribute grids.

  use tuvx_grid, only : grid_ptr

  implicit none

  private
  public :: grid_warehouse_t

  !> Grid warehouse type
  type :: grid_warehouse_t
    private
    type(grid_ptr), allocatable :: grids_(:) ! grid objects
  contains
    !> get a copy of a grid object
    procedure, private :: get_grid_char, get_grid_string
    generic :: get_grid => get_grid_char, get_grid_string
    !> checks if a grid is present in the warehouse
    procedure, private :: exists_char, exists_string
    generic :: exists => exists_char, exists_string
    !> Returns the number of bytes required to pack the warehouse onto a buffer
    procedure :: pack_size
    !> Packs the warehouse onto a character buffer
    procedure :: mpi_pack
    !> Unpacks a warehouse from a character buffer into the object
    procedure :: mpi_unpack
    !> Finalize the object
    final :: finalize
  end type grid_warehouse_t

  !> Grid warehouse_t constructor
  interface grid_warehouse_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config ) result( grid_warehouse )
    ! Grid warehouse constructor

    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use musica_string,                 only : string_t
    use tuvx_grid_factory,             only : grid_builder

    type(config_t), intent(inout) :: config ! grid configuration data
    class(grid_warehouse_t), pointer :: grid_warehouse ! This :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`

    !> local variables
    character(len=*), parameter :: Iam = "Grid warehouse constructor: "
    type(config_t)              :: grid_config
    class(iterator_t), pointer  :: iter
    class(grid_warehouse_t), pointer :: grid_warehouse_ptr
    type(grid_ptr)              :: grid_obj
    character(len=32)           :: keychar
    type(string_t)              :: aswkey

    allocate( grid_warehouse )
    allocate( grid_warehouse%grids_(0) )

    ! iterate over grids
    iter => config%get_iterator()
    do while( iter%next() )
      keychar = config%key(iter)
      aswkey  = keychar
      call config%get( iter, grid_config, Iam )
      call grid_config%add( 'name', aswkey, Iam )

      ! Build grid objects
      grid_obj%val_ => grid_builder( grid_config )
      grid_warehouse%grids_ = [ grid_warehouse%grids_, grid_obj ]
    end do
    deallocate( iter )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_grid_char( this, name, units ) result( a_grid_ptr )
    ! Get copy of a grid object

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t

    class(grid_warehouse_t), intent(inout) :: this     ! This :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    character(len=*),        intent(in)    :: name     ! The name of a grid, see :ref:`configuration-grids` for grid names
    character(len=*),        intent(in)    :: units    ! The units of the grid
    class(grid_t), pointer                 :: a_grid_ptr ! The :f:type:`~tuvx_grid/grid_t` which matches the name passed in

    integer :: ndx
    logical :: found

    found = .false.
    do ndx = 1, size( this%grids_ )
      if( name .eq. this%grids_( ndx )%val_%handle_ ) then
        found = .true.
        exit
      endif
    end do
    call assert_msg( 345804219, found, "Invalid grid name: '"//name//"'" )
    call assert_msg( 509243577,                                               &
                     units .eq. this%grids_( ndx )%val_%units( ),             &
                     "Grid '"//name//"' has units of '"//                     &
                     this%grids_( ndx )%val_%units( )//"' not '"//            &
                     units//"' as requested." )
    allocate( a_grid_ptr, source = this%grids_( ndx )%val_ )

  end function get_grid_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_grid_string( this, name, units ) result( a_grid_ptr )
    ! Get a copy of a grid object

    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t

    class(grid_warehouse_t), intent(inout) :: this ! This :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(string_t),          intent(in)    :: name ! The name of a grid, see :ref:`configuration-grids` for grid names
    type(string_t),          intent(in)    :: units ! The units of the grid
    class(grid_t), pointer                 :: a_grid_ptr ! The :f:type:`~tuvx_grid/grid_t` which matches the name passed in

    a_grid_ptr => this%get_grid_char( name%to_char( ), units%to_char( ) )

  end function get_grid_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function exists_char( this, name, units ) result( exists )
    ! checks if a grid exists in the warehouse

    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t

    class(grid_warehouse_t), intent(inout) :: this ! This :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    character(len=*),        intent(in)    :: name ! The name of a grid, see :ref:`configuration-grids` for grid names
    character(len=*),        intent(in)    :: units ! The units of the grid

    integer :: ndx

    exists = .false.
    do ndx = 1, size( this%grids_ )
      if( name .eq. this%grids_( ndx )%val_%handle_ ) then
        exists = .true.
        exit
      endif
    end do

  end function exists_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function exists_string( this, name, units ) result( exists )
    ! checks if a grid exists in the warehouse

    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t

    class(grid_warehouse_t), intent(inout) :: this ! This :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(string_t),          intent(in)    :: name ! The name of a grid, see :ref:`configuration-grids` for grid names
    type(string_t),          intent(in)    :: units ! The units of the grid

    exists = this%exists_char( name%to_char( ), units%to_char( ) )

  end function exists_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! Returns the number of bytes required to pack the warehouse onto a buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack_size
    use musica_string,                 only : string_t
    use tuvx_grid_factory,             only : grid_type_name

    class(grid_warehouse_t), intent(in) :: this ! warehouse to be packed
    integer, optional,       intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: i_grid
    type(string_t) :: type_name

    call assert( 880887280, allocated( this%grids_ ) )
    pack_size = musica_mpi_pack_size( size( this%grids_ ), comm )
    do i_grid = 1, size( this%grids_ )
    associate( grid => this%grids_( i_grid )%val_ )
      type_name = grid_type_name( grid )
      pack_size = pack_size +                                                 &
                  type_name%pack_size( comm ) +                               &
                  grid%pack_size( comm )
    end associate
    end do
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, position, comm )
    ! Packs the warehouse onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack
    use musica_string,                 only : string_t
    use tuvx_grid_factory,             only : grid_type_name

    class(grid_warehouse_t), intent(in)    :: this      ! warehouse to be packed
    character,               intent(inout) :: buffer(:) ! memory buffer
    integer,                 intent(inout) :: position  ! current buffer position
    integer, optional,       intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_grid
    type(string_t) :: type_name

    prev_pos = position
    call assert( 607322234, allocated( this%grids_ ) )
    call musica_mpi_pack( buffer, position, size( this%grids_ ), comm )
    do i_grid = 1, size( this%grids_ )
    associate( grid => this%grids_( i_grid )%val_ )
      type_name = grid_type_name( grid )
      call type_name%mpi_pack( buffer, position, comm )
      call grid%mpi_pack(      buffer, position, comm )
    end associate
    end do
    call assert( 379779066, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! Unpacks a warehouse from a character buffer into the object

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack
    use musica_string,                 only : string_t
    use tuvx_grid_factory,             only : grid_allocate

    class(grid_warehouse_t), intent(out)   :: this      ! warehouse to be unpacked
    character,               intent(inout) :: buffer(:) ! memory buffer
    integer,                 intent(inout) :: position  ! current buffer position
    integer, optional,       intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_grid, n_grids
    type(string_t) :: type_name

    prev_pos = position
    call musica_mpi_unpack( buffer, position, n_grids, comm )
    if( allocated( this%grids_ ) ) deallocate( this%grids_ )
    allocate( this%grids_( n_grids ) )
    do i_grid = 1, n_grids
    associate( grid => this%grids_( i_grid )%val_ )
      call type_name%mpi_unpack( buffer, position, comm )
      grid => grid_allocate( type_name )
      call grid%mpi_unpack( buffer, position, comm )
    end associate
    end do
    call assert( 459962920, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Finalize grid warehouse

    !> Arguments
    type(grid_warehouse_t), intent(inout) :: this ! This :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`

    integer :: ndx

    if( allocated( this%grids_ ) ) then
      do ndx = 1, size( this%grids_ )
        if( associated( this%grids_( ndx )%val_ ) ) then
          deallocate( this%grids_( ndx )%val_ )
        end if
      end do
      deallocate( this%grids_ )
    endif

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_grid_warehouse
