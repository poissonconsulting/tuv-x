! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_cross_section_warehouse
! A type to store all cross sections needed for a particular run
!

  use musica_constants,                only : musica_dk
  use musica_string,                   only : string_t
  use tuvx_cross_section,              only : cross_section_ptr

  implicit none

  private
  public :: cross_section_warehouse_t

  !> Radiative tranfser cross section type
  type cross_section_warehouse_t
    private
    type(cross_section_ptr), allocatable :: cross_sections_(:) ! A:f:type:`~tuvx_cross_section/cross_section_ptr`
    type(string_t), allocatable          :: handles_(:) ! cross section "handle"
  contains
    !> Get a copy of a specific cross section
    procedure :: get
    !> Returns the number of bytes required to pack the warehouse onto a buffer
    procedure :: pack_size
    !> Packs the warehouse onto a character buffer
    procedure :: mpi_pack
    !> Unpacks the warehouse from a character buffer into the object
    procedure :: mpi_unpack
    !> Finalize the object
    final :: finalize
  end type cross_section_warehouse_t

  !> cross_section_warehouse_t constructor
  interface cross_section_warehouse_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result( new_obj )
    ! Constructor of cross_section_warehouse_t objects

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use tuvx_cross_section_factory,    only : cross_section_builder
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : Profile_warehouse_t

    class(cross_section_warehouse_t), pointer :: new_obj ! New radiative :f:type:`~tuvx_cross_section_warehouse/cross_section_warehouse_t` object
    type(config_t),             intent(inout) :: config ! Cross section configuration data
    type(grid_warehouse_t),     intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(Profile_warehouse_t),  intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`

    ! local variables
    integer :: ndx
    character(len=*), parameter :: Iam = "Cross section constructor"
    type(config_t) :: cross_section_config
    class(iterator_t), pointer :: iter
    type(cross_section_ptr) :: cross_section_ptr
    character(len=32)           :: keychar
    type(string_t)              :: cross_section_name
    type(string_t), allocatable :: netcdfFiles(:)
    logical                     :: found

    allocate( new_obj )
    allocate( string_t :: new_obj%handles_(0) )
    allocate( new_obj%cross_sections_(0) )

    ! iterate over cross sections
    iter => config%get_iterator( )
    do while( iter%next( ) )
      call config%get( iter, cross_section_config, Iam )

      ! save the cross section name for lookups
      call cross_section_config%get( "name", cross_section_name, Iam )
      new_obj%handles_ = [ new_obj%handles_, cross_section_name ]

      ! build and store cross section object
      cross_section_ptr%val_ => cross_section_builder( cross_section_config,  &
                                           grid_warehouse, profile_warehouse )
      new_obj%cross_sections_ =                                         &
          [ new_obj%cross_sections_, cross_section_ptr ]
    end do
    deallocate( iter )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get( this, cross_section_name ) result( cross_section_ptr )
    ! Get a copy of a specific radiative transfer cross section object

    use musica_assert,                 only : die_msg
    use musica_string,                 only : string_t
    use tuvx_cross_section,            only : cross_section_t

    class(cross_section_t),           pointer       :: cross_section_ptr ! Pointer to a copy of the requested :f:type:`~tuvx_cross_section/cross_section_t`
    class(cross_section_warehouse_t), intent(inout) :: this ! A :f:type:`~tuvx_cross_section_warehouse/cross_section_warehouse_t`
    type(string_t),                   intent(in)    :: cross_section_name ! Name of the cross section to find

    ! Local variables
    character(len=*), parameter :: Iam = 'radXfer cross section warehouse get: '
    integer :: ndx
    logical :: found

    found = .false.
    do ndx = 1, size( this%handles_ )
      if( cross_section_name .eq. this%handles_( ndx ) ) then
        found = .true.
        exit
      endif
    end do

    if( .not. found ) then
      call die_msg( 980636382, "Invalid cross_section_name: '"//              &
                               cross_section_name%to_char()//"'" )
    endif
    allocate( cross_section_ptr,                                              &
              source = this%cross_sections_( ndx )%val_ )

  end function get

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! Returns the number of bytes required to pack the warehouse onto a buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack_size
    use tuvx_cross_section_factory,    only : cross_section_type_name

    class(cross_section_warehouse_t), intent(in) :: this ! warehouse to be packed
    integer,                          intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: i_cross_section
    type(string_t) :: type_name

    call assert( 161062144, allocated( this%cross_sections_ ) )
    call assert( 385698834, allocated( this%handles_ ) )
    call assert( 259945397,                                                   &
                 size( this%cross_sections_ ) .eq. size( this%handles_ ) )
    pack_size = musica_mpi_pack_size( size( this%cross_sections_ ), comm )
    do i_cross_section = 1, size( this%cross_sections_ )
    associate( cross_section => this%cross_sections_( i_cross_section )%val_ )
      type_name = cross_section_type_name( cross_section )
      pack_size = pack_size +                                                 &
                  type_name%pack_size( comm ) +                               &
                  cross_section%pack_size( comm ) +                           &
                  this%handles_( i_cross_section )%pack_size( comm )
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
    use tuvx_cross_section_factory,    only : cross_section_type_name

    class(cross_section_warehouse_t), intent(in) :: this ! warehouse to be packed
    character,                        intent(inout) :: buffer(:) ! memory buffer
    integer,                          intent(inout) :: position  ! current buffer position
    integer,                          intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_cross_section
    type(string_t) :: type_name

    prev_pos = position
    call assert( 442823572, allocated( this%cross_sections_ ) )
    call assert( 386890547, allocated( this%handles_        ) )
    call assert( 264043588,                                                   &
                 size( this%cross_sections_ ) .eq. size( this%handles_ ) )
    call musica_mpi_pack( buffer, position, size( this%cross_sections_ ),     &
                          comm )
    do i_cross_section = 1, size( this%cross_sections_ )
    associate( cross_section => this%cross_sections_( i_cross_section )%val_ )
      type_name = cross_section_type_name( cross_section )
      call type_name%mpi_pack( buffer, position, comm )
      call cross_section%mpi_pack( buffer, position, comm )
      call this%handles_( i_cross_section )%mpi_pack( buffer, position, comm )
    end associate
    end do
    call assert( 359568068, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! Unpacks the warehouse from a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack
    use tuvx_cross_section_factory,    only : cross_section_allocate

    class(cross_section_warehouse_t), intent(out)   :: this      ! warehouse to be unpacked
    character,                        intent(inout) :: buffer(:) ! memory buffer
    integer,                          intent(inout) :: position  ! current buffer position
    integer,                          intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_cross_section, n_cross_sections
    type(string_t) :: type_name

    prev_pos = position
    call musica_mpi_unpack( buffer, position, n_cross_sections, comm )
    if( allocated( this%cross_sections_ ) ) deallocate( this%cross_sections_ )
    if( allocated( this%handles_        ) ) deallocate( this%handles_        )
    allocate( this%cross_sections_( n_cross_sections ) )
    allocate( this%handles_(        n_cross_sections ) )
    do i_cross_section = 1, n_cross_sections
    associate( cross_section => this%cross_sections_( i_cross_section )%val_ )
      call type_name%mpi_unpack( buffer, position, comm )
      cross_section => cross_section_allocate( type_name )
      call cross_section%mpi_unpack( buffer, position, comm )
      call this%handles_( i_cross_section )%mpi_unpack( buffer, position, comm )
    end associate
    end do
    call assert( 659576600, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Finalize the cross section warehouse

    type(cross_section_warehouse_t), intent(inout) :: this ! This :f:type:`~tuvx_cross_section_warehouse/cross_section_warehouse_t`

    integer :: ndx

    if( allocated( this%cross_sections_ ) ) then
      do ndx = 1,size( this%cross_sections_ )
        if( associated( this%cross_sections_( ndx )%val_ ) ) then
          deallocate( this%cross_sections_( ndx )%val_ )
        endif
      enddo
      deallocate( this%cross_sections_ )
    end if

    if( allocated( this%handles_ ) ) then
      deallocate( this%handles_ )
    end if

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_cross_section_warehouse
