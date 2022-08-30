! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_radiator_warehouse
! A class that holds and gives out
! :f:type:`~tuvx_radiator/radiator_t`'s built by the
! :f:mod:`~tuvx_radiator_factory`.


  use musica_constants,                only : musica_dk
  use musica_iterator,                 only : iterator_t
  use musica_string,                   only : string_t
  use tuvx_radiator,                   only : radiator_ptr

  implicit none

  private
  public :: radiator_warehouse_t, warehouse_iterator_t

  type radiator_warehouse_t
    ! Radiator warehouse

    private
    type(radiator_ptr), allocatable :: radiators_(:) ! Radiators
    type(string_t),     allocatable :: handle_(:) ! Radiator "handles"
  contains
    !> @name Returns a pointer to a requested radiator
    !! @{
    procedure, private :: get_radiator_from_handle
    procedure, private :: get_radiator_from_iterator
    generic   :: get_radiator => get_radiator_from_handle, get_radiator_from_iterator
    !> Returns the index associated with a given radiator name
    procedure, public  :: get_radiator_ndx_from_handle
    !> @}
    !> Returns whether a radiator exists in the warehouse
    procedure :: in_warehouse
    !> Gets an iterator for the warehouse
    procedure :: get_iterator
    !> Accumulates the state of all radiators in the warehouse
    procedure :: accumulate_states
    !> Cleans up memory
    final     :: finalize
  end type radiator_warehouse_t

  type warehouse_iterator_t
    !  Radiator warehouse iterator

    type(radiator_warehouse_t), pointer :: warehouse_ => null( ) ! Pointer to the radiator warehouse
    integer :: id_ = 0 ! Current index in the data set
  contains
    !> Advances to the next key-value pair
    procedure :: next => iterator_next
    !> Resets the iterator
    procedure :: reset => iterator_reset
  end type warehouse_iterator_t

  interface radiator_warehouse_t
    ! radiator_warehouse_t constructor
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse ) result( radiator_warehouse )
    ! Constructs radiator_warehouse_t abjects

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use musica_string,                 only : string_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_radiator_factory,         only : radiator_builder

    type(config_t), intent(inout)        :: config ! Radiator configuration
    type(grid_warehouse_t), intent(inout):: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    class(radiator_warehouse_t), pointer :: radiator_warehouse ! A :f:type:`~tuvx_radiator_warehouse/radiator_warehouse_t`

    ! local variables
    character(len=*), parameter :: Iam = "Radiator warehouse constructor: "
    integer                     :: ndx
    type(config_t)              :: radiator_config
    class(iterator_t), pointer  :: iter
    type(radiator_ptr)          :: aRadiator

    allocate( radiator_warehouse )
    allocate( radiator_warehouse%radiators_(0) )
    allocate( radiator_warehouse%handle_(0) )

    iter => config%get_iterator()
    do while( iter%next() )
      call config%get( iter, radiator_config, Iam )

      ! build and store the radiator
      aRadiator%val_ => radiator_builder( radiator_config, grid_warehouse )
      radiator_warehouse%radiators_ =                                         &
          [ radiator_warehouse%radiators_, aRadiator ]
      radiator_warehouse%handle_ =                                            &
          [ radiator_warehouse%handle_, aRadiator%val_%handle_ ]
    end do
    deallocate( iter )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_radiator_from_handle( this, radiator_handle )                  &
      result( radiator )
    ! Returns a pointer to a requested radiator

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : string_t
    use tuvx_radiator,                 only : radiator_t

    class(radiator_warehouse_t), intent(inout) :: this ! This :f:type:`~tuvx_radiator_warehouse/radiator_warehouse_t`
    type(string_t),              intent(in)    :: radiator_handle ! Name associated with requested radiator
    class(radiator_t),           pointer       :: radiator ! Pointer to the requested radiator of type :f:type:`~tuvx_radiator/radiator_t`

    ! Local variables
    character(len=*), parameter :: Iam = 'radiator warehouse get radiator'
    integer :: ndx
    logical :: found

    found = .false.
    do ndx = 1, size( this%handle_ )
      if( radiator_handle .eq. this%handle_( ndx ) ) then
        found = .true.
        exit
      endif
    end do
    call assert_msg( 200578546, found,                                        &
                     "Invalid radiator handle: '"//                           &
                     radiator_handle%to_char()//"'" )
    radiator => this%radiators_(ndx)%val_

  end function get_radiator_from_handle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_radiator_ndx_from_handle( this, radiator_handle )              &
      result( index )
    ! Returns the index associated with a given radiator
    !
    ! If the radiator is not found, returns -1

    use tuvx_radiator,                 only : radiator_t
    use musica_string,                 only : string_t

    class(radiator_warehouse_t), intent(inout) :: this ! This :f:type:`~tuvx_radiator_warehouse/radiator_warehouse_t`
    type(string_t),              intent(in)    :: radiator_handle ! Requested radiator name
    integer                                    :: index ! Index of requested radiator in warehouse

    ! Local variables
    character(len=*), parameter :: Iam =                                      &
        'radiator warehouse get radiator index'
    logical :: found

    found = .false.
    do index = 1, size( this%handle_ )
      if( radiator_handle .eq. this%handle_( index ) ) then
        found = .true.
        exit
      endif
    end do
    if( .not. found ) index = -1

  end function get_radiator_ndx_from_handle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_radiator_from_iterator( this, iterator ) result( radiator )
    ! Returns a pointer to a radiator in the warehouse from an iterator

    use tuvx_radiator,                 only : radiator_t

    class(radiator_warehouse_t), intent(inout) :: this ! This :f:type:`~tuvx_radiator_warehouse/radiator_warehouse_t`
    type(warehouse_iterator_t),  intent(in)    :: iterator ! ! A :f:type:`~tuvx_radiator_warehouse/warehouse_iterator_t`
    class(radiator_t),           pointer       :: radiator ! Pointer to requested :f:type:`~tuvx_radiator/radiator_t`

    ! Local variables
    character(len=*), parameter :: Iam =                                      &
        'radiator warehouse get radiator from iterator'

    radiator => this%radiators_( iterator%id_ )%val_

  end function get_radiator_from_iterator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function in_warehouse( this, radiator_handle )
    ! Returns whether a radiator exists in the warehouse

    use musica_string,      only : string_t
    use musica_assert,      only : die_msg

    class(radiator_warehouse_t), intent(inout) :: this ! This :f:type:`~tuvx_radiator_warehouse/radiator_warehouse_t`
    type(string_t),              intent(in)    :: radiator_handle ! Name associated with requested radiator
    logical                                    :: in_warehouse ! Flag indicating whether the radiator exists in the warehouse

    ! Local variables
    character(len=*), parameter :: Iam = 'radiator in warehouse'
    integer :: ndx

    in_warehouse = .false.
    do ndx = 1, size( this%handle_ )
      if( radiator_handle == this%handle_( ndx ) ) then
        in_warehouse = .true.
        exit
      endif
    end do

  end function in_warehouse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function get_iterator( this )
    ! Gets an iterator for the radiator warehouse

    use musica_assert,                 only : assert

    class(warehouse_iterator_t), pointer            :: get_iterator ! Pointer to the :f:type:`~tuvx_radiator_warehouse/warehouse_iterator_t`
    class(radiator_warehouse_t), intent(in), target :: this ! This :f:type:`~tuvx_radiator_warehouse/radiator_warehouse_t`

    call assert( 753334333, allocated( this%radiators_ ) )
    allocate( warehouse_iterator_t :: get_iterator )
    select type( iter => get_iterator )
      type is( warehouse_iterator_t )
        iter%warehouse_ => this
    end select

  end function get_iterator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine accumulate_states( this, state )
    ! Accumulates the states off all radiators in the warehouse into a
    ! single representative state.

    use tuvx_radiator,                 only : radiator_state_t

    class(radiator_warehouse_t), intent(in)    :: this
    class(radiator_state_t),     intent(inout) :: state

    call state%accumulate( this%radiators_ )

  end subroutine accumulate_states

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function iterator_next( this ) result( continue )
    ! Advances the iterator
    !
    ! Returns false if the end of the collection has been reached

    class(warehouse_iterator_t), intent(inout) :: this ! A :f:type:`~tuvx_radiator_warehouse/warehouse_iterator_t`

    logical :: continue

    this%id_ = this%id_ + 1
    continue = this%id_ <= size( this%warehouse_%radiators_ )

  end function iterator_next

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine iterator_reset( this )
    ! Resets the iterator

    class(warehouse_iterator_t), intent(inout) :: this ! This :f:type:`~tuvx_radiator_warehouse/warehouse_iterator_t`

    this%id_ = 0

  end subroutine iterator_reset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Finalize the radiator warehouse

    type(radiator_warehouse_t), intent(inout) :: this ! This :f:type:`~tuvx_radiator_warehouse/radiator_warehouse_t`

    integer :: ndx
    character(len=*), parameter :: Iam = 'radiator_warehouse finalize: '

    if( allocated( this%radiators_ ) ) then
      do ndx = 1, size( this%radiators_ )
        if( associated( this%radiators_( ndx )%val_ ) ) then
          deallocate( this%radiators_( ndx )%val_ )
        end if
      end do
      deallocate( this%radiators_ )
    endif

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_radiator_warehouse
