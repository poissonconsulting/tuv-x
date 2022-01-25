! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The photolysis_component_set module

!> The component_set_t type and related functions
module photolysis_component_set

  use photolysis_component,            only : component_ptr
  use musica_constants,                only : musica_dk, musica_ik
  use musica_string,                   only : string_t

  implicit none
  private

  public :: component_set_t

  !> An ordered set of model components
  !!
  !! \todo add detailed description and example of component_set_t usage
  !!
  type :: component_set_t
    private
    !> Photolysis Components
    type(component_ptr), allocatable :: components_(:)
  contains
    !> Returns the index of a component in the set
    procedure :: index => component_set_index
    !> Returns the number of components in the set
    procedure :: size => component_set_size
    !> Adds a component to the set
    procedure :: add
    !> @name Returns a component
    !! @{
    procedure, private :: get_by_index
    procedure, private :: get_by_name
    generic :: get => get_by_index, get_by_name
  end type component_set_t

  !> Photolysis component set constructor
  interface component_set_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Creates a new component set
  function constructor( ) result( new_obj )

    !> New component set
    type(component_set_t), pointer :: new_obj

    allocate( new_obj )
    allocate( new_obj%components_( 0 ) )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the index of a component in the set
  integer(kind=musica_ik) function component_set_index( this, component )

    use musica_assert,                 only : assert, die_msg
    use photolysis_component,          only : component_t
    use musica_data_type,              only : data_type_t

    !> Photolysis component set
    class(component_set_t), intent(in) :: this
    !> Photolysis component to find
    class(component_t), intent(in) :: component

    integer(kind=musica_ik) :: i_comp

    component_set_index = 0
    call assert( 880659416, allocated( this%components_ ) )
    do i_comp = 1, size( this%components_ )
      component_set_index = component_set_index + 1
      if( this%components_( i_comp )%val_%name( ) .eq. component%name( ) ) then
        return
      end if
    end do
    call die_msg( 486178542, "Photolysis component set does not include component '"//     &
                  component%name( )//"'" )

  end function component_set_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the number of components in the set
  integer(kind=musica_ik) function component_set_size( this )

    use musica_assert,                 only : assert
    use musica_data_type,              only : data_type_t

    !> Photolysis component set
    class(component_set_t), intent(in) :: this

    call assert( 545922635, allocated( this%components_ ) )
    component_set_size = size( this%components_ )

  end function component_set_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Adds a component to the set
  subroutine add( this, new_component )

    use musica_assert,                 only : assert
    use photolysis_component,          only : component_t

    !> Photolysis component set
    class(component_set_t), intent(inout) :: this
    !> Photolysis component to add
    class(component_t), pointer, intent(in) :: new_component

    type(component_ptr) :: new_ptr

    call assert( 995196015, allocated( this%components_ ) )
    new_ptr%val_ => new_component
    this%components_ = [ this%components_, new_ptr ]
    new_ptr%val_ => null( )

  end subroutine add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Gets a component by index
  function get_by_index( this, index ) result( comp )

    use musica_assert,                 only : assert
    use photolysis_component,          only : component_t

    !> Photolysis component
    class(component_t), pointer :: comp
    !> Photolysis component set
    class(component_set_t), intent(in) :: this
    !> Index of the component in the set
    integer(kind=musica_ik) :: index

    call assert( 202250054, allocated( this%components_ ) )
    call assert( 314568399, index .ge. 1 .and.                                &
                            index .le. size( this%components_ ) )
    comp => this%components_( index )%val_

  end function get_by_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Gets a component by name
  function get_by_name( this, name, found ) result( comp )

    use musica_assert,                 only : assert, die_msg
    use photolysis_component,          only : component_t

    !> Photolysis component
    class(component_t), pointer :: comp
    !> Photolysis component set
    class(component_set_t), intent(in) :: this
    !> Photolysis component name
    character(len=*), intent(in) :: name
    !> Flag indicating whether the variable was found
    logical, intent(out), optional :: found

    integer(kind=musica_ik) :: i_comp

    if( present( found ) ) found = .false.
    comp => null( )
    call assert( 423527971, allocated( this%components_ ) )
    do i_comp = 1, size( this%components_ )
      if( this%components_( i_comp )%val_%name( ) .eq. name ) then
        if( present( found ) ) found = .true.
        comp => this%components_( i_comp )%val_
        return
      end if
      end do
    if( .not. present( found ) ) then
      call die_msg( 248106760, "Photolysis component '"//name//"' not found in "//        &
                    "component set." )
    end if

  end function get_by_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module photolysis_component_set
