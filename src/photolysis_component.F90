! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The photolysis component module

!> The abstract component_t type and related functions
module photolysis_component

  implicit none
  private

  public :: component_t, component_ptr

  !> Model component
  !!
  !! Model components calculate diagnostics for a given model state and/or
  !! advance the model state for a given timestep.
  !!
  !! \todo add full description and example usage for component_t
  !!
  type, abstract :: component_t
  contains
    !> Returns the name of the component
    procedure(component_name), deferred :: name
    !> Returns a description of the component purpose
    procedure(description), deferred :: description
    !> Update the component for a given timestep
    procedure(upDate), deferred :: upDate
  end type component_t

  !> Unique pointer for component_t objects
  type :: component_ptr
    class(component_t), pointer :: val_
  contains
    !> Finalize the pointer
    final :: component_ptr_finalize
  end type component_ptr

interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns the name of the component
  type(string_t) function component_name( this )
    use musica_string,                 only : string_t

    import component_t
    !> Model component
    class(component_t), intent(in) :: this
  end function component_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Returns a description of the component purpose
  type(string_t) function description( this )
    use musica_string,                 only : string_t

    import component_t
    !> Model component
    class(component_t), intent(in) :: this
  end function description

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update the component for a given timestep
  subroutine upDate( this, la_srb, SphericalGeom, GridWareHouse, ProfileWareHouse, radiationFld )
    use micm_grid_warehouse,           only : grid_warehouse_t
    use micm_Profile_warehouse,        only : Profile_warehouse_t
    use musica_constants,              only : dk => musica_dk
    use spherical_geom_type,           only : spherical_geom_t
    use la_srb_type,                   only : la_srb_t
    use abstract_radXfer_type,         only : radField_t

    import component_t
    !> Model component
    class(component_t), intent(inout)        :: this
    type(grid_warehouse_t), intent(inout)    :: GridWareHouse
    type(Profile_warehouse_t), intent(inout) :: ProfileWareHouse
    type(spherical_geom_t), intent(inout)    :: SphericalGeom
    type(la_srb_t), intent(inout)            :: la_srb
    class(radField_t), allocatable           :: radiationFld
  end subroutine upDate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the unique pointer
  elemental subroutine component_ptr_finalize( this )

    !> Component pointer
    type(component_ptr), intent(inout) :: this

    if( associated( this%val_ ) ) then
      deallocate( this%val_ )
      this%val_ => null( )
    end if

  end subroutine component_ptr_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module photolysis_component
