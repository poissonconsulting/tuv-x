! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The micm_radXfer_xsect_warehouse module

!> The radXfer_xsect_warehouse_t type and related functions
!!
module micm_radXfer_xsect_warehouse

  use micm_abs_cross_section_type,     only : abs_cross_section_ptr
  use musica_constants,                only : musica_dk, musica_ik
  use musica_string,                   only : string_t

  implicit none

  private
  public :: radXfer_xsect_warehouse_t

  !> Radiative xfer cross section type
  type :: radXfer_xsect_warehouse_t
    !> cross section calculators
    type(abs_cross_section_ptr), allocatable :: cross_section_objs_(:)
    !> Current cross section values
    real(kind=musica_dk), allocatable :: cross_section_values_(:,:)
    type(string_t), allocatable       :: reaction_key(:)
  contains
    !> Update the object for new environmental conditions
    procedure :: update_for_new_environmental_state
    !> Finalize the object
    final :: finalize
  end type radXfer_xsect_warehouse_t

  !> radXfer_xsect_warehouse_t constructor
  interface radXfer_xsect_warehouse_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of radXfer_xsect_warehouse_t objects
  function constructor( config,mdlLambdaEdge ) result( radXfer_xsect_obj )

    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use musica_constants,              only : musica_rk
    use musica_assert,                 only : die_msg
    use micm_cross_section_factory,    only : cross_section_builder

    real(musica_dk), intent(in)      :: mdlLambdaEdge(:)
    !> Kinetics configuration data
    type(config_t), intent(inout) :: config
    !> New radiative xfer cross section obj
    class(radXfer_xsect_warehouse_t), pointer :: radXfer_xsect_obj

    !> local variables
    integer :: nSize, ndx
    character(len=*), parameter :: Iam = "Radiative xfer xsect constructor: "
    type(config_t) :: reaction_set, reaction_config
    type(config_t) :: cross_section_config
    class(iterator_t), pointer :: iter
    type(abs_cross_section_ptr) :: cross_section_ptr
    character(:), allocatable   :: jsonkey
    character(len=32)           :: keychar
    type(string_t)              :: areaction_key
    type(string_t), allocatable :: netcdfFiles(:)

    allocate( radXfer_xsect_obj )

    associate(new_obj=>radXfer_xsect_obj)

    allocate( string_t :: new_obj%reaction_key(0) )

    allocate( new_obj%cross_section_objs_(0) )

    jsonkey = 'radiative xfer cross sections'
    call config%get( jsonkey, reaction_set, Iam )
    iter => reaction_set%get_iterator( )
!-----------------------------------------------------------------------------
!> iterate over cross sections
!-----------------------------------------------------------------------------
    do while( iter%next( ) )
      keychar = reaction_set%key(iter)
      areaction_key = keychar 
      write(*,*) ' '
      write(*,*) Iam,'key = ',trim(keychar)
      new_obj%reaction_key = [new_obj%reaction_key,areaction_key]
      call reaction_set%get( iter, reaction_config, Iam )
!-----------------------------------------------------------------------------
!> cross section first
!-----------------------------------------------------------------------------
      call reaction_config%get( "cross section", cross_section_config, Iam )
      cross_section_ptr%val_ => cross_section_builder( cross_section_config,mdlLambdaEdge )
      new_obj%cross_section_objs_ = [new_obj%cross_section_objs_,cross_section_ptr]
    end do

    deallocate( iter )

    nSize = size(new_obj%cross_section_objs_)
    write(*,*) ' '
    write(*,'(a,''There are '',i3,'' cross sections'')') Iam,nSize

!-----------------------------------------------------------------------------
!> setup cross section arrays
!-----------------------------------------------------------------------------
    allocate( new_obj%cross_section_values_(0,0) )

    end associate

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update the object for new environmental conditions
  subroutine update_for_new_environmental_state( this, environment, nwave )

    use micm_environment, only : environment_t
    use musica_assert,    only : die_msg

    !> Kinetics
    class(radXfer_xsect_warehouse_t), intent(inout) :: this
    integer(musica_ik), intent(in)         :: nwave
    !> Environmental conditions
    class(environment_t), intent(in) :: environment

    character(len=*), parameter :: Iam = 'update_for_new_environmental_state: '
    integer(kind=musica_ik) :: ndx
    real(musica_dk), allocatable :: a_cross_section(:)
    real(musica_dk), allocatable :: cross_section_tray(:)

    write(*,*) ' '
    write(*,*) Iam,'entering'

    allocate(cross_section_tray(0))
    do ndx = 1, size(this%cross_section_objs_)
      associate( calc_ftn => this%cross_section_objs_(ndx)%val_ )
        a_cross_section = calc_ftn%calculate( environment )
      end associate
      cross_section_tray = [cross_section_tray,a_cross_section]
    end do

    this%cross_section_values_ = reshape( cross_section_tray, &
                                          (/nwave,size(this%cross_section_objs_) /) )

    write(*,*) Iam,'size of cross section values = ',&
        size(this%cross_section_values_,dim=1), size(this%cross_section_values_,dim=2)

    write(*,*) Iam,'exiting'

  end subroutine update_for_new_environmental_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the radXfer xsect warehouse
  subroutine finalize( this )

    !> radXfer xsect warehouse
    type(radXfer_xsect_warehouse_t), intent(inout) :: this

    integer(kind=musica_ik) :: ndx
    character(len=*), parameter :: Iam = 'radXfer_xsect finalize: '

    write(*,*) Iam,'entering'

    if( allocated( this%cross_section_values_ ) ) then
      deallocate( this%cross_section_values_ )
    endif
    if( allocated( this%cross_section_objs_ ) ) then
      do ndx = 1,size(this%cross_section_objs_)
        if( associated( this%cross_section_objs_(ndx)%val_ ) ) then
          deallocate( this%cross_section_objs_(ndx)%val_ )
        endif
      enddo
      deallocate( this%cross_section_objs_ )
    end if

    if( allocated( this%reaction_key ) ) then
      deallocate( this%reaction_key )
    end if

    write(*,*) Iam,'exiting'

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_radXfer_xsect_warehouse
