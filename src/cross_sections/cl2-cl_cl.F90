! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> This cl2_cl_cl_cross_section module

!> The cl2+hv->cl_cl cross_section type and related functions
module tuvx_cross_section_cl2_cl_cl

  use tuvx_cross_section, only : cross_section_t, base_constructor

  implicit none

  private
  public :: cross_section_cl2_cl_cl_t

  !> Calculator for base_cross_section
  type, extends(cross_section_t) :: cross_section_cl2_cl_cl_t
  contains
    !> Calculate the cross section
    procedure :: calculate => run
  end type cross_section_cl2_cl_cl_t

  !> Constructor
  interface cross_section_cl2_cl_cl_t
    module procedure constructor
  end interface cross_section_cl2_cl_cl_t

contains

  !> Initialize the cross section
  function constructor( config, grid_warehouse, profile_warehouse, at_mid_point ) result ( this )

    use musica_config,    only : config_t
    use musica_constants, only : lk => musica_lk
    use tuvx_grid_warehouse,    only : grid_warehouse_t
    use tuvx_profile_warehouse, only : profile_warehouse_t


    !> Cross section calculator
    logical(lk), optional, intent(in)          :: at_mid_point
    class(cross_section_t), pointer  :: this
    type(config_t), intent(inout)              :: config
    type(grid_warehouse_t), intent(inout)      :: grid_warehouse
    type(profile_warehouse_t), intent(inout)   :: profile_warehouse

    allocate ( cross_section_cl2_cl_cl_t :: this )
    call base_constructor( this, config, grid_warehouse, profile_warehouse )
  end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate cross section for a given set of environmental conditions
  function run( this, grid_warehouse, profile_warehouse, at_mid_point ) result( cross_section )

    use musica_constants,           only : dk => musica_dk, ik => musica_ik, lk => musica_lk
    use tuvx_grid_warehouse,        only : grid_warehouse_t
    use tuvx_grid,               only : abs_1d_grid_t
    use tuvx_profile_warehouse,     only : profile_warehouse_t
    use tuvx_profile,               only : abs_profile_t
    use musica_string,              only : string_t

    !> Arguments
    class(cross_section_cl2_cl_cl_t), intent(in) :: this
    logical(lk), intent(in), optional            :: at_mid_point
    !> The warehouses
    type(grid_warehouse_t), intent(inout)        :: grid_warehouse
    type(profile_warehouse_t), intent(inout)     :: profile_warehouse
    !> Calculated cross section
    real(kind=dk), allocatable                   :: cross_section(:,:)

    !> Local variables
    character(len=*), parameter :: Iam  = 'cl2+hv->cl+cl cross section calculate: '
    integer(ik), parameter  :: iONE = 1_ik
    real(dk), parameter     :: rONE = 1.0_dk
    integer(ik) :: lambdaNdx, vertNdx, nzdim
    real(dk)    :: aa, bb, bbsq, alpha, ex1, ex2
    real(dk), allocatable   :: modelTemp(:)
    class(abs_1d_grid_t), pointer :: lambdaGrid
    class(abs_1d_grid_t), pointer :: zGrid
    class(abs_profile_t), pointer :: Temperature
    type(string_t) :: Handle

    write(*,*) Iam,'entering'

    Handle = 'Photolysis, wavelength' ; lambdaGrid => grid_warehouse%get_grid( Handle )
    Handle = 'Vertical Z'             ; zGrid => grid_warehouse%get_grid( Handle )
    Handle = 'Temperature'            ; Temperature => profile_warehouse%get_Profile( Handle )

    nzdim = zGrid%ncells_ + iONE
    if( present(at_mid_point) ) then
      if( at_mid_point ) then
        nzdim = nzdim - iONE
        modelTemp = Temperature%mid_val_
      else
        modelTemp = Temperature%edge_val_
      endif
    else
      modelTemp = Temperature%edge_val_
    endif

    allocate( cross_section(lambdaGrid%ncells_,zGrid%ncells_+iONE) )

    associate( wc => lambdaGrid%mid_ )
    do vertNdx = iONE,nzdim
      aa    = 402.7_dk/modelTemp(vertNdx)
      bb    = exp( aa )
      bbsq  = bb * bb
      alpha = (bbsq - rONE)/(bbsq + rONE)
      do lambdaNdx = iONE,lambdaGrid%ncells_
        ex1 = 27.3_dk  * exp(-99.0_dk * alpha * (log(329.5_dk/wc(lambdaNdx)))**2)
        ex2 =  .932_dk * exp(-91.5_dk * alpha * (log(406.5_dk/wc(lambdaNdx)))**2)
        cross_section(lambdaNdx,vertNdx) = 1.e-20_dk * sqrt(alpha) * (ex1 + ex2)
      enddo
    enddo
    end associate

    cross_section = transpose( cross_section )

    write(*,*) Iam,'exiting'

  end function run

end module tuvx_cross_section_cl2_cl_cl
