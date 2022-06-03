! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> This oclo_cross_section module

!> The oclo_cross_section type and related functions
module tuvx_cross_section_oclo

  use tuvx_cross_section, only : cross_section_t, base_constructor

  implicit none

  private
  public ::cross_section_oclo_t

  !> Calculator for oclo_cross_section
  type, extends(cross_section_t) ::cross_section_oclo_t
    !> The cross section array
  contains
    !> Calculate the cross section
    procedure :: calculate => run
  end type cross_section_oclo_t

  !> Constructor
  interface cross_section_oclo_t
    module procedure constructor
  end interface cross_section_oclo_t

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

    allocate (cross_section_oclo_t :: this )
    call base_constructor( this, config, grid_warehouse, profile_warehouse )
  end function constructor


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
    class(cross_section_oclo_t), intent(in)  :: this
    logical(lk), optional, intent(in)        :: at_mid_point
    !> The warehouses
    type(grid_warehouse_t), intent(inout)    :: grid_warehouse
    type(profile_warehouse_t), intent(inout) :: profile_warehouse
    !> Calculated cross section
    real(kind=dk), allocatable               :: cross_section(:,:)

    !> Local variables
    character(len=*), parameter :: Iam = 'oclo cross section calculate: '
    integer(ik), parameter :: iONE  = 1_ik
    real(dk), parameter    :: rZERO = 0.0_dk
    integer(ik) :: ndx, nParms
    integer(ik) :: lambdaNdx, vertNdx, nzdim
    real(dk)    :: Tfac
    real(dk), allocatable  :: wrkCrossSection(:)
    real(dk), allocatable  :: modelTemp(:)
    class(abs_1d_grid_t), pointer :: zGrid
    class(abs_1d_grid_t), pointer :: lambdaGrid
    class(abs_profile_t), pointer :: mdlTemperature
    type(string_t) :: Handle

    write(*,*) Iam,'entering'

    Handle = 'Vertical Z'             ; zGrid => grid_warehouse%get_grid( Handle )
    Handle = 'Photolysis, wavelength' ; lambdaGrid => grid_warehouse%get_grid( Handle )
    Handle = 'Temperature'            ; mdlTemperature => profile_warehouse%get_Profile( Handle )

    nzdim = zGrid%ncells_ + iONE
    if( present(at_mid_point) ) then
      if( at_mid_point ) then
        nzdim = nzdim - iONE
        modelTemp = mdlTemperature%mid_val_
      else
        modelTemp = mdlTemperature%edge_val_
      endif
    else
      modelTemp = mdlTemperature%edge_val_
    endif

    allocate( cross_section(lambdaGrid%ncells_,nzdim) )
    allocate( wrkCrossSection(lambdaGrid%ncells_) )
    cross_section = rZERO

    associate( Temp => modelTemp, Xsection => this%cross_section_parms )
    nParms = size(Xsection)
    do vertNdx = iONE,nzdim
      if( Temp(vertNdx) <= Xsection(1)%temperature(1) ) then
        wrkCrossSection = Xsection(1)%array(:,1)
      elseif( Temp(vertNdx) >= Xsection(nParms)%temperature(1) ) then
        wrkCrossSection = Xsection(nParms)%array(:,1)
      else
        do ndx = 2,nParms
          if( Xsection(ndx)%temperature(1) > Temp(vertNdx) ) then
            exit
          endif
        enddo
        ndx = ndx - 1
        Tfac = (Temp(vertNdx) - Xsection(ndx)%temperature(1)) &
               /(Xsection(ndx+1)%temperature(1) - Xsection(ndx)%temperature(1))
        wrkCrossSection = Xsection(ndx)%array(:,1) &
                        + Tfac*(Xsection(ndx+1)%array(:,1) - Xsection(ndx)%array(:,1))
      endif
      cross_section(:,vertNdx) = wrkCrossSection
    enddo
    end associate

    cross_section = transpose( cross_section )

    write(*,*) Iam,'exiting'

  end function run

end module tuvx_cross_section_oclo
