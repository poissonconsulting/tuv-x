! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_spectral_weight
  ! The spectral_weight type and related functions

  use musica_constants,                only : dk => musica_dk
  use musica_config,                   only : config_t
  use tuvx_grid_warehouse,             only : grid_warehouse_t
  use tuvx_profile_warehouse,          only : profile_warehouse_t

  implicit none

  private
  public :: spectral_weight_t, spectral_weight_parms_t, spectral_weight_ptr,  &
            base_constructor

  type spectral_weight_parms_t
    real(dk), allocatable :: temperature(:)
    real(dk), allocatable :: array(:,:)
  end type spectral_weight_parms_t

  type :: spectral_weight_t
    ! Calculator for spectral weights
    type(spectral_weight_parms_t), allocatable :: spectral_weight_parms(:)
  contains
    procedure :: calculate => run
    procedure, private :: add_points
    final     :: finalize
  end type spectral_weight_t

  !> Pointer type for building sets of dose rate constants
  type :: spectral_weight_ptr
    class(spectral_weight_t), pointer :: val_ => null( )
  end type spectral_weight_ptr

  interface spectral_weight_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result( new_spectral_weight )
    ! Instantiate the base spectral weight type

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : string_t

    class(spectral_weight_t),  pointer       :: new_spectral_weight  ! New spectral weight object
    type(config_t),            intent(inout) :: config               ! Spectral weight configuration data
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse       ! Grid warehouse
    type(profile_warehouse_t), intent(inout) :: profile_warehouse    ! Profile warehouse

    ! Local variables
    type(string_t) :: required_keys(1), optional_keys(4)

    required_keys(1) = 'type'
    optional_keys(1) = 'netcdf files'
    optional_keys(2) = 'lower extrapolation'
    optional_keys(3) = 'upper extrapolation'
    optional_keys(4) = 'name'
    call assert_msg( 124969901,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "base spectral wght." )

    allocate( new_spectral_weight )
    call base_constructor( new_spectral_weight, config, grid_warehouse,       &
                           profile_warehouse )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine base_constructor( this, config, grid_warehouse,                  &
      profile_warehouse )

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_netcdf,                   only : netcdf_t
    use tuvx_profile,                  only : profile_t
    use tuvx_util,                     only : inter2

    class(spectral_weight_t),  pointer       :: this              ! New spectral weight object
    type(config_t),            intent(inout) :: config            ! Spectral weight configuration data
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse    ! Grid warehouse
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! Profile warehouse

    ! Local variables
    real(dk), parameter    :: rZERO = 0.0_dk
    real(dk), parameter    :: rONE  = 1.0_dk
    character(len=*), parameter    :: Iam = 'base spectral weight initialize: '
    character(len=*), parameter    :: Hdr = 'spectral_weight_'

    integer :: retcode
    integer :: parmNdx, fileNdx
    integer :: nParms
    real(dk), allocatable :: data_lambda(:)
    real(dk), allocatable :: data_parameter(:)
    logical :: found
    type(netcdf_t), allocatable :: netcdf_obj
    type(string_t), allocatable :: netcdfFiles(:)
    class(grid_t), pointer      :: lambdaGrid => null()

    ! Get model wavelength grid
    lambdaGrid => grid_warehouse%get_grid( "wavelength", "nm" )

    ! Get spectral wght netcdf filespec
    call config%get( 'netcdf files', netcdfFiles, Iam, found = found )

    has_netcdf_file: if( found ) then
      allocate( this%spectral_weight_parms( size( netcdfFiles ) ) )
      file_loop: do fileNdx = 1, size( this%spectral_weight_parms )
        allocate( netcdf_obj )
        ! Read netcdf spectral wght parameters
        call netcdf_obj%read_netcdf_file(                                     &
                             file_path = netcdfFiles( fileNdx )%to_char( ),   &
                             variable_name = Hdr )
        nParms = size( netcdf_obj%parameters, dim = 2 )
        call assert_msg( 222975865, nParms > 0, Iam//'File: ' //              &
                         trim(netcdfFiles( fileNdx )%to_char( ) ) //          &
                         '  parameters array has < 1 column' )

        ! Interpolate from data to model wavelength grid
        if( allocated( netcdf_obj%wavelength ) ) then
          if( .not. allocated( this%spectral_weight_parms( fileNdx )%array ) )&
              then
            allocate( this%spectral_weight_parms( fileNdx                     &
                                       )%array( lambdaGrid%ncells_, nParms ) )
          endif
          do parmNdx = 1, nParms
            data_lambda    = netcdf_obj%wavelength
            data_parameter = netcdf_obj%parameters( :, parmNdx )
            call this%add_points( config, data_lambda, data_parameter )
            call inter2(                                                      &
              xto = lambdaGrid%edge_,                                         &
              yto = this%spectral_weight_parms( fileNdx )%array( :, parmNdx ),&
              xfrom = data_lambda,                                            &
              yfrom = data_parameter, ierr = retcode )
          enddo
        else
          this%spectral_weight_parms( fileNdx )%array = netcdf_obj%parameters
        endif
        if( allocated( netcdf_obj%temperature ) ) then
          this%spectral_weight_parms( fileNdx )%temperature =                 &
              netcdf_obj%temperature
        endif
        deallocate( netcdf_obj )
      enddo file_loop
    endif has_netcdf_file

    if( associated( lambdaGrid ) ) deallocate( lambdaGrid )

  end subroutine base_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function run( this, grid_warehouse, profile_warehouse )                     &
      result( spectral_weight )
    ! Calculate the spectral wght for a given set of environmental conditions

    use musica_string,                   only : string_t
    use tuvx_grid,                       only : grid_t

    class(spectral_weight_t),  intent(in)     :: this
    type(grid_warehouse_t),    intent(inout)  :: grid_warehouse
    type(profile_warehouse_t), intent(inout) :: profile_warehouse
    real(kind=dk), allocatable               :: spectral_weight(:)

    ! Local variables
    character(len=*), parameter :: Iam = 'spectral weight calculate: '

    spectral_weight = this%spectral_weight_parms(1)%array( :, 1 )

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine add_points( this, config, data_lambda, data_parameter )
    ! Adds points to the spectral wght grid based on configuration data

    use musica_string,                 only : string_t
    use musica_assert,                 only : assert_msg, die_msg
    use tuvx_util,                     only : addpnt

    class(spectral_weight_t), intent(in)    :: this
    type(config_t),           intent(inout) :: config
    real(dk), allocatable,    intent(inout) :: data_lambda(:)
    real(dk), allocatable,    intent(inout) :: data_parameter(:)

    ! Local variables
    real(dk), parameter :: rZERO = 0.0_dk
    real(dk), parameter :: rONE  = 1.0_dk
    real(dk), parameter :: deltax = 1.e-5_dk
    character(len=*), parameter :: Iam = 'spectral_weight; add_points'

    integer  :: nRows
    real(dk) :: lowerLambda, upperLambda
    real(dk) :: addpnt_val_lower, addpnt_val_upper
    type(string_t)  :: addpnt_type
    type(config_t)  :: extrap_config
    logical         :: found
    character(len=:), allocatable :: number
    type(string_t) :: required_keys(1), optional_keys(1)

    required_keys(1) = "type"
    optional_keys(1) = "value"

    ! add endpoints to data arrays; first the lower bound
    nRows = size( data_lambda )
    lowerLambda = data_lambda(1) ; upperLambda = data_lambda( nRows )

    addpnt_val_lower = rZERO
    call config%get( 'lower extrapolation', extrap_config, Iam, found = found )
    if( found ) then
      call assert_msg( 671608111,                                             &
                       extrap_config%validate( required_keys, optional_keys ),&
                       "Bad format for extrapolation in base spectral wght." )
      call extrap_config%get( "type", addpnt_type, Iam )
      if( addpnt_type == "boundary" ) then
        addpnt_val_lower = data_parameter(1)
      elseif( addpnt_type == "constant" ) then
        call extrap_config%get( "value", addpnt_val_lower, Iam )
      else
        call die_msg( 316405972,                                              &
                      "Bad extrapolation type: '"//addpnt_type//"'" )
      endif
    endif

    ! add endpoints to data arrays; now the upper bound
    addpnt_val_upper = rZERO
    call config%get( 'upper extrapolation', extrap_config, Iam, found = found )
    if( found ) then
      call assert_msg( 918590096,                                             &
                       extrap_config%validate( required_keys, optional_keys ),&
                       "Bad format for extrapolation in base spectral wght." )
      call extrap_config%get( "type", addpnt_type, Iam )
      if( addpnt_type == "boundary" ) then
        addpnt_val_upper = data_parameter( nRows )
      elseif( addpnt_type == "constant" ) then
        call extrap_config%get( "value", addpnt_val_upper, Iam )
      else
        call die_msg( 302970880,                                              &
                      "Bad extrapolation type: '"//addpnt_type//"'" )
      endif
    endif

    call addpnt( x = data_lambda, y = data_parameter,                         &
                 xnew = ( rONE - deltax ) * lowerLambda,                      &
                 ynew = addpnt_val_lower )
    call addpnt( x = data_lambda, y = data_parameter, xnew = rZERO,           &
                 ynew = addpnt_val_lower )
    call addpnt( x = data_lambda, y = data_parameter,                         &
                 xnew = ( rONE + deltax ) * upperLambda,                      &
                 ynew = addpnt_val_upper )
    call addpnt( x = data_lambda, y = data_parameter, xnew = 1.e38_dk,        &
                 ynew = addpnt_val_upper )

  end subroutine add_points

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Finalizes the spectral wght type

    type(spectral_weight_t), intent(inout) :: this

    ! Local variables
    integer :: ndx

    if( allocated( this%spectral_weight_parms ) ) then
      do ndx = 1,size( this%spectral_weight_parms )
        if( allocated( this%spectral_weight_parms( ndx )%array ) ) then
          deallocate( this%spectral_weight_parms( ndx )%array )
        endif
        if( allocated( this%spectral_weight_parms( ndx )%temperature ) ) then
          deallocate( this%spectral_weight_parms( ndx )%temperature )
        endif
      enddo
      deallocate( this%spectral_weight_parms )
    endif

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_spectral_weight
