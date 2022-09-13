! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_quantum_yield
  ! This base quantum yield module
  ! The base quantum yield type and related functions

  use musica_constants,                only : dk => musica_dk

  implicit none

  private
  public :: quantum_yield_t, quantum_yield_ptr, base_constructor

  type quantum_yield_parms_t
    real(dk), allocatable :: temperature(:) ! temperature in Kelvin
    real(dk), allocatable :: array(:,:) ! Parameters for calculating quantum yields (wavelength, parameter)
  contains
    ! returns the number of bytes needed to pack the parameters onto a buffer
    procedure :: pack_size => parms_pack_size
    ! packs the parameters onto a character buffer
    procedure :: mpi_pack => parms_mpi_pack
    ! unpacks parameters from a charcater buffer
    procedure :: mpi_unpack => parms_mpi_unpack
  end type quantum_yield_parms_t

  type quantum_yield_t
    ! Calculator for base quantum yield
    type(quantum_yield_parms_t), allocatable :: quantum_yield_parms(:)
  contains
    procedure :: calculate => run
    procedure :: add_points
    ! returns the number of bytes required to pack the object onto a buffer
    procedure :: pack_size
    ! packs the object onto a character buffer
    procedure :: mpi_pack
    ! unpacks an object from a character buffer
    procedure :: mpi_unpack
  end type quantum_yield_t

  type quantum_yield_ptr
    ! Pointer type for building sets of quantum yields
    class(quantum_yield_t), pointer :: val_ => null( )
  end type quantum_yield_ptr

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine base_constructor( this, config, grid_warehouse,                  &
      profile_warehouse )
    ! Performs default initialization for
    ! :f:type:`~tuvx_quantum_yield/quantum_yield_t` object
    ! Should only be called by sub-class constructors. Sub-classes can decide
    ! whether to call this function during construction to load standard
    ! NetCDF files and configuration options.
    !
    ! Reads NetCDF files specified in configuration array 'netcdf files'.
    ! Data from each NetCDF file will be loaded into an element of the
    ! quantum_yield_parms data member. If a NetCDF variable named
    ! quantum_yield_parameters is present, it will be used to populate
    ! the array data member of the :f:type:`~tuvx_quantum_yield/quantum_yield_parms_t` object.
    ! If a NetCDF variable named temperature is present, it will be
    ! used to populate the temperature data member of the
    ! :f:type:`~tuvx_quantum_yield/quantum_yield_parms_t` object.

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_interpolate,              only : interpolator_conserving_t
    use tuvx_netcdf,                   only : netcdf_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(quantum_yield_t),    intent(inout) :: this ! This :f:type:`~tuvx_quantum_yield/quantum_yield_t` calculator
    type(config_t),            intent(inout) :: config ! Quantum yield configuration data
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`

    ! Local variables
    character(len=*), parameter :: Iam = 'base quantum yield constructor'
    character(len=*), parameter :: Hdr = 'quantum_yield_'
    real(dk), parameter    :: rZERO = 0.0_dk
    real(dk), parameter    :: rONE  = 1.0_dk

    integer :: parmNdx, fileNdx
    integer :: nParms
    real(dk)    :: quantum_yield_constant
    real(dk), allocatable :: data_lambda(:)
    real(dk), allocatable :: data_parameter(:)
    logical :: found
    character(len=:), allocatable :: msg
    type(netcdf_t), allocatable   :: netcdf_obj
    type(string_t), allocatable   :: netcdfFiles(:)
    class(grid_t),  pointer       :: lambdaGrid => null( )
    type(interpolator_conserving_t) :: interpolator

    ! Get model wavelength grid
    lambdaGrid => grid_warehouse%get_grid( "wavelength", "nm" )

    ! get quantum yield netcdf filespec
    call config%get( 'netcdf files', netcdfFiles, Iam, found = found )
has_netcdf_file: &
    if( found ) then
      allocate( this%quantum_yield_parms( size( netcdfFiles ) ) )
file_loop: &
      do fileNdx = 1, size( netcdfFiles )
        allocate( netcdf_obj )
        ! read netcdf file quantum yield data
        call netcdf_obj%read_netcdf_file(                                     &
                     file_path = netcdfFiles( fileNdx )%to_char( ),           &
                     variable_name = Hdr )
        nParms = size( netcdf_obj%parameters, dim = 2 )
        if( nParms < 1 ) then
          write(msg,*) Iam//'File: ',                                         &
              trim( netcdfFiles( fileNdx )%to_char( ) ),                      &
              ' parameters array has < 1 parameter'
          call die_msg( 493253966, msg )
        endif
        ! interpolate from data to model wavelength grid
        if( allocated( netcdf_obj%wavelength ) ) then
          if( .not. allocated( this%quantum_yield_parms( fileNdx )%array ) )  &
              then
            allocate( this%quantum_yield_parms( fileNdx )%array(              &
                                                lambdaGrid%ncells_, nParms ) )
          endif
          do parmNdx = 1, nParms
            data_lambda    = netcdf_obj%wavelength
            data_parameter = netcdf_obj%parameters(:,parmNdx)
            call this%add_points( config, data_lambda, data_parameter )
            this%quantum_yield_parms( fileNdx )%array( :, parmNdx ) =         &
                interpolator%interpolate( x_target = lambdaGrid%edge_,        &
                                          x_source = data_lambda,             &
                                          y_source = data_parameter )
          enddo
        else
          this%quantum_yield_parms( fileNdx )%array = netcdf_obj%parameters
        endif
        if( allocated( netcdf_obj%temperature ) ) then
          this%quantum_yield_parms( fileNdx )%temperature =                   &
              netcdf_obj%temperature
        endif
        deallocate( netcdf_obj )
      enddo file_loop
    else has_netcdf_file
      ! check for quantum yield constant
      call config%get( 'constant value', quantum_yield_constant, Iam,         &
                       found = found )
      if( found ) then
        allocate( this%quantum_yield_parms(1) )
        allocate( this%quantum_yield_parms(1)%array( lambdaGrid%ncells_, 1 ) )
        this%quantum_yield_parms(1)%array(:,1) = quantum_yield_constant
      endif
    endif has_netcdf_file

    deallocate( lambdaGrid )

  end subroutine base_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function run( this, grid_warehouse, profile_warehouse )                     &
      result( quantum_yield )
    ! Calculates the quantum yield
    !
    ! Uses the interpolated first quantum yield parameter as the quantum
    ! yield.

    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(quantum_yield_t),    intent(in) :: this ! This :f:type:`~tuvx_quantum_yield/quantum_yield_t` calculator
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    real(dk), allocatable                    :: quantum_yield(:,:) ! Calculated quantum yield (height, wavelength) [unitless]

    ! Local variables
    character(len=*), parameter :: Iam = 'base quantum yield calculate'
    integer                     :: vertNdx
    class(grid_t),  pointer     :: zGrid => null( )
    real(dk),       allocatable :: wrkQuantumYield(:,:)

    zGrid => grid_warehouse%get_grid( "height", "km" )

    allocate( wrkQuantumYield(                                                &
      size( this%quantum_yield_parms(1)%array, dim = 1 ), zGrid%ncells_ + 1 ) )

    ! Just copy the lambda interpolated array
    do vertNdx = 1, zGrid%ncells_ + 1
      wrkQuantumYield( :, vertNdx ) =                                         &
          this%quantum_yield_parms(1)%array( :, 1 )
    enddo

    quantum_yield = transpose( wrkQuantumYield )

    deallocate( zGrid )

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine add_points( this, config, data_lambda, data_parameter )
    ! Add points to the cross-section gridded data based on configuration
    ! options

    use musica_assert,                 only : assert_msg, die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_util,                     only : add_point

    class(quantum_yield_t), intent(in)    :: this ! This :f:type:`~tuvx_quantum_yield/quantum_yield_t` calculator
    type(config_t),         intent(inout) :: config ! Quantum yield configuration data
    real(dk), allocatable,  intent(inout) :: data_lambda(:) ! Wavelength grid
    real(dk), allocatable,  intent(inout) :: data_parameter(:) ! Parameters (wavelength)

    real(dk), parameter :: rZERO = 0.0_dk
    real(dk), parameter :: rONE  = 1.0_dk
    real(dk), parameter :: deltax = 1.e-5_dk
    character(len=*), parameter :: Iam = 'cross_section; addpnts: '

    integer :: nRows
    real(dk) :: lowerLambda, upperLambda
    real(dk) :: addpnt_val
    type(string_t)  :: addpnt_type
    type(config_t)  :: extrap_config
    logical         :: found
    character(len=:), allocatable :: number
    type(string_t) :: required_keys(1), optional_keys(1)

    required_keys(1) = "type"
    optional_keys(1) = "value"

    nRows = size( data_lambda )
    lowerLambda = data_lambda(1) ; upperLambda = data_lambda( nRows )

    ! add endpoints to data arrays; first the lower bound
    addpnt_val = rZERO
    call config%get( 'lower extrapolation', extrap_config, Iam, found = found )
    if( found ) then
      call assert_msg( 231583250,                                             &
                       extrap_config%validate( required_keys, optional_keys ),&
                       "Bad format for extrapolation in base quantum yield." )
      call extrap_config%get( "type", addpnt_type, Iam )
      if( addpnt_type == 'boundary' ) then
        addpnt_val = data_parameter(1)
      elseif( addpnt_type == 'constant' ) then
        call extrap_config%get( "value", addpnt_val, Iam )
      else
        call die_msg( 163668372,                                              &
                      "Bad extrapolation type: '"//addpnt_type//"'" )
      endif
    endif
    call add_point( x = data_lambda, y = data_parameter, xnew = rZERO,        &
                    ynew = addpnt_val )
    call add_point( x = data_lambda, y = data_parameter,                      &
                    xnew = ( rONE - deltax ) * lowerLambda, ynew = addpnt_val )

    ! add endpoints to data arrays; now the upper bound
    addpnt_val = rZERO
    call config%get( 'upper extrapolation', extrap_config, Iam, found = found )
    if( found ) then
      call assert_msg( 104376574,                                             &
                       extrap_config%validate( required_keys, optional_keys ),&
                       "Bad format for extrapolation in base quantum yield." )
      call extrap_config%get( "type", addpnt_type, Iam )
      if( addpnt_type == 'boundary' ) then
        addpnt_val = data_parameter( nRows )
      elseif( addpnt_type == 'constant' ) then
        call extrap_config%get( "value", addpnt_val, Iam )
      else
        call die_msg( 216694919,                                              &
                      "Bad extrapolation type: '"//addpnt_type//"'" )
      endif
    endif
    call add_point( x = data_lambda, y = data_parameter,                      &
                    xnew = ( rONE + deltax ) * upperLambda, ynew = addpnt_val )
    call add_point( x = data_lambda, y = data_parameter, xnew = 1.e38_dk, &
                    ynew = addpnt_val )

  end subroutine add_points

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! Returns the number of bytes required to pack the object onto a buffer

    use musica_mpi,                    only : musica_mpi_pack_size

    class(quantum_yield_t), intent(in) :: this ! quantum yield to be packed
    integer, optional,      intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: i_param

    pack_size =                                                               &
      musica_mpi_pack_size( allocated( this%quantum_yield_parms ), comm )
    if( allocated( this%quantum_yield_parms ) ) then
      pack_size = pack_size +                                                 &
        musica_mpi_pack_size( size( this%quantum_yield_parms ), comm )
      do i_param = 1, size( this%quantum_yield_parms )
        pack_size = pack_size +                                               &
          this%quantum_yield_parms( i_param )%pack_size( comm )
      end do
    end if
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, position, comm )
    ! Packs the quantum yield onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack

    class(quantum_yield_t), intent(in)    :: this      ! quantum yield to pack
    character,              intent(inout) :: buffer(:) ! memory buffer
    integer,                intent(inout) :: position  ! current buffer position
    integer, optional,      intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_param

    prev_pos = position
    call musica_mpi_pack( buffer, position,                                   &
                          allocated( this%quantum_yield_parms ), comm )
    if( allocated( this%quantum_yield_parms ) ) then
      call musica_mpi_pack( buffer, position,                                 &
                            size( this%quantum_yield_parms ), comm )
      do i_param = 1, size( this%quantum_yield_parms )
        call this%quantum_yield_parms( i_param )%mpi_pack( buffer, position,  &
                                                           comm )
      end do
    end if
    call assert( 165656641, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! Unpacks a quantum yield from a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    class(quantum_yield_t), intent(out)   :: this      ! quantum yield to be unpacked
    character,              intent(inout) :: buffer(:) ! memory buffer
    integer,                intent(inout) :: position  ! current buffer position
    integer, optional,      intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_param, n_params
    logical :: alloced

    prev_pos = position
    call musica_mpi_unpack( buffer, position, alloced, comm )
    if( allocated( this%quantum_yield_parms ) )                               &
        deallocate( this%quantum_yield_parms )
    if( alloced ) then
      call musica_mpi_unpack( buffer, position, n_params, comm )
      allocate( this%quantum_yield_parms( n_params ) )
      do i_param = 1, n_params
        call this%quantum_yield_parms( i_param )%mpi_unpack( buffer,          &
                                                             position, comm )
      end do
    end if
    call assert( 865270779, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function parms_pack_size( this, comm ) result( pack_size )
    ! Returns the number of bytes required to pack the parameters onto a
    ! buffer.

    use musica_mpi,                    only : musica_mpi_pack_size

    class(quantum_yield_parms_t), intent(in) :: this ! parameters to pack
    integer, optional,            intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    pack_size = musica_mpi_pack_size( this%temperature, comm )                &
                + musica_mpi_pack_size( this%array, comm )
#else
    pack_size = 0
#endif

  end function parms_pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine parms_mpi_pack( this, buffer, position, comm )
    ! Packs the parameters onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack

    class(quantum_yield_parms_t), intent(in)    :: this      ! quantum yield to pack
    character,                    intent(inout) :: buffer(:) ! memory buffer
    integer,                      intent(inout) :: position  ! current buffer position
    integer, optional,            intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call musica_mpi_pack( buffer, position, this%temperature, comm )
    call musica_mpi_pack( buffer, position, this%array,       comm )
    call assert( 220419984, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine parms_mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine parms_mpi_unpack( this, buffer, position, comm )
    ! Unpacks parameters from a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    class(quantum_yield_parms_t), intent(out)   :: this      ! parameter to be unpacked
    character,                    intent(inout) :: buffer(:) ! memory buffer
    integer,                      intent(inout) :: position  ! current buffer position
    integer, optional,            intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call musica_mpi_unpack( buffer, position, this%temperature, comm )
    call musica_mpi_unpack( buffer, position, this%array,       comm )
    call assert( 203626119, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine parms_mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_quantum_yield
