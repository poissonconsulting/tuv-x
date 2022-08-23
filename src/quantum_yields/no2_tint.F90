! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_quantum_yield_no2_tint
  !> The no2 tint quantum yield type and related functions

  use musica_constants,                only : dk => musica_dk
  use tuvx_quantum_yield,              only : quantum_yield_t, base_constructor

  implicit none

  private
  public :: quantum_yield_no2_tint_t

  type quantum_yield_data_t
    real(dk), allocatable :: temperature(:) ! Temperature grid [K]
    real(dk), allocatable :: deltaT(:)      ! Temperature difference between grid points [K]
    real(dk), allocatable :: array(:,:)     ! Quantum yield parameters (wavelength, temperature)
  end type quantum_yield_data_t

  type, extends(quantum_yield_t) :: quantum_yield_no2_tint_t
    ! Calculator for tint quantum yield
    type(quantum_yield_data_t), allocatable :: quantum_yield(:)
  contains
    procedure :: calculate => run
    final     :: finalize
  end type quantum_yield_no2_tint_t

  interface quantum_yield_no2_tint_t
    module procedure constructor
  end interface quantum_yield_no2_tint_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result( this )
    ! Constructor

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_interpolate,              only : interpolator_conserving_t
    use tuvx_netcdf,                   only : netcdf_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    type(quantum_yield_no2_tint_t),  pointer :: this ! This :f:type:`~tuvx_quantum_yield_no2_tint/quantum_yield_no2_tint_t`
    type(config_t),            intent(inout) :: config ! Quantum yield configuration data
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`

    ! Local variables
    character(len=*), parameter :: Iam = 'no2 tint quantum yield constructor'
    character(len=*), parameter :: Hdr = 'quantum_yield_'
    real(dk), parameter    :: rZERO = 0.0_dk
    real(dk), parameter    :: rONE  = 1.0_dk

    integer     :: nmdlLambda
    integer     :: nTemps, nParms
    integer     :: parmNdx, fileNdx, Ndxl, Ndxu
    real(dk)    :: tmp
    real(dk)    :: quantum_yield_constant
    real(dk), allocatable :: data_lambda(:)
    real(dk), allocatable :: data_parameter(:)
    logical     :: found, monopos
    character(len=:), allocatable :: msg
    type(netcdf_t),   allocatable :: netcdf_obj
    type(string_t),   allocatable :: netcdfFiles(:)
    class(grid_t),    pointer     :: lambdaGrid => null( )
    type(interpolator_conserving_t) :: interpolator

    allocate( this )

    !> Get model wavelength grid
    lambdaGrid => grid_warehouse%get_grid( "wavelength", "nm" )

    ! get quantum yield netcdf filespec
    call config%get( 'netcdf files', netcdfFiles, Iam, found = found )
    call assert_msg( 252847407, found,                                        &
                     Iam//'must have at least one netcdf input file' )
    allocate( this%quantum_yield( size( netcdfFiles ) ) )
    file_loop: do fileNdx = 1, size( netcdfFiles )
      allocate( netcdf_obj )
      ! read netcdf file quantum yield data
      call netcdf_obj%read_netcdf_file(                                       &
                               file_path = netcdfFiles( fileNdx )%to_char( ), &
                               variable_name = Hdr )
      nParms = size( netcdf_obj%parameters, dim = 2 )
      call assert_msg( 235314124, nParms >= 2, Iam//'File: '//                &
                       trim( netcdfFiles( fileNdx )%to_char( ) )//            &
                       ' array must have 2 or more parameters' )
      associate( Qyield => this%quantum_yield( fileNdx ) )
      ! interpolation temperatures must be in netcdf file
      call assert_msg( 264376965, allocated( netcdf_obj%temperature ),        &
                       Iam//'File: '//                                        &
                       trim( netcdfFiles( fileNdx )%to_char( ) )//            &
                       ' does not have interpolation temperatures' )
      Qyield%temperature = netcdf_obj%temperature
      nTemps = size( Qyield%temperature )
      ! must have two or more interpolation temperatures
      call assert_msg( 393489175, nTemps >= 2, Iam//'File: '//                &
                       trim( netcdfFiles( fileNdx )%to_char( ) )//            &
                       ' temperature array has < 2 entries' )
      call assert_msg( 167399246, nTemps >= nParms, Iam//'File: '//           &
                       trim( netcdfFiles( fileNdx )%to_char( ) )//            &
                       ' temperature array < number parameters' )
      Qyield%deltaT = Qyield%temperature( 2 : nParms ) -                      &
                        Qyield%temperature( 1 : nParms - 1 )
      monopos = all( Qyield%deltaT > rZERO )
      if( .not. monopos ) then
        call assert_msg( 606144012, .not. any( Qyield%deltaT > rZERO ),       &
                         Iam//'File: '//                                      &
                         trim( netcdfFiles( fileNdx )%to_char( ) )//          &
                         ' temperature array not monotonic' )
        do Ndxl = 1, nParms / 2
          Ndxu = nParms - Ndxl + 1
          tmp = Qyield%temperature( Ndxl )
          Qyield%temperature( Ndxl ) = Qyield%temperature( Ndxu )
          Qyield%temperature( Ndxu ) = tmp
          data_parameter = netcdf_obj%parameters(:,Ndxl)
          netcdf_obj%parameters( :, Ndxl ) =                                  &
              netcdf_obj%parameters( :, Ndxu )
          netcdf_obj%parameters( :, Ndxu ) = data_parameter
        enddo
        Qyield%deltaT = Qyield%temperature( 2 : nParms ) -                    &
                          Qyield%temperature( 1 : nParms - 1 )
      endif
      ! interpolate from data to model wavelength grid
      if( allocated( netcdf_obj%wavelength ) ) then
        if( .not. allocated( this%quantum_yield( fileNdx )%array) ) then
          allocate( this%quantum_yield( fileNdx )%array(                      &
                                              lambdaGrid%ncells_, nParms ) )
        endif
        do parmNdx = 1, nParms
          data_lambda    = netcdf_obj%wavelength
          data_parameter = netcdf_obj%parameters( :, parmNdx )
          call this%add_points( config, data_lambda, data_parameter )
          this%quantum_yield( fileNdx )%array( :, parmNdx ) =                 &
                interpolator%interpolate( x_target = lambdaGrid%edge_,        &
                                          x_source = data_lambda,             &
                                          y_source = data_parameter )
        enddo
      else
        this%quantum_yield( fileNdx )%array = netcdf_obj%parameters
      endif
      end associate
      deallocate( netcdf_obj )
    enddo file_loop

    deallocate( lambdaGrid )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the quantum yield for the environmental conditions
  function run( this, grid_warehouse, profile_warehouse )                     &
      result( quantum_yield )

    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile,                  only : profile_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(quantum_yield_no2_tint_t), intent(in)    :: this ! This :f:type:`~tuvx_quantum_yield_no2_tint/quantum_yield_no2_tint_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    real(kind=dk), allocatable               :: quantum_yield(:,:) ! Calculated quantum_yield

    !> Local variables
    character(len=*), parameter :: Iam = 'no2 tint quantum yield calculate'
    integer     :: nTemp
    integer     :: fileNdx, tNdx, vertNdx
    real(dk)    :: Tadj, Tstar
    real(dk),         allocatable :: WrkQuantumYield(:,:)
    class(grid_t),    pointer     :: zGrid => null( )
    class(grid_t),    pointer     :: lambdaGrid => null( )
    class(profile_t), pointer     :: mdlTemperature => null( )

    zGrid => grid_warehouse%get_grid( "height", "km" )
    lambdaGrid => grid_warehouse%get_grid( "wavelength", "nm" )
    mdlTemperature => profile_warehouse%get_profile( "temperature", "K" )

    allocate( wrkQuantumYield( lambdaGrid%ncells_, zGrid%ncells_ + 1 ) )
    wrkQuantumYield = 0.0_dk

    do fileNdx = 1, size( this%quantum_yield )
      associate( Temp => this%quantum_yield( fileNdx )%temperature,           &
                 wrkQyield => this%quantum_yield( fileNdx ) )
      nTemp = size( Temp )
      do vertNdx = 1, zGrid%ncells_ + 1
        Tadj  = mdlTemperature%edge_val_( vertNdx )
        do tNdx = 2, nTemp
          if( Tadj <= Temp( tNdx ) ) then
            exit
          endif
        enddo
        tndx = min( nTemp, tNdx ) - 1
        Tstar = ( Tadj - Temp( tNdx ) ) / wrkQyield%deltaT( tNdx )
        WrkQuantumYield( :, vertNdx ) = WrkQuantumYield( :, vertNdx ) +       &
                    wrkQyield%array( :, tNdx ) +                              &
                    Tstar * ( wrkQyield%array( :, tNdx + 1 ) -             &
                              wrkQyield%array( :, tNdx ) )
      enddo
      end associate
    enddo

    quantum_yield = transpose( max( WrkQuantumYield, 0.0_dk ) )

    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( mdlTemperature )

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! finalize the quantum yield type

    type(quantum_yield_no2_tint_t), intent(inout) :: this

    integer     :: ndx

    if( allocated( this%quantum_yield ) ) then
      do ndx = 1, size( this%quantum_yield )
        associate( Qyield => this%quantum_yield( ndx ) )
        if( allocated( Qyield%array ) ) then
          deallocate( Qyield%array )
        endif
        if( allocated( Qyield%temperature ) ) then
         deallocate( Qyield%temperature )
        endif
        if( allocated( Qyield%deltaT ) ) then
          deallocate( Qyield%deltaT )
        endif
        end associate
      enddo
      deallocate( this%quantum_yield )
    endif

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_quantum_yield_no2_tint
