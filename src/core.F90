! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> The photolysis_core_t and related procedures
module tuvx_core

  use musica_config,                   only : config_t
  use musica_string,                   only : string_t
  use musica_assert,                   only : assert
  use musica_constants,                only : dk => musica_dk
  use tuvx_grid_warehouse,             only : grid_warehouse_t
  use tuvx_grid,                       only : grid_t
  use tuvx_profile_warehouse,          only : profile_warehouse_t
  use tuvx_spherical_geometry,         only : spherical_geometry_t
  use tuvx_la_sr_bands,                only : la_sr_bands_t
  use tuvx_radiative_transfer,         only : radiative_transfer_t
  use tuvx_photolysis_rates,           only : photolysis_rates_t
  use tuvx_dose_rates,                 only : dose_rates_t

  implicit none

  private
  public :: photolysis_core_t

  type :: photolysis_core_t
    type(grid_warehouse_t),      pointer :: grid_warehouse_ => null()
    type(profile_warehouse_t),   pointer :: profile_warehouse_ => null()
    type(spherical_geometry_t),  pointer :: spherical_geometry_ => null()
    type(la_sr_bands_t),         pointer :: la_sr_bands_ => null()
    type(string_t),          allocatable :: diagnostics_(:)
    type(radiative_transfer_t),  pointer :: radiative_transfer_ => null()
    type(photolysis_rates_t),    pointer :: photolysis_rates_ => null()
    type(dose_rates_t),          pointer :: dose_rates_ => null()
  contains
    procedure :: run
    procedure :: output
    final     :: finalize
  end type photolysis_core_t

  interface photolysis_core_t
    module procedure constructor
  end interface photolysis_core_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config ) result( new_core )

    use musica_assert,                 only : assert_msg
    use musica_iterator,               only : iterator_t
    use musica_string,                 only : string_t
    use tuvx_diagnostic_util,          only : diagout
    use tuvx_profile,                  only : profile_t

    !> Configuration data
    type(string_t), intent(in) :: config

    class(photolysis_core_t), pointer :: new_core

    ! Local variables
    character(len=*), parameter :: Iam = 'Photolysis core constructor: '
    character(len=32)           :: keyChar
    logical                     :: found
    type(string_t)              :: keyString
    type(config_t)              :: core_config, child_config
    class(iterator_t), pointer  :: iter
    class(profile_t),  pointer  :: aprofile
    type(string_t)              :: required_keys(4), optional_keys(2)

    call core_config%from_file( config%to_char() )

    ! Check json configuration file for basic structure, integrity
    required_keys(1) = "radiative transfer"
    required_keys(2) = "grids"
    required_keys(3) = "profiles"
    required_keys(4) = "O2 absorption"
    optional_keys(1) = "photolysis reactions"
    optional_keys(2) = "dose rates"
    call assert_msg( 255400232,                                               &
                     core_config%validate( required_keys, optional_keys ),    &
                     "Bad configuration data format for tuv-x core." )

    ! Instantiate photolysis core
    allocate( new_core )

    ! Instantiate and initialize grid warehouse
    call core_config%get( "grids", child_config, Iam )
    new_core%grid_warehouse_ => grid_warehouse_t( child_config )

    ! Instantiate and initialize profile warehouse
    call core_config%get( "profiles", child_config, Iam )
    new_core%profile_warehouse_ =>                                            &
       profile_warehouse_t( child_config, new_core%grid_warehouse_ )

    ! Diagnostics for testing
    aprofile => new_core%profile_warehouse_%get_profile( "temperature", "K" )
    call diagout( 'vptmp.new', aprofile%edge_val_ )
    deallocate( aprofile )

    aprofile => new_core%profile_warehouse_%get_profile( "air",               &
                                                         "molecule cm-3" )
    call diagout( 'vpair.new', aprofile%edge_val_ )
    deallocate( aprofile )

    if( new_core%profile_warehouse_%exists( "O3", "molecule cm-3" ) ) then
      aprofile => new_core%profile_warehouse_%get_profile( "O3",              &
                                                           "molecule cm-3" )
      call diagout( 'vpco3.new', aprofile%layer_dens_ )
      deallocate( aprofile )
    end if

    ! Set up radiative transfer calculator
    call core_config%get( "radiative transfer", child_config, Iam )
    new_core%radiative_transfer_ => &
        radiative_transfer_t( child_config,                                   &
                                  new_core%grid_warehouse_,                   &
                                  new_core%profile_warehouse_ )

    ! get optical depth diagnostics to output
    !> \todo this should be moved out of the radiative transfer config if it
    !!       is owned by the core
    call child_config%get( "diagnostics", new_core%diagnostics_,              &
                           Iam, found = found )
    if( .not. found ) allocate( new_core%diagnostics_( 0 ) )

    ! photolysis rate constants
    call core_config%get( "photolysis reactions", child_config, Iam,          &
                          found = found )
    if( found ) then
      new_core%photolysis_rates_ => &
          photolysis_rates_t( child_config,                                   &
                              new_core%grid_warehouse_,                       &
                              new_core%profile_warehouse_ )
    end if

    ! dose rates
    call core_config%get( "dose rates", child_config, Iam, found = found )
    if( found ) then
      new_core%dose_rates_ => &
          dose_rates_t( child_config, new_core%grid_warehouse_,               &
                        new_core%profile_warehouse_ )
    end if

    ! instantiate and initialize spherical geometry type
    new_core%spherical_geometry_ =>                                           &
        spherical_geometry_t( new_core%grid_warehouse_ )

    ! instantiate and initialize lyman alpha, srb type
    call core_config%get( "O2 absorption", child_config, Iam )
    new_core%la_sr_bands_ => la_sr_bands_t( child_config,                     &
                                            new_core%grid_warehouse_ )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Perform calculations for specified photolysis and dose rates for a
  !! given set of conditions
  subroutine run( this )

    use tuvx_profile,                    only : profile_t
    use tuvx_radiator_warehouse,         only : radiator_warehouse_t
    use tuvx_radiator,                   only : radiator_t
    use tuvx_radiative_transfer_solver,  only : radiation_field_t
    use tuvx_diagnostic_util,            only : diagout

    !> Photolysis core
    class(photolysis_core_t), intent(inout)  :: this

    ! Local variables
    character(len=*), parameter       :: Iam = 'Photolysis core run: '
    integer                           :: i_ndx, i_diag
    ! photolysis rate constants (time, vertical level, reaction)
    real(dk), allocatable             :: all_photo_rates(:,:,:)
    real(dk), allocatable             :: photo_rates(:,:)
    real(dk), allocatable             :: dose_rates(:,:)
    character(len=2)                  :: number
    class(profile_t),         pointer :: solar_zenith_angles => null( )
    class(radiator_t),        pointer :: radiator => null()
    class(radiation_field_t), pointer :: radiation_field => null( )

    ! get the solar zenith angles
    solar_zenith_angles =>                                                    &
        this%profile_warehouse_%get_profile( "solar zenith angle", "degrees" )

    ! calculate the radiation field
    sza_loop: do i_ndx = 1,size(solar_zenith_angles%edge_val_)
      if( associated( this%spherical_geometry_ ) ) then
        call this%spherical_geometry_%set_parameters(                         &
            solar_zenith_angles%edge_val_(i_ndx), this%grid_warehouse_ )
      endif
      call this%radiative_transfer_%calculate( this%la_sr_bands_,             &
                                               this%spherical_geometry_,      &
                                               this%grid_warehouse_,          &
                                               this%profile_warehouse_,       &
                                               radiation_field )
      write(number,'(i2.2)') i_ndx
      call diagout( 'radField.' // number // '.new',                          &
                    radiation_field%fdr_ + radiation_field%fup_ +             &
                      radiation_field%fdn_ )
      if( associated( this%photolysis_rates_ ) ) then
        if( allocated( photo_rates ) ) then
          deallocate( photo_rates )
        endif
        call this%photolysis_rates_%get( this%la_sr_bands_,                   &
                                         this%spherical_geometry_,            &
                                         this%grid_warehouse_,                &
                                         this%profile_warehouse_,             &
                                         radiation_field,                     &
                                         photo_rates,                         &
                                         number )
        if( .not. allocated( all_photo_rates ) ) then
          allocate( all_photo_rates( size( solar_zenith_angles%edge_val_ ),   &
                                     size( photo_rates, 1 ),                  &
                                     size( photo_rates, 2 ) ) )
        end if
        all_photo_rates( i_ndx, :, : ) = photo_rates(:,:)
      elseif( associated(this%dose_rates_) ) then
        if( allocated(dose_rates) ) then
          deallocate(dose_rates)
        endif
        call this%dose_rates_%get( this%grid_warehouse_,                      &
                                   this%profile_warehouse_,                   &
                                   radiation_field,                           &
                                   dose_rates,                                &
                                   number )
      endif
      deallocate( radiation_field )
    enddo sza_loop

    ! output photolysis rate constants
    if( associated( this%photolysis_rates_ ) ) then
      call this%output( all_photo_rates, "photolysis_rate_constants.nc" )
    end if

    ! diagnostic output
    do i_diag = 1, size( this%diagnostics_ )
      associate( diagnostic => this%diagnostics_( i_diag ) )
        radiator =>                                                           &
            this%radiative_transfer_%radiator_warehouse_%get_radiator(        &
                                                                  diagnostic )
        ! Diagnostics for testing
        if( diagnostic == 'air' ) then
          call diagout( 'dtrl.new', radiator%state_%layer_OD_ )
        elseif( diagnostic == 'Aerosols' ) then
          call diagout( 'dtaer.new', radiator%state_%layer_OD_ )
        elseif( diagnostic == 'O3' ) then
          call diagout( 'dto3.new', radiator%state_%layer_OD_ )
        elseif( diagnostic == 'O2' ) then
          call diagout( 'dto2.new', radiator%state_%layer_OD_ )
        endif
      end associate
    end do

    deallocate( solar_zenith_angles )

  end subroutine run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Outputs calculated photolysis rate constants
  subroutine output( this, values, file_path )

    use musica_assert,                 only : assert
    use nc4fortran,                    only : netcdf_file
    use tuvx_grid,                     only : grid_t
    use tuvx_profile,                  only : profile_t

    !> TUV-x core
    class(photolysis_core_t), intent(in) :: this
    !> Photolysis rate constants (time, vertical level, reaction)
    real(dk),                 intent(in) :: values(:,:,:)
    !> File path to output to
    character(len=*),         intent(in) :: file_path

    type(netcdf_file)           :: output_file
    integer                     :: i_rxn
    type(string_t), allocatable :: rxn_names(:)
    class(profile_t),   pointer :: sza
    class(grid_t),      pointer :: time, vertical

    call assert( 337750978, associated( this%photolysis_rates_ ) )
    sza => this%profile_warehouse_%get_profile( "solar zenith angle",         &
                                                "degrees" )
    time => this%grid_warehouse_%get_grid( "time", "hours" )
    vertical => this%grid_warehouse_%get_grid( "height", "km" )
    rxn_names = this%photolysis_rates_%labels( )
    call assert( 182934700,                                                   &
                 size( sza%edge_val_ ) .eq. size( time%edge_ ) )
    call assert( 394136298,                                                   &
                 size( values, 1 ) .eq. size( time%edge_ ) )
    call assert( 664629694,                                                   &
                 size( values, 2 ) .eq. size( vertical%edge_ ) )
    call assert( 266929622,                                                   &
                 size( values, 3 ) .eq. size( rxn_names ) )
    call output_file%open( file_path, action='w' )
    call output_file%write( "altitude", vertical%edge_,                       &
                            (/ "vertical_level" /) )
    call output_file%write_attribute( "altitude", "units", "km" )
    call output_file%write( "time", time%edge_, (/ "time" /) )
    call output_file%write_attribute( "time", "units", "hr" )
    call output_file%write( "solar zenith angle", sza%edge_val_, (/ "time" /) )
    call output_file%write_attribute( "solar zenith angle", "units",          &
                                      "degrees" )
    do i_rxn = 1, size( rxn_names )
      call output_file%write( rxn_names( i_rxn )%val_, values( :, :, i_rxn ), &
                              (/ "time          ", "vertical_level" /) )
    end do
    call output_file%close( )
    deallocate( sza )
    deallocate( time )
    deallocate( vertical )

  end subroutine output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the photolysis core
  subroutine finalize( this )

    !> Photolysis core
    type(photolysis_core_t), intent(inout) :: this

    if( associated( this%grid_warehouse_ ) ) then
      deallocate( this%grid_warehouse_ )
    end if
    if( associated( this%profile_warehouse_ ) ) then
      deallocate( this%profile_warehouse_ )
    end if
    if( associated( this%spherical_geometry_ ) ) then
      deallocate( this%spherical_geometry_ )
    end if
    if( associated( this%la_sr_bands_ ) ) then
      deallocate( this%la_sr_bands_ )
    end if
    if( associated( this%radiative_transfer_ ) ) then
      deallocate( this%radiative_transfer_ )
    end if
    if( associated( this%photolysis_rates_ ) ) then
      deallocate( this%photolysis_rates_ )
    end if
    if( associated( this%dose_rates_ ) ) then
      deallocate( this%dose_rates_ )
    end if

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_core
