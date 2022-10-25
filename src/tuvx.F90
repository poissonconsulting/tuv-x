! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
program tuvx
  ! Driver for stand-alone TUV-x and regression tests
  !
  ! If MPI support is compiled in, core driver constructs the core
  ! on the primary process, passes it to the other processes, and
  ! performs calculations on process 1. This allows MPI functions to
  ! be tested as part of the regression tests.

  use musica_assert,                   only : assert
  use musica_constants,                only : dk => musica_dk
  use musica_string,                   only : string_t
  use musica_mpi
  use tuvx_core,                       only : core_t

  implicit none

  class(core_t), pointer :: core

  ! Command-line options
  character(len=256) :: argument
  type(string_t)     :: configFileSpec
  character, allocatable :: buffer(:)
  integer :: pos, pack_size
  integer, parameter :: comm = MPI_COMM_WORLD

  call musica_mpi_init( )

  ! Get the model configuration file and options from the command line
  if( command_argument_count() /= 1 ) then
    call fail_run( )
  endif
  call get_command_argument( 1, argument )

  configFileSpec = argument

  ! instatiate and initialize photolysis core object on the
  ! primary MPI process
  if( musica_mpi_rank( comm ) == 0 ) then
    core => core_t( configFileSpec )
    pack_size = core%pack_size( comm )
    allocate( buffer( pack_size ) )
    pos = 0
    call core%mpi_pack( buffer, pos, comm )
    call assert( 115315474, pos <= pack_size )
  end if

  ! broadcast the core data to other MPI processes
  call musica_mpi_bcast( pack_size, comm )
  if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
  call musica_mpi_bcast( buffer, comm )

  ! Unpack the core on other MPI processes
  if( musica_mpi_rank( comm ) .ne. 0 ) then
    pos = 0
    allocate( core )
    call core%mpi_unpack( buffer, pos, comm )
    call assert( 501938283, pos <= pack_size )
  end if

  ! Perform photolysis calculations on MPI process 1
  ! if MPI support is compiled in
  if( musica_mpi_size( comm ) == 1 .or. musica_mpi_rank( comm ) == 1 ) then
    call run_tuvx( core )
  end if

  deallocate( buffer )
  deallocate( core )

  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine run_tuvx( core )
    ! Runs TUV-x under SZA conditions provided in the configuration

    use tuvx_grid,                     only : grid_t
    use tuvx_profile,                  only : profile_t

    type(core_t), intent(inout) :: core

    integer                   :: i_sza
    class(grid_t),    pointer :: height
    class(profile_t), pointer :: sza
    real(dk), allocatable     :: photo_rates(:,:,:) ! (time, vertical level, reaction) [s-1]
    real(dk), allocatable     :: dose_rates(:,:,:)  ! (time, vertical level, dose rate) [?]
    type(string_t)            :: file_path
    character(len=2)          :: diagnostic_label

    height => core%get_grid( "height", "km" )
    sza => core%get_profile( "solar zenith angle", "degrees" )

    allocate( photo_rates( sza%ncells_ + 1,                                   &
                           height%ncells_ + 1,                                &
                           core%number_of_photolysis_reactions( ) ) )
    allocate( dose_rates(           sza%ncells_ + 1,                          &
                                    height%ncells_ + 1,                       &
                                    core%number_of_dose_rates( ) ) )
    ! calculate photolysis and dose rates
    do i_sza = 1, sza%ncells_ + 1
      write(diagnostic_label,'(i2.2)') i_sza
      call core%run( sza%edge_val_( i_sza ),                                  &
                     photolysis_rate_constants = photo_rates( i_sza, :, : ),  &
                     dose_rates = dose_rates( i_sza, :, : ),                  &
                     diagnostic_label = diagnostic_label )
    end do

    ! output results
    if( core%number_of_photolysis_reactions( ) > 0 ) then
      file_path = "photolysis_rate_constants.nc"
      call output_photolysis_rate_constants( core, photo_rates, file_path )
    end if
    if( core%number_of_dose_rates( ) > 0 ) then
      file_path = "dose_rates.nc"
      call output_dose_rates( core, dose_rates, file_path )
    end if

    deallocate( height )
    deallocate( sza    )

  end subroutine run_tuvx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine output_photolysis_rate_constants( core, values, file_path )
    ! Outputs calculated photolysis rate constants

    use musica_assert,                 only : assert
    use musica_io,                     only : io_t
    use musica_io_netcdf,              only : io_netcdf_t
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_netcdf,                   only : clean_string
    use tuvx_profile,                  only : profile_t

    type(core_t),   intent(in) :: core          ! TUV-x core
    real(dk),       intent(in) :: values(:,:,:) ! Photolysis rate constants (time, vertical level, reaction)
    type(string_t), intent(in) :: file_path     ! File path to output to

    character(len=*), parameter :: Iam = "photolysis output"
    class(io_t),        pointer :: out_file
    integer                     :: i_rxn
    type(string_t)              :: var_name, dim_names(2), units
    type(string_t), allocatable :: rxn_names(:)
    class(profile_t),   pointer :: sza
    class(grid_t),      pointer :: time, vertical
    integer                     :: stat

    sza => core%get_profile( "solar zenith angle", "degrees" )
    time => core%get_grid( "time", "hours" )
    vertical => core%get_grid( "height", "km" )
    rxn_names = core%photolysis_reaction_labels( )
    call assert( 182934700, size( sza%edge_val_ ) .eq. size( time%edge_ ) )
    call assert( 394136298, size( values, 1 ) .eq. size( time%edge_ ) )
    call assert( 664629694, size( values, 2 ) .eq. size( vertical%edge_ ) )
    call assert( 266929622, size( values, 3 ) .eq. size( rxn_names ) )

    ! Remove any existing file with the same name
    open( unit = 16, iostat = stat, file = file_path%to_char( ),              &
          status = 'old' )
    if( stat == 0 ) close( 16, status = 'delete' )

    out_file => io_netcdf_t( file_path )

    var_name = "altitude"
    dim_names(1) = "vertical_level"
    units = "km"
    call out_file%write( var_name, dim_names(1), vertical%edge_, Iam )
    call out_file%set_variable_units( var_name, units, Iam )

    var_name = "time"
    dim_names(1) = "time"
    units = "hr"
    call out_file%write( var_name, dim_names(1), time%edge_, Iam )
    call out_file%set_variable_units( var_name, units, Iam )

    var_name = "solar zenith angle"
    dim_names(1) = "time"
    units = "degrees"
    call out_file%write( var_name, dim_names(1), sza%edge_val_, Iam )
    call out_file%set_variable_units( var_name, units, Iam )

    dim_names(1) = "time"
    dim_names(2) = "vertical_level"
    units = "s-1"
    do i_rxn = 1, size( rxn_names )
      var_name = clean_string( rxn_names( i_rxn ) )
      call out_file%write( var_name, dim_names, values( :, :, i_rxn ), Iam )
      call out_file%set_variable_units( var_name, units, Iam )
    end do

    deallocate( sza )
    deallocate( time )
    deallocate( vertical )
    deallocate( out_file )

  end subroutine output_photolysis_rate_constants

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine output_dose_rates( core, values, file_path )
    ! Outputs calculated dose rates

    use musica_assert,                 only : assert
    use musica_io,                     only : io_t
    use musica_io_netcdf,              only : io_netcdf_t
    use tuvx_grid,                     only : grid_t
    use tuvx_netcdf,                   only : clean_string
    use tuvx_profile,                  only : profile_t

    type(core_t),   intent(in) :: core          ! TUV-x core
    real(dk),       intent(in) :: values(:,:,:) ! Dose rates (time, vertical level, dose rate type)
    type(string_t), intent(in) :: file_path     ! File path to output to

    character(len=*), parameter :: Iam = "photolysis output"
    class(io_t),        pointer :: out_file
    integer                     :: i_rate
    type(string_t)              :: var_name, dim_names(2), units
    type(string_t), allocatable :: rate_names(:)
    class(profile_t),   pointer :: sza
    class(grid_t),      pointer :: time, vertical
    integer                     :: stat

    sza => core%get_profile( "solar zenith angle", "degrees" )
    time => core%get_grid( "time", "hours" )
    vertical => core%get_grid( "height", "km" )
    rate_names = core%dose_rate_labels( )
    call assert( 182934700, size( sza%edge_val_ ) .eq. size( time%edge_ ) )
    call assert( 394136298, size( values, 1 ) .eq. size( time%edge_ ) )
    call assert( 664629694, size( values, 2 ) .eq. size( vertical%edge_ ) )
    call assert( 266929622, size( values, 3 ) .eq. size( rate_names ) )

    ! Remove any existing file with the same name
    open( unit = 16, iostat = stat, file = file_path%to_char( ),              &
          status = 'old' )
    if( stat == 0 ) close( 16, status = 'delete' )

    out_file => io_netcdf_t( file_path )

    var_name = "altitude"
    dim_names(1) = "vertical_level"
    units = "km"
    call out_file%write( var_name, dim_names(1), vertical%edge_, Iam )
    call out_file%set_variable_units( var_name, units, Iam )

    var_name = "time"
    dim_names(1) = "time"
    units = "hr"
    call out_file%write( var_name, dim_names(1), time%edge_, Iam )
    call out_file%set_variable_units( var_name, units, Iam )

    var_name = "solar zenith angle"
    dim_names(1) = "time"
    units = "degrees"
    call out_file%write( var_name, dim_names(1), sza%edge_val_, Iam )
    call out_file%set_variable_units( var_name, units, Iam )

    dim_names(1) = "time"
    dim_names(2) = "vertical_level"
    do i_rate = 1, size( rate_names )
      var_name = clean_string( rate_names( i_rate ) )
      call out_file%write( var_name, dim_names, values( :, :, i_rate ), Iam )
    end do

    deallocate( sza )
    deallocate( time )
    deallocate( vertical )
    deallocate( out_file )

  end subroutine output_dose_rates

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fail_run( )
    ! Fail run and print usage info

    write(*,*) "Usage: ./photolysis configuration_file.json"
    stop 3

  end subroutine fail_run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program tuvx
