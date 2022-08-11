! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests the photo_decomp Profile module

!> Test module for the Profile_t type
program test_Profile

  use musica_string,    only : string_t

  implicit none

  !> Command-line options
  character(len=256) :: argument
  type(string_t)     :: configFileSpec
  !> Command-line argument index
  integer :: i_arg

  !> Get the model configuration file and options from the command line
  argument = 'test/data/profile.test.config.json'

  configFileSpec = argument
  call test_Profile_t( configFileSpec )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Tests the Profile_t type
  subroutine test_Profile_t( config_flsp )

    use musica_config,    only : config_t
    use musica_string,    only : string_t
    use musica_assert,    only : assert, almost_equal
    use musica_constants, only : ik => musica_ik, dk => musica_dk
    use tuvx_diagnostic_util,          only : prepare_diagnostic_output
    use tuvx_grid_warehouse, only : grid_warehouse_t
    use tuvx_grid,    only : grid_t
    use tuvx_profile_warehouse, only : Profile_warehouse_t
    use tuvx_profile,           only : profile_t

    !> Arguments
    type(string_t), intent(in) :: config_flsp
    !> Local variables
    character(len=*), parameter :: Iam = 'test_Profile: '
    type(config_t)              :: tst_config, child_config
    type(grid_warehouse_t), pointer :: theGridWarehouse
    class(grid_t), pointer   :: zGrid, lambdaGrid
    type(Profile_warehouse_t), pointer :: theProfileWarehouse
    class(profile_t), pointer      :: aProfile
    class(profile_t), pointer      :: AirProfile, TemperatureProfile
    class(profile_t), pointer      :: O3Profile
    type(string_t)                  :: Handle

    write(*,*) Iam // 'entering'

    call prepare_diagnostic_output( )

    !> master configuration -> config type
    call tst_config%from_file( config_flsp%to_char() )

    !> Initialize grid warehouse
    call tst_config%get( "grids", child_config, Iam )
    theGridWarehouse => grid_warehouse_t( child_config )

    !> Get copy of grid
    zGrid => theGridWarehouse%get_grid( "height", "km" )
    call assert( 412238768, zGrid%ncells_ .eq. 120_ik )
    call assert( 412238769, all( zGrid%delta_ .eq. 1._dk ) )

    !> Get copy of wavelength grid
    lambdaGrid => theGridWarehouse%get_grid( "wavelength", "nm" )

    !> Initialize profile warehouse
    call tst_config%get( "profiles", child_config, Iam )
    theProfileWarehouse =>                                                    &
        Profile_warehouse_t( child_config, theGridWareHouse )

    !> Get copy of the Air Profile
    AirProfile => theProfileWarehouse%get_profile( "air", "molecule cm-3" )
    call assert( 412238771, all( AirProfile%delta_val_ < 0._dk ) )

    !> Get copy of the temperature Profile
    TemperatureProfile => theProfileWarehouse%get_profile( "temperature", "K" )
    call assert( 412238772, all( TemperatureProfile%edge_val_ < 400._dk ) )
    call assert( 412238772, all( TemperatureProfile%edge_val_ > 150._dk ) )
    call assert( 412238773, all( abs(TemperatureProfile%delta_val_) < 20._dk ) )

    !> Get copy of the Unit test
    aProfile => theProfileWarehouse%get_profile( "UnitTest", "foos" )
    write(*,*) Iam // 'UnitTest'
    write(*,'(1p10g15.7)') aProfile%mid_val_
    call assert( 412238770, all( aProfile%mid_val_ .eq. 7.7e11_dk ) )

    deallocate( aProfile )

    !> Get copy of O2 profile
    aProfile => theProfileWarehouse%get_profile( "O2", "molecule cm-3" )

    deallocate( aProfile )

    !> Get copy of ozone profile
    O3Profile => theProfileWarehouse%get_profile( "O3", "molecule cm-3" )

    write(*,*) Iam // 'Handle = ',O3Profile%handle_
    write(*,*) ' '
    write(*,*) ' '
    write(*,*) Iam // 'O3 on model z grid'
    write(*,'(1p10g15.7)') O3Profile%edge_val_

    deallocate( O3Profile )

    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( AirProfile )
    deallocate( TemperatureProfile )
    deallocate( theGridWarehouse )
    deallocate( theProfileWarehouse )
    write(*,*) Iam // 'leaving'

  end subroutine test_Profile_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Fail run and print usage info
  subroutine fail_run( )

    write(*,*) "Usage: ./Profile_test configuration_file.json"
    stop 3

  end subroutine fail_run

end program test_Profile
