
! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

program test_spectral_weight_mpi_pack
  use musica_constants,       only : dk => musica_dk
  use musica_config,          only : config_t
  use tuvx_grid_warehouse,    only : grid_warehouse_t
  use tuvx_profile_warehouse, only : profile_warehouse_t
  use tuvx_spectral_weight,   only : spectral_weight_t
  use musica_assert,          only : assert
  use musica_mpi

  implicit none

  character(len=*), parameter :: conf = 'test/data/spectral_weights.json'
  type(config_t)                  :: grid_config, config, profile_config
  type(config_t)                  :: weights_config
  class(grid_warehouse_t), pointer :: grid_warehouse => null()
  class(spectral_weight_t), pointer :: spectral_weight => null()
  class(profile_warehouse_t), pointer :: profile_warehouse => null()

  call config%from_file( conf )

  call config%get( "grids", grid_config, "" )
  call config%get( "profiles", profile_config, "" )
  call config%get( "weights", weights_config, "" )

  grid_warehouse => grid_warehouse_t( grid_config )
  profile_warehouse => profile_warehouse_t( profile_config, grid_warehouse )
  spectral_weight => spectral_weight_t( weights_config, grid_warehouse,       &
    profile_warehouse)
  
  call musica_mpi_init( )

  call test_mpi_pack( spectral_weight )

  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_mpi_pack( the_weight )
    use tuvx_test_utils,  only : check_values
    use musica_constants, only : dk => musica_dk

    class(spectral_weight_t), pointer :: the_weight
    class(spectral_weight_t), pointer :: unpacked

    character, allocatable :: buffer(:)
    integer :: pos, pack_size

    ! Get copy of the rayliegh radiator and test MPI functions
    if( musica_mpi_rank( ) == 0 ) then
      pack_size = the_weight%pack_size( )

      allocate( buffer( pack_size ) )
      pos = 0
      call the_weight%mpi_pack( buffer, pos )
      call assert( 192131787, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size )
    if( musica_mpi_rank( ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer )

    if( musica_mpi_rank( ) .ne. 0 ) then
      pos = 0
      allocate( unpacked )
      call unpacked%mpi_unpack( buffer, pos )
      call assert( 862314213, pos > 0 )
      call check_values(                                                      &
        unpacked%spectral_weight_parms(1)%array,                              &
        the_weight%spectral_weight_parms(1)%array,                            &
        0.01_dk )
    end if


  end subroutine test_mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_spectral_weight_mpi_pack
