! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

program test_profile_warehouse
  ! Test module for the profile warehouse

  use musica_mpi,                      only : musica_mpi_init,                &
                                              musica_mpi_finalize

  implicit none

  call musica_mpi_init( )
  call profile_warehouse_test( )
  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine profile_warehouse_test()
    use musica_assert,              only : assert, almost_equal
    use musica_config,              only : config_t
    use musica_mpi
    use tuvx_grid_warehouse,        only : grid_warehouse_t
    use tuvx_profile_warehouse,     only : profile_warehouse_t
    use tuvx_profile,               only : profile_t

    character(len=*), parameter :: grid_config = 'test/data/grid.test.config.json'
    character(len=*), parameter :: profile_config = 'test/data/profile.temperature.config.json'
    type(config_t) :: grid_tst_config
    type(config_t) :: profile_tst_config
    class(grid_warehouse_t), pointer :: grid_warehouse => null()
    class(profile_warehouse_t), pointer :: profile_warehouse => null()
    class(profile_t),           pointer :: profile_ptr => null()
    character, allocatable :: buffer(:)
    integer :: pos, pack_size
    integer, parameter :: comm = MPI_COMM_WORLD

    call grid_tst_config%from_file( grid_config )
    grid_warehouse => grid_warehouse_t( grid_tst_config )

    if( musica_mpi_rank( comm ) == 0 ) then
      call profile_tst_config%from_file( profile_config )
      profile_warehouse =>                                                    &
          profile_warehouse_t( profile_tst_config, grid_warehouse )
      pack_size = profile_warehouse%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call profile_warehouse%mpi_pack( buffer, pos , comm )
      call assert( 893633789, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size , comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer , comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      allocate( profile_warehouse )
      call profile_warehouse%mpi_unpack( buffer, pos , comm )
      call assert( 436189624, pos <= pack_size )
    end if
    deallocate( buffer )

    profile_ptr => profile_warehouse%get_profile( "temperature", "K" )

    call assert( 418832741, associated(profile_ptr) )

    deallocate( profile_ptr )
    deallocate( grid_warehouse )
    deallocate( profile_warehouse )

  end subroutine profile_warehouse_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_profile_warehouse
