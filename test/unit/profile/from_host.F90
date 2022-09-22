! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
program test_profile_from_host

  use musica_mpi,                      only : musica_mpi_init,                &
                                              musica_mpi_finalize

  implicit none

  call musica_mpi_init( )
  call test_profile_from_host_t( )
  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_profile_from_host_t( )

    use musica_assert,                 only : assert, almost_equal, die
    use musica_constants,              only : dk => musica_dk
    use musica_string,                 only : string_t
    use musica_mpi
    use tuvx_profile,                  only : profile_t
    use tuvx_profile_from_host,        only : profile_from_host_t,            &
                                              profile_updater_t
    use tuvx_profile_factory,          only : profile_type_name,              &
                                              profile_allocate

    class(profile_t), pointer :: my_profile
    type(profile_updater_t) :: my_updater
    character, allocatable :: buffer(:)
    integer :: pos, pack_size
    type(string_t) :: type_name
    integer, parameter :: comm = MPI_COMM_WORLD

    if( musica_mpi_rank( comm ) == 0 ) then
      my_profile => profile_from_host_t( "foo", "bars", 3 )
      type_name = profile_type_name( my_profile )
      pack_size = type_name%pack_size( comm ) + my_profile%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call type_name%mpi_pack( buffer, pos, comm )
      call my_profile%mpi_pack( buffer, pos, comm )
      call assert( 102131176, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size, comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer, comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      call type_name%mpi_unpack( buffer, pos, comm )
      my_profile => profile_allocate( type_name )
      call my_profile%mpi_unpack( buffer, pos, comm )
      call assert( 997319163, pos <= pack_size )
    end if
    deallocate( buffer )

    call assert( 204373202, my_profile%handle_ == "foo" )
    call assert( 941386138, my_profile%units( ) == "bars" )
    call assert( 995865924, size( my_profile%mid_val_    ) == 3 )
    call assert( 992507151, size( my_profile%edge_val_   ) == 4 )
    call assert( 199561190, size( my_profile%delta_val_  ) == 3 )
    call assert( 371623628, size( my_profile%layer_dens_ ) == 3 )

    select type( my_profile )
    class is( profile_from_host_t )
      my_updater = profile_updater_t( my_profile )
    class default
      call die( 870112487 )
    end select

    call my_updater%update( mid_point_values = (/ 1.0_dk, 12.3_dk, 32.4_dk /),&
                         edge_values = (/ 0.5_dk, 9.8_dk, 15.4_dk, 45.0_dk /),&
                         layer_densities = (/ 94.3_dk, 0.52_dk, -12.3_dk /) )
    call assert( 792121290, my_profile%mid_val_(1) ==  1.0_dk )
    call assert( 339489137, my_profile%mid_val_(2) == 12.3_dk )
    call assert( 786856983, my_profile%mid_val_(3) == 32.4_dk )
    call assert( 399233230, my_profile%edge_val_(1) ==  0.5_dk )
    call assert( 846601076, my_profile%edge_val_(2) ==  9.8_dk )
    call assert( 393968923, my_profile%edge_val_(3) == 15.4_dk )
    call assert( 841336769, my_profile%edge_val_(4) == 45.0_dk )
    call assert( 106229367,                                                   &
                 almost_equal( my_profile%delta_val_(1),  9.8_dk -  0.5_dk ) )
    call assert( 901080862,                                                   &
                 almost_equal( my_profile%delta_val_(2), 15.4_dk -  9.8_dk ) )
    call assert( 448448709,                                                   &
                 almost_equal( my_profile%delta_val_(3), 45.0_dk - 15.4_dk ) )
    call assert( 613341306, my_profile%layer_dens_(1) == 94.3_dk )
    call assert( 725659651, my_profile%layer_dens_(2) == 0.52_dk )
    call assert( 338035898, my_profile%layer_dens_(3) == -12.3_dk )

    deallocate( my_profile )

  end subroutine test_profile_from_host_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_profile_from_host
