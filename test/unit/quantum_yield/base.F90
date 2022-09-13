! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
program test_quantum_yield

  use musica_mpi,                      only : musica_mpi_init,                &
                                              musica_mpi_finalize
  use tuvx_quantum_yield

  implicit none

  call musica_mpi_init( )
  call test_quantum_yield_t( )
  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_quantum_yield_t( )
    ! Test functionality of the :f:type:`~tuvx_quantum_yield/quantum_yield_t`
    ! class.
    !
    ! This test only checks the MPI functions currently. Additional test will
    ! be added as part of issue #177

    use musica_assert,                 only : assert
    use musica_constants,              only : dk => musica_dk
    use musica_mpi
    use musica_string,                 only : string_t
    use tuvx_quantum_yield_factory,    only : quantum_yield_type_name,        &
                                              quantum_yield_allocate
    use tuvx_test_utils,               only : check_values

    class(quantum_yield_t), pointer :: quantum_yield
    character, allocatable :: buffer(:)
    type(string_t) :: type_name
    integer :: pos, pack_size

    real(dk) :: temperature1(3), temperature2(4)
    real(dk) :: array1(3,2), array2(2,1)

    temperature1(:) = (/  12.5_dk,  16.3_dk, -290.4_dk /)
    array1(:,1)     = (/  42.3_dk, 132.4_dk,   13.4_dk /)
    array1(:,2)     = (/ 132.4_dk,  0.43_dk,   2.34_dk /)

    temperature2(:) = (/ -123.4_dk, 41.2_dk, 0.053_dk, 1.2e-7_dk /)
    array2(:,1)     = (/ 12.34_dk, -142.3_dk /)

    if( musica_mpi_rank( ) == 0 ) then
      allocate( quantum_yield )
      allocate( quantum_yield%quantum_yield_parms(2) )
      quantum_yield%quantum_yield_parms(1)%temperature = temperature1
      quantum_yield%quantum_yield_parms(1)%array = array1
      quantum_yield%quantum_yield_parms(2)%temperature = temperature2
      quantum_yield%quantum_yield_parms(2)%array = array2
      type_name = quantum_yield_type_name( quantum_yield )
      pack_size = type_name%pack_size( ) + quantum_yield%pack_size( )
      allocate( buffer( pack_size ) )
      pos = 0
      call type_name%mpi_pack(     buffer, pos )
      call quantum_yield%mpi_pack( buffer, pos )
      call assert( 209765802, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size )
    if( musica_mpi_rank( ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer )

    if( musica_mpi_rank( ) .ne. 0 ) then
      pos = 0
      call type_name%mpi_unpack( buffer, pos )
      quantum_yield => quantum_yield_allocate( type_name )
      call quantum_yield%mpi_unpack( buffer, pos )
      call assert( 264697883, pos <= pack_size )
    end if
    deallocate( buffer )

    call assert( 962406487, associated( quantum_yield ) )
    call assert( 964312021, allocated( quantum_yield%quantum_yield_parms ) )
    call assert( 401267057, size( quantum_yield%quantum_yield_parms ) == 2 )
    associate( params => quantum_yield%quantum_yield_parms( 1 ) )
      call assert( 784078798, allocated( params%temperature ) )
      call assert( 891132836, allocated( params%array ) )
      call check_values( params%temperature, temperature1, 1.0e-6_dk )
      call check_values( params%array,       array1,       1.0e-6_dk )
    end associate
    associate( params => quantum_yield%quantum_yield_parms( 2 ) )
      call assert( 784078798, allocated( params%temperature ) )
      call assert( 891132836, allocated( params%array ) )
      call check_values( params%temperature, temperature2, 1.0e-6_dk )
      call check_values( params%array,       array2,       1.0e-6_dk )
    end associate

    deallocate( quantum_yield )

  end subroutine test_quantum_yield_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_quantum_yield
