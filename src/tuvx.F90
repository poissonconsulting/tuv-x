! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
program tuvx
  ! Driver for stand-alone TUV-x and regression tests
  !
  ! If MPI support is compiled in, this driver constructs the core
  ! on the primary process, passes it to the other processes, and
  ! performs calculations on process 1. This allows MPI functions to
  ! be tested as part of the regression tests.

  use musica_assert,                   only : assert
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
    call core%run()
  end if

  deallocate( buffer )
  deallocate( core )

  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fail_run( )
    ! Fail run and print usage info

    write(*,*) "Usage: ./photolysis configuration_file.json"
    stop 3

  end subroutine fail_run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program tuvx
