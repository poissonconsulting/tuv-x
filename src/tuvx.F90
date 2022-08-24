! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
program tuvx
  ! Driver for stand-alone TUV-x

  use musica_string,                   only : string_t
  use tuvx_core,                       only : core_t
  use tuvx_diagnostic_util,            only : prepare_diagnostic_output

  implicit none

  class(core_t), pointer :: core

  ! Command-line options
  character(len=256) :: argument
  type(string_t)     :: configFileSpec

  ! Get the model configuration file and options from the command line
  if( command_argument_count() /= 1 ) then
    call fail_run( )
  endif
  call get_command_argument( 1, argument )

  configFileSpec = argument

  ! set up diagnostic output
  call prepare_diagnostic_output( )

  ! instatiate and initialize photolysis core object
  core => core_t( configFileSpec )

  ! run photolysis
  call core%run()

  deallocate( core )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fail_run( )
    ! Fail run and print usage info

    write(*,*) "Usage: ./photolysis configuration_file.json"
    stop 3

  end subroutine fail_run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program tuvx
