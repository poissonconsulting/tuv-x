! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_diagnostic_util
  ! Diagnostic utilities

  use musica_constants, only : dk => musica_dk

  implicit none

  interface diagout
    module procedure :: diagnostic_1d
    module procedure :: diagnostic_1d_dk
    module procedure :: diagnostic_2d
    module procedure :: diagnostic_2d_dk
  end interface diagout

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine prepare_diagnostic_output( )
    ! Creates the folder to hold diagnostic output if it doesn't exist

    call system( "mkdir -p output" )

  end subroutine prepare_diagnostic_output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine diagnostic_1d( filename, variable )
    ! Output 1D float diagnostics to a specified file

    character(len=*), intent(in) :: filename    ! File path to output to
    real, intent(in)             :: variable(:) ! Diagnostics to output

    integer :: ios

    open( unit = 44, file = 'output/'//filename, form = 'unformatted',        &
          iostat = ios)
    if( ios /= 0 ) then
      write(*,*) 'diagnostic_1d: failed to open ', filename, '; error = ', ios
      stop 3
    endif
    write( unit = 44, iostat = ios ) variable
    if( ios /= 0 ) then
      write(*,*) 'diagnostic_1d: failed to write ', filename, '; error = ', ios
      stop 3
    endif
    close( unit = 44 )

  end subroutine diagnostic_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine diagnostic_1d_dk( filename, variable )
    ! Output 1D double diagnostics to a specified file

    character(len=*), intent(in) :: filename    ! File path to output to
    real(dk), intent(in)         :: variable(:) ! Diagnostics to output

    character(len=*), parameter  :: Iam = 'diagnostic_1d_dk: '

    integer :: ios
    character(len=256) :: iomsg

    open( unit = 44, file = 'output/'//filename, form = 'unformatted',        &
          iostat = ios, iomsg = iomsg )
    if( ios /= 0 ) then
      write(*,*) 'diagnostic_1d: failed to open ', filename, '; error = ', ios
      write(*,*) trim( iomsg )
      stop 3
    endif
    write( unit = 44, iostat = ios, iomsg = iomsg ) variable
    if( ios /= 0 ) then
      write(*,*) 'diagnostic_1d: failed to write ', filename, '; error = ', ios
      write(*,*) trim( iomsg )
      stop 3
    endif
    close(unit=44)

  end subroutine diagnostic_1d_dk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine diagnostic_2d( filename, variable )
    ! Ouptut 2D float diagnostics to a specified file

    character(len=*), intent(in) :: filename      ! File path to output to
    real, intent(in)             :: variable(:,:) ! Diagnostics to output

    integer :: ios

    open( unit = 44, file = 'output/'//filename, form = 'unformatted',        &
          iostat = ios )
    if( ios /= 0 ) then
      write(*,*) 'diagnostic_2d: failed to open ', filename, '; error = ', ios
      stop 3
    endif
    write( unit = 44, iostat = ios ) variable
    if( ios /= 0 ) then
      write(*,*) 'diagnostic_2d: failed to write ', filename, '; error = ', ios
      stop 3
    endif
    close( unit = 44 )

  end subroutine diagnostic_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine diagnostic_2d_dk( filename, variable )
    ! Output 2D double diagnostics to a specified file

    character(len=*), intent(in) :: filename      ! File path to output to
    real(dk), intent(in)         :: variable(:,:) ! Diagnostics to output

    integer :: ios
    character(len=512) :: iomsg

    open( unit = 44, file = 'output/'//filename, form = 'unformatted',        &
          iostat = ios, iomsg = iomsg )
    if( ios /= 0 ) then
      write(*,*) 'diagnostic_2d_dk: failed to open ', filename, '; error = ', &
                 ios
      write(*,*) trim( iomsg )
      stop 3
    endif
    write( unit=44, iostat = ios, iomsg = iomsg ) variable
    if( ios /= 0 ) then
      write(*,*) 'diagnostic_2d_dk: failed to write ', filename, '; error = ',&
                  ios
      write(*,*) trim( iomsg )
      stop 3
    endif
    close( unit = 44 )

   end subroutine diagnostic_2d_dk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_diagnostic_util
