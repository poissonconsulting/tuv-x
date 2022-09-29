! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!

module tuvx_diagnostic_util
  ! Diagnostic utilities

  use musica_constants, only : dk => musica_dk
  use musica_string,    only : string_t

  implicit none

  private
  public :: diagout, prepare_diagnostic_output

  interface diagout
    module procedure :: diagnostic_1d
    module procedure :: diagnostic_1d_dk
    module procedure :: diagnostic_2d
    module procedure :: diagnostic_2d_dk
    module procedure :: diagnostic_array_string_t
  end interface diagout

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function output_enabled( enable_output )
    ! Returns true if the output policy enables diagnostic output

    logical, optional :: enable_output ! Enables diagnostic output

    if( .not. present( enable_output ) ) then
      output_enabled = .true.
      return
    end if

    output_enabled = enable_output
  end function output_enabled

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine prepare_diagnostic_output( enable_output )
    ! Creates the folder to hold diagnostic output if it doesn't exist
    logical, optional :: enable_output ! Enables diagnostic output

    if (output_enabled( enable_output )) then
      call execute_command_line( "mkdir -p output"  )
    endif

  end subroutine prepare_diagnostic_output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine diagnostic_1d( filename, variable, enable_output )
    ! Output 1D float diagnostics to a specified file

    character(len=*), intent(in) :: filename      ! File path to output to
    real, intent(in)             :: variable(:)   ! Diagnostics to output
    logical, optional            :: enable_output ! Enables diagnostic output

    integer :: ios

    if (.not. output_enabled( enable_output )) return

    open( unit = 44, file = 'output/' // filename, form = 'unformatted',      &
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

  subroutine diagnostic_1d_dk( filename, variable, enable_output )
    ! Output 1D double diagnostics to a specified file

    character(len=*), intent(in) :: filename    ! File path to output to
    real(dk), intent(in)         :: variable(:) ! Diagnostics to output
    logical, optional            :: enable_output ! Enables diagnostic output

    character(len=*), parameter  :: Iam = 'diagnostic_1d_dk: '

    integer :: ios
    character(len=256) :: iomsg

    if (.not. output_enabled( enable_output )) return

    open( unit = 44, file = 'output/' // filename, form = 'unformatted',      &
      iostat = ios)
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

  subroutine diagnostic_2d( filename, variable, enable_output )
    ! Ouptut 2D float diagnostics to a specified file

    character(len=*), intent(in) :: filename      ! File path to output to
    real, intent(in)             :: variable(:,:) ! Diagnostics to output
    logical, optional            :: enable_output ! Enables diagnostic output

    integer :: ios

    if (.not. output_enabled( enable_output )) return

    open( unit = 44, file = 'output/' // filename, form = 'unformatted',      &
      iostat = ios)
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

  subroutine diagnostic_2d_dk( filename, variable, enable_output )
    ! Output 2D double diagnostics to a specified file

    character(len=*), intent(in) :: filename      ! File path to output to
    real(dk), intent(in)         :: variable(:,:) ! Diagnostics to output
    logical, optional            :: enable_output ! Enables diagnostic output

    integer :: ios
    character(len=512) :: iomsg

    if (.not. output_enabled( enable_output )) return

    open( unit = 44, file = 'output/' // filename, form = 'unformatted',      &
      iostat = ios)
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

  subroutine diagnostic_array_string_t( filename, variable, enable_output )
    ! Output 2D double diagnostics to a specified file

    character(len=*), intent(in)              :: filename      ! File path to output to
    type(string_t), allocatable, intent(in)   :: variable(:) ! Diagnostics to output
    logical, optional                         :: enable_output ! Enables diagnostic output

    integer :: ios, idx
    character(len=512) :: iomsg

    if (.not. output_enabled( enable_output )) return

    open( unit = 44, file = 'output/' // filename, iostat = ios)
    if( ios /= 0 ) then
      write(*,*) 'diagnostic_array_string_t: failed to open ', filename, '; error = ', &
                 ios
      write(*,*) trim( iomsg )
      stop 3
    endif
    do idx = 1, size(variable)
      write(44,*) variable( idx )%to_char( )
    enddo
    if( ios /= 0 ) then
      write(*,*) 'diagnostic_array_string_t: failed to write ', filename, '; error = ',&
                  ios
      write(*,*) trim( iomsg )
      stop 3
    endif
    close( unit = 44 )

   end subroutine diagnostic_array_string_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_diagnostic_util
