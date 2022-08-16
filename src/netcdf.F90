! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_netcdf
  ! NetCDF I/O class

  use nc4fortran,                      only : netcdf_file
  use musica_constants,                only : musica_dk

  implicit none

  private
  public :: netcdf_t, clean_string

  type netcdf_t
    ! NetCDF I/O
    real(musica_dk), allocatable :: wavelength(:)
    real(musica_dk), allocatable :: temperature(:)
    real(musica_dk), allocatable :: parameters(:,:)
  contains
    procedure :: read_netcdf_file => run
    final     :: finalize
  end type netcdf_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine run( this, file_path, variable_name )
    ! Reads data from a NetCDF file

    class(netcdf_t), intent(inout) :: this
    character(len=*), intent(in)   :: file_path
    character(len=*), intent(in)   :: variable_name

    integer, parameter :: noErr = 0
    character(len=*), parameter :: Iam = "read_netcdf_file: "

    integer :: stat, nLambda
    integer, allocatable :: dims(:)
    character(:), allocatable :: varName
    type(netcdf_file) :: ncObj

    !  open the netcdf file
    call ncObj%initialize( file_path, ierr = stat, action = 'r' )
    if( stat /= noErr ) then
      write(*,*) Iam, 'retcode from initialize = ', stat
      stop 3
    endif

    !  parameter array must be in netcdf file, read it
    varName = trim( variable_name ) // 'parameters'
    if( ncObj%exist( varName ) ) then
      call ncObj%shape( varName, dims )
      nLambda = dims(1)
      ! read input parameters array
      allocate( this%parameters( dims(1), dims(2) ) )
      call ncObj%read( varName, this%parameters )
    else
      write(*,*) Iam, ' variable ', trim( varName ), ' not in netcdf file'
      stop 3
    endif

    !  if it exists, read wavelength array
    if( ncObj%exist( 'wavelength' ) ) then
      call ncObj%shape( 'wavelength', dims )
      if( dims(1) /= nLambda ) then
        write(*,*) Iam, ' wavelength, parameters array size mismatch'
        stop 3
      endif
      allocate( this%wavelength( dims(1) ) )
      call ncObj%read( 'wavelength', this%wavelength )
    endif

    !  if it exists, read temperature array
    if( ncObj%exist( 'temperature' ) ) then
      call ncObj%shape( 'temperature', dims )
      allocate( this%temperature( dims(1) ) )
      call ncObj%read( 'temperature', this%temperature )
    endif

    !  close the netcdf file
    call ncObj%finalize( )

  end subroutine run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Cleans up memory of the NetCDF class

    type(netcdf_t), intent(inout) :: this

    if( allocated( this%wavelength ) ) then
      deallocate( this%wavelength )
    endif
    if( allocated( this%temperature ) ) then
      deallocate( this%temperature )
    endif
    if( allocated( this%parameters ) ) then
      deallocate( this%parameters )
    endif

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(string_t) function clean_string( from )

    use musica_string,                 only : string_t

    type(string_t), intent(in) :: from

    clean_string = from%replace( "/", "_" )

  end function clean_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_netcdf
