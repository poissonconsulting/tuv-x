! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_util
  ! Utility functions used in TUV-x

  use musica_constants,                only : musica_dk
  use musica_assert,                   only : die_msg

  implicit none

  private
  public :: inter2, inter4, addpnt

  real(musica_dk), parameter :: rZERO = 0.0_musica_dk

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inter2( xto, yto, xfrom, yfrom, ierr, debug )
    !  Map input data given on single, discrete points onto a set of target
    !  bins.
    !  The original input data are given on single, discrete points of an
    !  arbitrary grid and are being linearly interpolated onto a specified set
    !  of target bins.  In general, this is the case for most of the weighting
    !  functions (action spectra, molecular cross section, and quantum yield
    !  data), which have to be matched onto the specified wavelength intervals.
    !  The average value in each target bin is found by averaging the trapezoi-
    !  dal area underneath the input data curve (constructed by linearly connec-
    !  ting the discrete input values).
    !  Some caution should be used near the endpoints of the grids.  If the
    !  input data set does not span the range of the target grid, an error
    !  message is printed and the execution is stopped, as extrapolation of the
    !  data is not permitted.
    !  If the input data does not encompass the target grid, use ADDPNT to
    !  expand the input array.
    !
    ! \todo this function could use a more descriptive name

    real(musica_dk),              intent(in)  :: xto(:)   ! Target grid
    real(musica_dk),              intent(out) :: yto(:)   ! Regridded data
    real(musica_dk),              intent(in)  :: xfrom(:) ! Source grid
    real(musica_dk),              intent(in)  :: yfrom(:) ! Source data
    integer,                      intent(out) :: ierr     ! Error code (0 = success; >0 = error)
    logical, optional, intent(in)  :: debug    ! \todo needs description

    character(len=*), parameter :: Iam = 'Interpolation scheme 2'
    integer :: nto, nfrom
    integer :: ntom1, nfromm1
    integer :: i, k
    real(musica_dk) :: area, xtol, xtou
    real(musica_dk) :: slope
    real(musica_dk) :: a1, a2, b1, b2
    logical :: debugging

    if( present(debug) ) then
      debugging = debug
    else
      debugging = .false.
    endif

    ierr   = 0
    nfrom  = size( xfrom )
    nto    = size( xto )
    nfromm1 = nfrom - 1
    ntom1   = nto - 1

    !  check data grid for monotonicity
    if( any( xfrom( 2 : nfrom ) <= xfrom( 1 : nfromm1 ) ) ) then
      call die_msg( 860490498, Iam//'data grid not monotonically increasing' )
    endif

    !  do model grid x values lie completley inside data grid x values?
    if ( ( xfrom(1) > xto(1) ) .or. ( xfrom( nfrom ) < xto( nto ) ) ) then
      call die_msg( 127288630, Iam//'Data do not span grid; Use ADDPNT to '// &
                               'expand data and re-run.' )
    endif

    !  find the integral of each grid interval and use this to
    !  calculate the average y value for the interval
    !  xtol and xtou are the lower and upper limits of the to grid interval
    to_interval_loop: do i = 1, ntom1
      xtol = xto( i )
      xtou = xto( i + 1 )
      a1 = max( xfrom(1), xtol )
      a2 = min( xfrom( nfrom ), xtou )
      if( a2 > a1 ) then
        area = rZERO
        from_interval_loop: do k = 1, nfromm1
          a1 = max( xfrom( k ), xtol )
          a2 = min( xfrom( k + 1 ), xtou )
          if( a2 > a1 ) then
            slope = ( yfrom( k + 1 ) - yfrom( k ) )                           &
                    / ( xfrom( k + 1 ) - xfrom( k ) )
            b1 = yfrom( k ) + slope * ( a1 - xfrom( k ) )
            b2 = yfrom( k ) + slope * ( a2 - xfrom( k ) )
            area = area + .5_musica_dk * ( a2 - a1 ) * ( b2 + b1 )
          endif
        enddo from_interval_loop
        yto( i ) = area / ( xtou - xtol )
      else
        yto( i ) = rZERO
      endif
    enddo to_interval_loop

  end subroutine inter2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine inter4( xto, yto, xfrom, yfrom, fold_in )
    ! Map input data given on a set of bins onto a different set of target
    ! bins.
    ! The input data are given on a set of bins (representing the integral
    ! of the input quantity over the range of each bin) and are being matched
    ! onto another set of bins (target grid).  A typical example would be an
    ! input data set spcifying the extra-terrestrial flux on wavelength inter-
    ! vals, that has to be matched onto the working wavelength grid.
    ! The resulting area in a given bin of the target grid is calculated by
    ! simply adding all fractional areas of the input data that cover that
    ! particular target bin.
    ! Some caution should be used near the endpoints of the grids.  If the
    ! input data do not span the full range of the target grid, the area in
    ! the "missing" bins will be assumed to be zero.  If the input data extend
    ! beyond the upper limit of the target grid, the user has the option to
    ! integrate the "overhang" data and fold the remaining area back into the
    ! last target bin.  Using this option is recommended when re-gridding
    ! vertical profiles that directly affect the total optical depth of the
    ! model atmosphere.
    !
    ! \todo this function could use a more descriptive name

    real(musica_dk), intent(in)  :: xto(:)   ! Target grid
    real(musica_dk), intent(out) :: yto(:)   ! Regridded data
    real(musica_dk), intent(in)  :: xfrom(:) ! Source grid
    real(musica_dk), intent(in)  :: yfrom(:) ! Source data
    integer,         intent(in)  :: fold_in  ! 0 -> No folding of "overhang" data
    ! 1 -> Integrate "overhang" data and fold back into last target bin

    character(len=*), parameter :: Iam = 'Interpolation scheme 4'
    integer :: nfrom, nto
    integer :: nfromm1, ntom1
    real(musica_dk) :: a1, a2, sum
    real(musica_dk) :: tail
    integer :: jstart, i, j, k

    ! check whether flag given is legal
    if ( fold_in /= 0 .and. fold_in /= 1 ) then
      call die_msg( 306355737, Iam//'Value for fold_in invalid. '//         &
                               'Must be 0 or 1' )
    endif
    nfrom   = size( xfrom )
    nfromm1 = nfrom - 1
    nto     = size( xto )
    ntom1   = nto - 1

    ! do interpolation
    jstart = 1
    do i = 1, ntom1
      yto( i ) = rZERO
      sum      = rZERO
      j = jstart
      if( j <= nfromm1 ) then
        do while( xfrom( j + 1 ) < xto( i ) )
          jstart = j
          j = j + 1
          if( j > nfromm1 ) then
            exit
          endif
        enddo

        do while( ( xfrom( j ) <= xto( i + 1 ) ) .and. ( j <= nfromm1 ) )
          a1 = max( xfrom( j ), xto( i ) )
          a2 = min( xfrom( j + 1 ), xto( i + 1 ) )
          sum = sum + yfrom( j ) * ( a2 - a1 )
          j = j + 1
        enddo
        yto( i ) = sum / ( xto( i + 1 ) - xto( i ) )
      endif
    enddo

    ! if requested, integrate data "overhang" and fold back into last bin
    if( fold_in == 1 ) then
      j = j - 1
      a1 = xto( nto )       ! upper limit of last interpolated bin
      a2 = xfrom( j + 1 )   ! upper limit of last input bin considered

      ! do fold_ing only if grids don't match up and there is more input
      if( a2 > a1 .or. j + 1 < nfrom ) then
        tail = yfrom( j ) * ( a2 - a1 ) / ( xfrom( j + 1 ) - xfrom( j ) )
        do k = j + 1, nfromm1
          tail = tail + yfrom( k ) * ( xfrom( k + 1 ) - xfrom( k ) )
        enddo
        yto( ntom1 ) = yto( ntom1 ) + tail
      endif
    endif

  end subroutine inter4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine addpnt ( x, y, xnew, ynew )
    ! Add a point <xnew,ynew> to a set of data pairs <x,y>.  x must be in
    ! ascending order

    real(musica_dk), allocatable, intent(inout) :: x(:) ! Grid
    real(musica_dk), allocatable, intent(inout) :: y(:) ! Data
    real(musica_dk),              intent(in)    :: xnew ! Grid point to add
    real(musica_dk),              intent(in)    :: ynew ! Data point to add

    character(len=*), parameter :: Iam = 'Add point to gridded data'
    integer :: n
    integer :: insertNdx
    real(musica_dk), allocatable    :: wrk(:)
    logical            :: found

    n = size( x )

    !  check data grid for monotonicity
    if( any( x( 2 : n ) <= x( 1 : n - 1 ) ) ) then
      write(*,*) Iam, 'grid not monotonically increasing'
      stop 3
    endif

    !  does xnew == any x value?
    if( any( x(:) == xnew ) ) then
      write(*,*) Iam, 'xnew exactly matches a grid x value'
      stop 3
    endif

    ! find the index at which xnew needs to be inserted into x
    found = .true.
    if( xnew < x(1) ) then
      insertNdx = 1
    else if( xnew > x( n ) ) then
      insertNdx = n
    else
      found = .false.
      do insertNdx = 2, n
        if ( x( insertNdx ) > xnew ) then
          found = .true.
          exit
        endif
      enddo
    endif
    if( .not. found ) then
      write(*,*) Iam, 'something really wrong; all stop'
      stop 3
    endif

    ! increment x,y arrays, then insert xnew,ynew
    if( insertNdx == 1 ) then
      x = [ xnew, x ]
      y = [ ynew, y ]
    elseif( insertNdx == n ) then
      x = [ x, xnew ]
      y = [ y, ynew ]
    else
      wrk = [ x( : insertNdx - 1 ), xnew ]
      x   = [ wrk, x( insertNdx : ) ]
      wrk = [ y( : insertNdx - 1 ), ynew ]
      y   = [ wrk, y( insertNdx : ) ]
    endif

  end subroutine addpnt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_util
