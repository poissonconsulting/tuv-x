! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_interpolate
 ! Classes that provide interpolation functions

  use musica_constants,                only : dk => musica_dk

  implicit none

  private
  public :: interpolator_t, interp1_t, interp2_t, interp3_t, interp4_t

  real(dk), parameter :: rZERO = 0.0_dk

  type, abstract :: interpolator_t
    ! General interpolator interface
    character(len=10) :: handle_
  contains
    procedure(interpolate), deferred :: interpolate
  end type interpolator_t


  type, extends(interpolator_t) :: interp1_t
    ! Type 1 interpolator
  contains
    procedure :: interpolate => inter1
  end type interp1_t

  type, extends(interpolator_t) :: interp2_t
    ! Type 2 interpolator
  contains
    procedure :: interpolate => inter2
  end type interp2_t

  type, extends(interpolator_t) :: interp3_t
    ! Type 3 interpolator
  contains
    procedure :: interpolate => inter3
  end type interp3_t

  type, extends(interpolator_t) :: interp4_t
    ! Type 4 interpolator
  contains
    procedure :: interpolate => inter4
  end type interp4_t

interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function interpolate( this, x_target, x_source, y_source, fold_in )         &
      result( y_target )
    ! Interpolates data in y_source on the x_source axis, to y_target on the
    ! x_target axis

    use musica_constants, only : dk => musica_dk

    import interpolator_t

    class(interpolator_t), intent(inout) :: this        ! Interpolator
    real(dk),              intent(in)    :: x_target(:) ! Target axis
    real(dk),              intent(in)    :: x_source(:) ! Source axis
    real(dk),              intent(in)    :: y_source(:) ! Source data
    integer, optional,     intent(in)    :: fold_in     ! Switch for folding option of "overhang" data

    real(dk), allocatable :: y_target(:) ! Target data

  end function interpolate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function inter1(this, x_target, x_source, y_source, fold_in)                &
      result( y_target )
    !  Map input data given on single, discrete points, onto a discrete target
    !  grid.
    !  The original input data are given on single, discrete points of an
    !  arbitrary grid and are being linearly interpolated onto a specified
    !  discrete target grid.  A typical example would be the re-gridding of a
    !  given data set for the vertical temperature profile to match the speci-
    !  fied altitude grid.
    !  Some caution should be used near the end points of the grids.  If the
    !  input data set does not span the range of the target grid, the remaining
    !  points will be set to zero, as extrapolation is not permitted.
    !  If the input data does not encompass the target grid, use ADDPNT to
    !  expand the input array.

    class(interp1_t),  intent(inout) :: this        ! Interpolator
    real(dk),          intent(in)    :: x_target(:) ! Target axis
    real(dk),          intent(in)    :: x_source(:) ! Source axis
    real(dk),          intent(in)    :: y_source(:) ! Source data
    integer, optional, intent(in)    :: fold_in     ! Not used

    real(dk), allocatable :: y_target(:)

    integer :: n
    integer :: jsave, i, j
    real(dk)    :: slope

    allocate( y_target( size( x_target ) ) )

    n = size( x_source )
    jsave = 1
    y_target = rZERO
    do i = 1, size( x_target )
      j = jsave
      do
        if( x_target( i ) < x_source( j ) .or.                                &
            x_target( i ) >= x_source( j + 1 ) ) then
          j = j+1
          if ( j >= n ) then
            exit
          endif
        else
          slope = ( y_source( j + 1 ) - y_source( j ) )                       &
                  / (x_source( j + 1 ) - x_source( j ) )
          y_target( i ) = y_source( j )                                       &
                          + slope * ( x_target( i ) - x_source( j ) )
          jsave = j
          exit
        endif
      enddo
    enddo

  end function inter1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function inter2(this, x_target, x_source, y_source, fold_in)                &
      result( y_target )
    !  Map input data given on single, discrete points onto a set of target
    !  bins.
    !  The original input data are given on single, discrete points of an
    !  arbitrary grid and are being linearly interpolated onto a specified set
    !  of target bins.  In general, this is the case for most of the weighting
    !  functions (action spectra, molecular cross section, and quantum yield
    !  data), which have to be matched onto the specified wavelength intervals.
    !  The average value in each target bin is found by averaging the trapezoi-
    !  dal area underneath the input data curve (constructed by linearly
    !  connectting the discrete input values).
    !  Some caution should be used near the endpoints of the grids.  If the
    !  input data set does not span the range of the target grid, an error
    !  message is printed and the execution is stopped, as extrapolation of the
    !  data is not permitted.
    !  If the input data does not encompass the target grid, use ADDPNT to
    !  expand the input array.

    use musica_assert, only : die_msg

    class(interp2_t),  intent(inout) :: this        ! Interpolator
    real(dk),          intent(in)    :: x_target(:) ! Target axis
    real(dk),          intent(in)    :: x_source(:) ! Source axis
    real(dk),          intent(in)    :: y_source(:) ! Source data
    integer, optional, intent(in)    :: fold_in     ! Not used

    real(dk), allocatable   :: y_target(:)

    integer :: ng, n, ngintv
    integer :: i, k, jstart
    real(dk) :: area, xgl, xgu
    real(dk) :: darea, slope
    real(dk) :: a1, a2, b1, b2

    n  = size( x_source )
    ng = size( x_target )

    allocate( y_target( ng - 1 ) )

    ! test for correct ordering of data, by increasing value of x
    if( any( x_source( 1 : n - 1 ) >= x_source( 2 : n ) ) ) then
      call die_msg( 996313430,'src grid must be monotonically increasing' )
    endif

    if( any( x_target( 1 : ng - 1 ) >= x_target( 2 : ng ) ) ) then
      call die_msg( 208631776,'target grid must be monotonically increasing' )
    endif

    ! check for xg-values outside the x-range
    if( ( x_source(1) > x_target(1) ) .or.                                    &
        ( x_source( n ) < x_target( ng ) ) ) then
      call die_msg( 320950121,'src and target grid do not overlap' )
    endif

    y_target = rZERO
    !  find the integral of each grid interval and use this to
    !  calculate the average y value for the interval
    !  xgl and xgu are the lower and upper limits of the grid interval

    jstart = 1
    ngintv = ng - 1
    do i = 1, ngintv
      ! initalize:
      area = rZERO
      xgl = x_target( i )
      xgu = x_target( i + 1 )
      ! discard data before the first grid interval and after the
      ! last grid interval
      ! for internal grid intervals, start calculating area by interpolating
      ! between the last point which lies in the previous interval and the
      ! first point inside the current interval
      k = jstart
      if ( k < n ) then
        ! if both points are before the first grid, go to the next point
        do
          if( x_source( k + 1 ) <= xgl ) then
            jstart = k - 1
            k = k + 1
            if( k >= n ) then
              exit
            endif
          else
            exit
          endif
        enddo
        ! if the last point is beyond the end of the grid, complete and go
        ! to the next grid
        do while( k < n .and. x_source( k ) < xgu )
          jstart = k - 1
          ! compute x-coordinates of increment
          a1 = max( x_source( k ), xgl )
          a2 = min( x_source( k + 1 ), xgu )
          ! if points coincide, contribution is zero
          if( x_source( k + 1 ) == x_source( k ) ) then
            darea = rZERO
          else
            slope = ( y_source( k + 1 ) - y_source( k ) )                     &
                    / ( x_source( k + 1 ) - x_source( k ) )
            b1 = y_source( k ) + slope * ( a1 - x_source( k ) )
            b2 = y_source( k ) + slope * ( a2 - x_source( k ) )
            darea = .5_dk * ( a2 - a1 ) * ( b2 + b1 )
          endif
          ! find the area under the trapezoid from a1 to a2
          area = area + darea
          ! go to next point
          k = k + 1
        enddo
      endif
      ! calculate the average y after summing the areas in the interval
      y_target( i ) = area / ( xgu - xgl )
    enddo

  end function inter2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function inter3(this, x_target, x_source, y_source, fold_in)                &
      result( y_target )
    !  Map input data given on a set of bins onto a different set of target
    !  bins.
    !  The input data are given on a set of bins (representing the integral
    !  of the input quantity over the range of each bin) and are being matched
    !  onto another set of bins (target grid).  A typical example would be an
    !  input data set spcifying the extra-terrestrial flux on wavelength inter-
    !  vals, that has to be matched onto the working wavelength grid.
    !  The resulting area in a given bin of the target grid is calculated by
    !  simply adding all fractional areas of the input data that cover that
    !  particular target bin.
    !  Some caution should be used near the endpoints of the grids.  If the
    !  input data do not span the full range of the target grid, the area in
    !  the "missing" bins will be assumed to be zero.  If the input data extend
    !  beyond the upper limit of the target grid, the user has the option to
    !  integrate the "overhang" data and fold the remaining area back into the
    !  last target bin.  Using this option is recommended when re-gridding
    !  vertical profiles that directly affect the total optical depth of the
    !  model atmosphere.
    !
    !  fold_in - Switch for folding option of "overhang" data
    !            fold_in = 0 -> No folding of "overhang" data
    !            fold_in = 1 -> Integerate "overhang" data and fold back into
    !                           last target bin

    use musica_assert, only : die_msg

    class(interp3_t),  intent(inout) :: this        ! Interpolator
    real(dk),          intent(in)    :: x_target(:) ! Target axis
    real(dk),          intent(in)    :: x_source(:) ! Source axis
    real(dk),          intent(in)    :: y_source(:) ! Source data
    integer, optional, intent(in)    :: fold_in     ! Switch for folding option of "overhang" data

    real(dk), allocatable   :: y_target(:)

    integer :: nfrom, nto
    integer :: jstart, i, j, k, ntobins
    real(dk)    :: a1, a2, sum
    real(dk)    :: tail
    logical :: do_fold_in

    ! check whether flag given is legal
    if( present( fold_in ) ) then
      if( fold_in /= 0 .and. fold_in /= 1 ) then
        call die_msg( 433268466, 'fold_in must be 0 or 1' )
      endif
      do_fold_in = fold_in == 1
    else
      call die_msg( 828062060, 'fold_in argument not present' )
    endif

    allocate( y_target( size( x_target ) - 1 ) )

    associate( xto => x_target, yto => y_target,                              &
               xfrom => x_source, yfrom => y_source )

    nto   = size( xto )
    nfrom = size( xfrom )

    jstart = 1

    yto     = rZERO
    ntobins = nto - 1
    do  i = 1, ntobins
      sum = rZERO
      j = jstart

      if( j < nfrom ) then
        do
          if( xfrom( j + 1 ) < xto( i ) ) then
            jstart = j
            j = j+1
            if( j >= nfrom ) then
              exit
            endif
          else
            exit
          endif
        enddo

        do while( ( xfrom( j ) <= xto( i + 1 ) ) .and. ( j < nfrom ) )
          a1 = max( xfrom( j ), xto( i ) )
          a2 = min( xfrom( j + 1 ), xto( i + 1 ) )
          sum = sum + yfrom( j ) * ( a2 - a1 )                                &
                      / ( xfrom( j + 1 ) - xfrom( j ) )
          j = j + 1
        enddo
        yto( i ) = sum
      endif
    enddo

    ! integrate data "overhang" and fold back into last target bin
    if ( do_fold_in ) then
       j = j - 1
       a1 = xto( nto )           ! upper limit of last interpolated bin
       a2 = xfrom( j + 1 )    ! upper limit of last input bin considered
       ! do folding only if grids don't match up and there is more input
       if( a2 > a1 .or. j + 1 < nfrom ) then
          tail = yfrom( j ) * ( a2 - a1 ) / ( xfrom( j + 1 ) - xfrom( j ) )
          do k = j + 1, nfrom - 1
             tail = tail + yfrom( k ) * ( xfrom( k + 1 ) - xfrom( k ) )
          enddo
          yto( ntobins ) = yto( ntobins ) + tail
       endif
    endif

    end associate

  end function inter3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function inter4(this, x_target, x_source, y_source, fold_in)                &
      result( y_target )
    !  Map input data given on a set of bins onto a different set of target
    !  bins.
    !  The input data are given on a set of bins (representing the integral
    !  of the input quantity over the range of each bin) and are being matched
    !  onto another set of bins (target grid).  A typical example would be an
    !  input data set spcifying the extra-terrestrial flux on wavelength inter-
    !  vals, that has to be matched onto the working wavelength grid.
    !  The resulting area in a given bin of the target grid is calculated by
    !  simply adding all fractional areas of the input data that cover that
    !  particular target bin.
    !  Some caution should be used near the endpoints of the grids.  If the
    !  input data do not span the full range of the target grid, the area in
    !  the "missing" bins will be assumed to be zero.  If the input data extend
    !  beyond the upper limit of the target grid, the user has the option to
    !  integrate the "overhang" data and fold the remaining area back into the
    !  last target bin.  Using this option is recommended when re-gridding
    !  vertical profiles that directly affect the total optical depth of the
    !  model atmosphere.
    !
    !  fold_in - Switch for folding option of "overhang" data
    !            fold_in = 0 -> No folding of "overhang" data
    !            fold_in = 1 -> Integerate "overhang" data and fold back into
    !                           last target bin

    use musica_assert, only : die_msg

    class(interp4_t),  intent(inout) :: this        ! Interpolator
    real(dk),          intent(in)    :: x_target(:) ! Target axis
    real(dk),          intent(in)    :: x_source(:) ! Source axis
    real(dk),          intent(in)    :: y_source(:) ! Source data
    integer, optional, intent(in)    :: fold_in     ! Switch for folding option of "overhang" data

    real(dk), allocatable :: y_target(:)

    integer :: n, ng
    integer :: jstart, i, j, k
    real(dk)    :: a1, a2, sum
    real(dk)    :: tail
    logical :: do_fold_in

    ! check whether flag given is legal
    if( present( fold_in ) ) then
      if( fold_in /= 0 .and. fold_in /= 1 ) then
        call die_msg( 322855655, 'fold_in must be 0 or 1' )
      endif
      do_fold_in = fold_in == 1
    else
      call die_msg( 152698751, 'fold_in argument not present' )
    endif

    n  = size( x_source )
    ng = size( x_target )
    allocate( y_target( ng - 1 ) )
    jstart  = 1
    y_target = rZERO

    do i = 1, ng - 1
       sum = rZERO
       j = jstart
       if (j < n) then
         do
           if ( x_source( j + 1 ) < x_target( i ) ) then
             jstart = j
             j = j + 1
             if (j >= n) then
               exit
             endif
           else
             exit
           endif
         enddo

         do while( x_source( j ) <= x_target( i + 1 ) .and. j < n )
            a1 = max( x_source( j ), x_target( i ) )
            a2 = min( x_source( j + 1 ), x_target( i + 1 ) )
            sum = sum + y_source( j ) * ( a2 - a1 )
            j = j + 1
         enddo
         y_target( i ) = sum / ( x_target( i + 1 ) - x_target( i ) )
      endif
    enddo
    ! integrate data "overhang" and fold back into last bin
    if( do_fold_in ) then
       j = j - 1
       a1 = x_target( ng )       ! upper limit of last interpolated bin
       a2 = x_source( j + 1 ) ! upper limit of last input bin considered

       ! do folding only if grids don't match up and there is more input
       if ( ( a2 > a1 ) .or. ( j + 1 < n ) ) then
         tail = y_source( j ) * ( a2 - a1 )                                   &
                / ( x_source( j + 1 ) - x_source( j ) )
         do k = j + 1, n - 1
            tail = tail + y_source( k )                                       &
                          * ( x_source( k + 1 ) - x_source( k ) )
         enddo
         y_target( ng - 1 ) = y_target( ng - 1 ) + tail
       endif
    endif

  end function inter4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_interpolate
