! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_delta_eddington

   use tuvx_radiative_transfer_solver, only : radiative_transfer_solver_t,    &
                                              radiation_field_t
   use musica_constants,               only : dk => musica_dk
   use tuvx_constants,                 only : pi

   implicit none

   private
   public :: delta_eddington_t

   type, extends(radiative_transfer_solver_t) :: delta_eddington_t
     ! Radiative flux calculator that applies the delta-Eddington Approximation
     ! (Joseph and Wiscombe, J. Atmos. Sci., 33, 2453-2459, 1976)
   contains
     procedure :: update_radiation_field
   end type delta_eddington_t

   real(dk), parameter :: rZERO = 0.0_dk
   real(dk), parameter :: rONE  = 1.0_dk
   real(dk), parameter :: rTWO  = 2.0_dk
   real(dk), parameter :: d2r   = pi/180._dk

   contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function update_radiation_field( this, solar_zenith_angle, n_streams,       &
      n_layers, spherical_geometry, grid_warehouse, profile_warehouse,        &
      radiator_warehouse ) result( radiation_field )

    use musica_string,                 only : string_t
    use tuvx_diagnostic_util,          only : diagout
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_linear_algebra_linpack,   only : linalgebra_t
    use tuvx_profile,                  only : profile_t
    use tuvx_profile_warehouse,        only : Profile_warehouse_t
    use tuvx_radiator,                 only : radiator_t, radiator_state_t
    use tuvx_radiator_warehouse,       only : radiator_warehouse_t
    use tuvx_radiator_warehouse,       only : warehouse_iterator_t
    use tuvx_spherical_geometry,       only : spherical_geometry_t

    class(delta_eddington_t), intent(inout) :: this ! Delta-Eddington solver

    integer,                    intent(in)    :: n_streams ! not used in delta eddington
    integer,                    intent(in)    :: n_layers  ! number of vertical layers
    real(dk),                   intent(in)    :: solar_zenith_angle ! solar zenith angle [degrees]
    type(grid_warehouse_t),     intent(inout) :: grid_warehouse
    type(Profile_warehouse_t),  intent(inout) :: profile_warehouse
    type(radiator_warehouse_t), intent(inout) :: radiator_warehouse
    type(spherical_geometry_t), intent(inout) :: spherical_geometry

    class(radiation_field_t),   pointer       :: radiation_field

    ! Local variables
    character(len=*), parameter :: Iam = 'Update radiation field: '
    real(dk) :: sum
    real(dk) :: mu
    real(dk) :: tausla( 0 : n_layers ), tauc( 0 : n_layers )
    real(dk) :: mu2( 0 : n_layers )

    ! internal coefficients and matrix
    integer     :: row
    real(dk)    :: lam( n_layers ), taun( n_layers ), bgam( n_layers )
    real(dk)    :: e1( n_layers ), e2( n_layers )
    real(dk)    :: e3( n_layers ), e4( n_layers )
    real(dk)    :: cup( n_layers ), cdn( n_layers )
    real(dk)    :: cuptn( n_layers ), cdntn( n_layers )
    real(dk)    :: mu1( n_layers )
    real(dk)    :: a( 2 * n_layers ), b( 2 * n_layers ), d( 2 * n_layers )
    real(dk)    :: e( 2 * n_layers ), y( 2 * n_layers )

    real(dk) :: pifs, fdn0, surfem, tempg
    real(dk) :: f, g, om
    real(dk) :: gam1, gam2, gam3, gam4
    real(dk) :: gi(n_layers), omi(n_layers)

    integer     :: mrows, lev
    integer     :: i, j
    real(dk) :: expon, expon0, expon1, divisr, temp, up, dn
    real(dk) :: ssfc

    ! Linear algebra package, radiation field type
    type(linalgebra_t) :: linpack

    ! Local variables
    real(dk), parameter                  :: largest = 1.e36_dk
    real(dk), parameter                  :: kfloor = rONE/largest
    real(dk), parameter                  :: precis = 1.e-7_dk
    real(dk), parameter                  :: eps    = 1.e-3_dk

    integer                              :: nlambda, lambdaNdx
    real(dk), allocatable                :: dscat(:,:)
    real(dk), allocatable                :: dscat_accum(:,:)
    real(dk), allocatable                :: dabs_accum(:,:)
    real(dk), allocatable                :: asym_accum(:,:)
    type(warehouse_iterator_t), pointer  :: iter => null( )
    class(radiator_t),          pointer  :: aRadiator
    type(radiator_state_t)               :: atmRadiatorState
    class(grid_t),    pointer            :: zGrid => null( )
    class(grid_t),    pointer            :: lambdaGrid => null( )
    class(profile_t), pointer            :: surfaceAlbedo => null( )

    allocate( radiation_field )

    zGrid => grid_warehouse%get_grid( "height", "km" )
    lambdaGrid => grid_warehouse%get_grid( "wavelength", "nm" )
    surfaceAlbedo => profile_warehouse%get_profile( "surface albedo", "none" )

    nlambda = lambdaGrid%ncells_
    allocate( dscat(n_layers,nlambda) )
    allocate( dscat_accum, mold=dscat )
    allocate( dabs_accum,  mold=dscat )
    allocate( asym_accum,  mold=dscat )
    allocate( radiation_field%edr_(n_layers+1,nlambda) )
    allocate( radiation_field%edn_, mold=radiation_field%edr_ )
    allocate( radiation_field%eup_, mold=radiation_field%edr_ )
    allocate( radiation_field%fdr_, mold=radiation_field%edr_ )
    allocate( radiation_field%fdn_, mold=radiation_field%edr_ )
    allocate( radiation_field%fup_, mold=radiation_field%edr_ )
    ! initialize the accumulators
    dscat_accum = rZERO
    dabs_accum  = rZERO
    asym_accum  = rZERO

    ! iterate over radiators accumulating radiative properties
    iter => radiator_warehouse%get_iterator()
    do while( iter%next() )
      aRadiator => radiator_warehouse%get_radiator( iter )
      associate( OD  => aRadiator%state_%layer_OD_,                           &
                 SSA => aRadiator%state_%layer_SSA_,                          &
                 G   => aRadiator%state_%layer_G_ )
        dscat       = OD * SSA
        dscat_accum = dscat_accum + dscat
        dabs_accum  = dabs_accum + OD*(rONE - SSA)
        asym_accum  = asym_accum + G*dscat
      end associate
      nullify( aRadiator )
    enddo
    deallocate( iter )

    ! set atmosphere radiative properties
    dscat_accum = max( dscat_accum, kfloor )
    dabs_accum  = max( dabs_accum, kfloor )

    atmRadiatorState%layer_OD_ = dscat_accum + dabs_accum
    allocate( atmRadiatorState%layer_SSA_,                                    &
              mold = atmRadiatorState%layer_OD_ )
    where( dscat_accum == kfloor )
      atmRadiatorState%layer_SSA_ = kfloor
    elsewhere
      atmRadiatorState%layer_SSA_ = dscat_accum/atmRadiatorState%layer_OD_
    endwhere

    atmRadiatorState%layer_G_ = asym_accum/dscat_accum

    ! MU = cosine of solar zenith angle
    ! RSFC = surface albedo
    ! TAUU =  unscaled optical depth of each layer
    ! OMU  =  unscaled single scattering albedo
    ! GU   =  unscaled asymmetry factor
    ! N_LAYERS = number of layers in the atmosphere
    ! N_LEVELS = nlayer + 1 = number of levels

    mu = cos( solar_zenith_angle*d2r )
    associate( nid  => spherical_geometry%nid_,                               &
               dsdh => spherical_geometry%dsdh_ )

    wavelength_loop: do lambdaNdx = 1, nlambda
      associate( rsfc => surfaceAlbedo%mid_val_( lambdaNdx ),                 &
             tauu => atmRadiatorState%layer_OD_( n_layers:1:-1, lambdaNdx ),  &
             omu  => atmRadiatorState%layer_SSA_( n_layers:1:-1, lambdaNdx ), &
             gu   => atmRadiatorState%layer_G_( n_layers:1:-1, lambdaNdx ) )

      ! initial conditions:  pi*solar flux = 1;  diffuse incidence = 0
      pifs = rONE
      fdn0 = rZERO
      ! emission at surface (for night light pollution, set pifs = 0, surfem = 1.)
      surfem = rZERO
      !************* compute coefficients for each layer:
      ! GAM1 - GAM4 = 2-stream coefficients, different for different approximations
      ! expON0 = calculation of e when TAU is zero
      ! expON1 = calculation of e when TAU is TAUN
      ! CUP and CDN = calculation when TAU is zero
      ! CUPTN and CDNTN = calc. when TAU is TAUN
      ! DIVISR = prevents division by zero
      tauc   = rZERO
      tausla = rZERO
      mu2    = rONE / sqrt( largest )
      ! delta-scaling. Has to be done for delta-Eddington approximation,
      ! delta discrete ordinate, Practical Improved Flux Method, delta function,
      ! and Hybrid modified Eddington-delta function methods approximations

      do i = 1, n_layers
        f         = gu( i ) * gu( i )
        gi( i )   = ( gu( i ) - f ) / ( rONE - f )
        omi( i )  = ( rONE - f ) * omu( i ) / ( rONE - omu( i ) * f )
        taun( i ) = ( rONE - omu( i ) * f ) * tauu( i )
      end do

      if( lambdaNdx == 1 ) then
        call diagout( 'tauu.new', tauu )
        call diagout( 'gu.new', gu )
        call diagout( 'omu.new', omu )
        call diagout( 'taun.new', taun )
      endif

      ! calculate slant optical depth at the top of the atmosphere when zen>90.
      ! in this case, higher altitude of the top layer is recommended which can
      ! be easily changed in gridz.f.

      if( mu < rZERO ) then
        if( nid(0) < 0 ) then
          tausla(0) = largest
        else
          sum = rZERO
          do j = 1, nid(0)
           sum = sum + rTWO * taun( j ) * dsdh( 0, j )
          end do
          tausla(0) = sum
        end if
      end if

      layer_loop: do i = 1, n_layers

        g  = gi( i )
        om = omi( i )
        tauc( i ) = tauc( i - 1 ) + taun( i )

        ! stay away from 1 by precision.  For g, also stay away from -1

        tempg = min( abs( g ), rONE - precis )
        g = sign( tempg, g )
        om = min( om, rONE - precis )

        ! calculate slant optical depth

        if( nid( i ) < 0 ) then
          tausla( i ) = largest
        else
          sum = rZERO
          do j = 1, min( nid( i ), i )
            sum = sum + taun( j ) * dsdh( i, j )
          enddo
          do j = min( nid( i ), i ) + 1, nid( i )
            sum = sum + rTWO * taun( j ) * dsdh( i, j )
          enddo
          tausla( i ) = sum
          if( tausla( i ) == tausla( i - 1 ) ) then
            mu2( i ) = sqrt( largest )
          else
            mu2( i ) = ( tauc( i ) - tauc( i - 1 ) )                          &
                       / ( tausla( i ) - tausla( i - 1 ) )
            mu2( i ) = sign( max( abs( mu2( i ) ), rONE / sqrt( largest ) ),  &
                             mu2( i ) )
          end if
        end if

        !** the following gamma equations are from pg 16,289, Table 1
        !** save mu1 for each approx. for use in converting irradiance to actinic flux
        ! Eddington approximation(Joseph et al., 1976, JAS, 33, 2452):

        gam1 =   ( 7._dk - om * ( 4._dk + 3._dk * g ) ) / 4._dk
        gam2 = - ( rONE - om * ( 4._dk - 3._dk * g ) ) / 4._dk
        gam3 =   ( rTWO - 3._dk * g * mu ) / 4._dk
        gam4 =   rONE - gam3
        mu1( i ) = 0.5_dk

        lam( i ) = sqrt( gam1 * gam1 - gam2 * gam2 )

        if( gam2 /= rZERO) then
          bgam( i ) = ( gam1 - lam( i ) ) / gam2
        else
          bgam( i ) = rZERO
        endif

        expon = exp( - lam( i ) * taun( i ) )

        ! e1 - e4 = pg 16,292 equation 44

        e1( i ) = rONE + bgam( i ) * expon
        e2( i ) = rONE - bgam( i ) * expon
        e3( i ) = bgam( i ) + expon
        e4( i ) = bgam( i ) - expon

        ! the following sets up for the C equations 23, and 24
        ! found on page 16,290
        ! prevent division by zero (if LAMBDA=1/MU, shift 1/MU^2 by EPS = 1.E-3
        ! which is approx equiv to shifting MU by 0.5*EPS* (MU)**3

        expon0 = exp( -tausla( i - 1 ) )
        expon1 = exp( -tausla( i ) )

        divisr = lam( i ) * lam( i ) - rONE / ( mu2( i ) * mu2( i ) )
        temp   = max( eps, abs( divisr ) )
        divisr = sign( temp, divisr )

        up = om * pifs * ( ( gam1 - rONE / mu2( i ) ) * gam3 + gam4 * gam2 ) &
             / divisr
        dn = om * pifs * ( ( gam1 + rONE / mu2( i ) ) * gam4 + gam2 * gam3 ) &
             / divisr

        ! cup and cdn are when tau is equal to zero
        ! cuptn and cdntn are when tau is equal to taun

        cup(i) = up*expon0
        cdn(i) = dn*expon0
        cuptn(i) = up*expon1
        cdntn(i) = dn*expon1

      enddo layer_loop

      if( lambdaNdx == 1 ) then
        call diagout( 'e1.new',e1 )
        call diagout( 'e2.new',e2 )
        call diagout( 'e3.new',e3 )
        call diagout( 'e4.new',e4 )
        call diagout( 'cup.new',cup )
        call diagout( 'cdn.new',cdn )
        call diagout( 'cuptn.new',cuptn )
        call diagout( 'cdntn.new',cdntn )
        call diagout( 'lam.new',lam )
        call diagout( 'mu2.new',mu2 )
      endif

      !**************** set up matrix ******
      ! ssfc = pg 16,292 equation 37  where pi Fs is one (unity).

      ssfc = rsfc * mu * exp( -tausla( n_layers ) ) * pifs + surfem

      ! MROWS = the number of rows in the matrix

      mrows = 2 * n_layers

      ! the following are from pg 16,292  equations 39 - 43.
      ! set up first row of matrix:

      a(1) = rZERO
      b(1) = e1(1)
      d(1) = -e2(1)
      e(1) = fdn0 - cdn(1)

      ! set up odd rows 3 thru (MROWS - 1):

      i = 0
      do row = 3, mrows - 1, 2
         i = i + 1
         a( row ) = e2( i ) * e3( i ) - e4( i ) * e1( i )
         b( row ) = e1( i ) * e1( i + 1 ) - e3( i ) * e3( i + 1 )
         d( row ) = e3( i ) * e4( i + 1 ) - e1( i ) * e2( i + 1 )
         e( row ) = e3( i ) * ( cup( i + 1 ) - cuptn( i ) )                   &
                    + e1( i ) * ( cdntn( i ) - cdn( i + 1 ) )
      enddo

      ! set up even rows 2 thru (MROWS - 2):

      i = 0
      do row = 2, mrows - 2, 2
         i = i + 1
         a( row ) = e2( i + 1 ) * e1( i ) - e3( i ) * e4( i + 1 )
         b( row ) = e2( i ) * e2( i + 1 ) - e4( i ) * e4( i + 1 )
         d( row ) = e1( i + 1 ) * e4( i + 1 ) - e2( i + 1) * e3( i + 1 )
         e( row ) = ( cup( i + 1 ) - cuptn( i ) ) * e2( i + 1 )               &
                    - ( cdn( i + 1 ) - cdntn( i ) ) * e4( i + 1 )
      enddo

      ! set up last row of matrix at MROWS:

      a( mrows ) = e1( n_layers ) - rsfc * e3( n_layers )
      b( mrows ) = e2( n_layers ) - rsfc * e4( n_layers )
      d( mrows ) = rZERO
      e( mrows ) = ssfc - cuptn( n_layers ) + rsfc * cdntn( n_layers )

      if( lambdaNdx == 1 ) then
        call diagout( 'a.new', a )
        call diagout( 'b.new', b )
        call diagout( 'd.new', d )
        call diagout( 'e.new', e )
        call diagout( 'tausla.new', tausla )
      endif

      ! solve tri-diagonal system:

      y = linpack%tridiag( a, b, d, e )

      !*** unfold solution of matrix, compute output fluxes:
      ! the following equations are from pg 16,291  equations 31 & 32

      associate( edr => radiation_field%edr_( :, lambdaNdx ),                 &
                 eup => radiation_field%eup_( :, lambdaNdx ),                 &
                 edn => radiation_field%edn_( :, lambdaNdx ),                 &
                 fdr => radiation_field%fdr_( :, lambdaNdx ),                 &
                 fup => radiation_field%fup_( :, lambdaNdx ),                 &
                 fdn => radiation_field%fdn_( :, lambdaNdx ) )
      fdr(1) = pifs * exp( -tausla(0) )
      edr(1) = mu * fdr(1)
      edn(1) = fdn0
      eup(1) =  y(1) * e3(1) - y(2) * e4(1) + cup(1)
      fdn(1) = edn(1) / mu1(1)
      fup(1) = eup(1) / mu1(1)

      j   = 1
      row = 1
      do lev = 2, n_layers + 1
         fdr( lev ) = pifs * exp( -tausla( lev - 1 ) )
         edr( lev ) = mu * fdr( lev )
         edn( lev ) = y( row ) * e3( j ) + y( row + 1 ) * e4( j ) + cdntn( j )
         eup( lev ) = y( row ) * e1( j ) + y( row + 1 ) * e2( j ) + cuptn( j )
         fdn( lev ) = edn( lev ) / mu1( j )
         fup( lev ) = eup( lev ) / mu1( j )

         row = row + 2
         j = j + 1
      enddo
      ! transform from top-down to buttom-up
      fdr = fdr( n_layers + 1 : 1 : -1 )
      fup = fup( n_layers + 1 : 1 : -1 )
      fdn = fdn( n_layers + 1 : 1 : -1 )
      edr = edr( n_layers + 1 : 1 : -1 )
      eup = eup( n_layers + 1 : 1 : -1 )
      edn = edn( n_layers + 1 : 1 : -1 )

      end associate

    end associate

    enddo wavelength_loop

    end associate

    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( surfaceAlbedo )

  end function update_radiation_field

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_delta_eddington
