      module DISORD_SUBS  
  
      use musica_constants, only : ik => musica_ik, dk => musica_dk, lk => musica_lk
      use tuv_params, only : PI  
  
      IMPLICIT NONE  
  
      private  
      public :: PSNDO  
  
      integer(ik), PARAMETER :: MXCLY = 151_ik  
      integer(ik), PARAMETER :: MXULV = 151_ik  
      integer(ik), PARAMETER :: MXCMU = 32_ik
      integer(ik), PARAMETER :: MXUMU = 32_ik
      integer(ik), PARAMETER :: MXPHI = 3_ik  
      real(dk), PARAMETER    :: rZERO = 0.0_dk
      real(dk), PARAMETER    :: ONEHALF  = 0.5_dk
      real(dk), PARAMETER    :: rONE  = 1.0_dk
      real(dk), PARAMETER    :: rTWO  = 2.0_dk
      real(dk), PARAMETER    :: TWOPI  = rTWO*PI
      real(dk), PARAMETER    :: FOURPI  = 4._dk*PI
! Discrete ordinate constants:  
! For pseudo-spherical DISORT, PLANK, USRTAU and USRANG must be .FALSE.;  
! ONLYFL must be .TRUE.; FBEAM = 1.; FISOT = 0.; IBCND = 0  
      integer(ik), parameter :: iZERO = 0_ik  
      integer(ik), parameter :: iONE  = 1_ik  
      integer(ik), parameter :: iTWO  = 2_ik
      integer(ik), parameter :: iTHREE  = 3_ik
      integer(ik), parameter :: NPHI  = iZERO
      integer(ik), parameter :: IBCND = iZERO
      real(dk), parameter    :: ACCUR = 0.0001_dk
      real(dk), parameter    :: FBEAM = rONE
      real(dk), parameter    :: FISOT = rZERO  
      real(dk), parameter    :: PHI0  = rZERO  
      logical(lk), parameter :: DELTAM = .true._lk  
      logical(lk), parameter :: LAMBER = .true._lk  
      logical(lk), parameter :: PLANK = .FALSE._lk  
      logical(lk), parameter :: USRANG = .FALSE._lk  
      logical(lk), parameter :: USRTAU = .FALSE._lk  
      logical(lk), parameter :: ONLYFL = .TRUE._lk  
      logical(lk), parameter :: PRNT(7) = .FALSE._lk  
  
      real(dk) :: DITHER  
  
      contains  
  
      SUBROUTINE PSNDO( dsdh, nid, &
                        NLYR, DTAUC, SSALB, PMOM, &
                        ALBEDO, NSTR, UMU0, &
                        RFLDIR, RFLDN, FLUP, &
                        uavgso, uavgup, uavgdn )  
  
      integer(ik), intent(in) :: NLYR  
      integer(ik), intent(in) :: NSTR  
      integer(ik), intent(in) :: nid(0:)  
      real(dk), intent(in)    :: ALBEDO  
      real(dk), intent(in)    :: UMU0  
      real(dk), intent(in)    :: dsdh(0:,:)  
      real(dk), intent(in)    :: DTAUC(:)  
      real(dk), intent(inout) :: SSALB(:)  
      real(dk), intent(inout) :: PMOM(0:,:)  
      real(dk), intent(out)   :: RFLDIR(:), RFLDN(:), FLUP(:)  
      real(dk), intent(out)   :: uavgso(:)  
      real(dk), intent(out)   :: uavgup(:)  
      real(dk), intent(out)   :: uavgdn(:)  
  
      real(dk), PARAMETER     :: RPD  = PI / 180.0_dk
!     spherical geometry  
      real(dk) tausla(0:NLYR), tauslau(0:NLYR), mu2(0:NLYR)  
!     ..  
  
!     local variables  
      real(dk) :: &
          DFDT(NLYR+iONE), HL(iZERO:NSTR), PHI(MXPHI), &
          TRNMED(NSTR), U0U(nstr,nlyr+iONE), UAVG(NLYR+iONE), &
          UMU(NSTR), CWT(NSTR), UTAU(NLYR+iONE), UU(NSTR,NLYR,MXPHI)  
!     ..  
!     .. Local Scalars ..  
  
      logical(lk) :: COMPAR, LYRCUT, PASS1  
      integer(ik) :: IQ, IU, J, KCONV, L, LC, LEV, LU, MAZIM, NAZ, NCOL, &
                     NCOS, NCUT, NN  
      real(dk)    :: AZERR, AZTERM, BPLANK, COSPHI, DELM0, DUM, SGN, TPLANK  
!     ..  
!     .. Local Arrays ..  
  
      integer(ik) :: NTAU, NUMU  
      integer(ik) :: IPVT(NSTR*NLYR), LAYRU(NLYR+iONE)  
  
      real(dk) ::   ANGCOS(iONE)  
      real(dk) ::   AMB(NSTR/iTWO,NSTR/iTWO), APB(NSTR/iTWO,NSTR/iTWO), ARRAY(NSTR,NSTR), &
                B(NSTR*NLYR), BDR(NSTR/iTWO,0:NSTR/iTWO), BEM(NSTR/iTWO), &
                CBAND(9_ik*(NSTR/iTWO)-iTWO,NSTR*NLYR), CC(NSTR,NSTR), &
                CMU(NSTR), DTAUCP(NLYR), &
                EMU(NSTR), EVAL(NSTR/iTWO), EVECC(NSTR, NSTR), &
                EXPBEA(iZERO:NLYR), FLDIR(NLYR+iONE), FLDN(NLYR+iONE), &
                FLYR(NLYR), GC(NSTR,NSTR,NLYR), &
                GL(iZERO:NSTR,NLYR), GU(NSTR,NSTR,NLYR), &
                HLPR(iZERO:NSTR), KK(NSTR,NLYR), LL(NSTR,NLYR), &
                OPRIM(NLYR), PHIRAD(MXPHI), PKAG(iZERO:NLYR), &
                PSI(NSTR), RMU(NSTR,iZERO:NSTR/iTWO), TAUC(iZERO:NLYR), &
                TAUCPR(iZERO:NLYR), U0C(NSTR,NLYR+iONE), UTAUPR(NLYR+iONE), &
                UUM(NSTR,NLYR), WK(NSTR), XR0(NLYR), &
                XR1(NLYR), YLM0(iZERO:NSTR,iONE), YLMC(iZERO:NSTR,NSTR), &
                YLMU(iZERO:NSTR,NSTR), Z(NSTR*NLYR), Z0(NSTR), &
                Z0U(NSTR,NLYR), Z1(NSTR), Z1U(NSTR,NLYR), &
                ZBEAM(NSTR,NLYR), ZJ(NSTR), &
                ZPLK0(NSTR,NLYR), ZPLK1(NSTR,NLYR), ZZ(NSTR,NLYR)  
  
      real(dk) :: sindir(nlyr+iONE), sinup(nlyr+iONE), sindn(nlyr+iONE)  
  
!gy added glsave and dgl to allow adjustable dimensioning in SOLVEC  
      real(dk) GLSAVE(iZERO:nstr), DGL(iZERO:nstr)  
  
      real(dk) :: AAD(NSTR/iTWO,NSTR/iTWO), EVECCD(NSTR/iTWO,NSTR/iTWO), &
                  EVALD(NSTR/iTWO), WKD(NSTR)  
  
      real(dk)    :: PLKAVG  
  
      SAVE      PASS1  
      DATA      PASS1 / .TRUE._lk /  
  
      IF( PASS1 ) THEN  
!        DITHER = 10._dk*R1MACH( 4 )  
         DITHER = 10._dk*EPSILON( rONE )
!                            ** Must dither more on Cray (14-digit prec)  
         IF( DITHER < 1.E-10_dk ) DITHER = 10._dk*DITHER  
         PASS1 = .FALSE._ik
      END IF  
   
!     ** Calculate cumulative optical depth  
!     and dither single-scatter albedo  
!     to improve numerical behavior of  
!     eigen{value/vector} computation  
  
      TAUC = rZERO  
  
      DO LC = iONE, NLYR  
        IF( SSALB(LC) == rONE ) THEN  
           SSALB(LC) = rONE - DITHER  
        ENDIF  
        TAUC(LC) = TAUC(LC - iONE) + DTAUC(LC)  
      ENDDO  

!                                ** Check input dimensions and variables  
      CALL CHEKIN( NLYR, DTAUC, SSALB, PMOM, &
                   NTAU, UTAU, NSTR, NUMU, UMU, NPHI, &
                   PHI, UMU0, FISOT, ALBEDO, &
                   HL, TAUC )  
  
!                                 ** Zero internal and output arrays  
      CALL  ZEROAL( EXPBEA(iONE:), FLYR, OPRIM, TAUCPR(iONE:), XR0, XR1, &
                    CMU, CWT, PSI, WK, Z0, Z1, ZJ, &
                    HLPR, YLM0(:,1), ARRAY, CC, EVECC, &
                    GL, YLMC, YLMU, &
                    KK, LL, ZZ, ZPLK0, ZPLK1, &
                    GC, LAYRU, UTAUPR, &
                    GU, Z0U, Z1U, ZBEAM, &
                    EVAL, AMB, APB, IPVT, Z, &
                    RFLDIR, RFLDN, FLUP, UAVG, DFDT, &
                    TRNMED, U0U, UU )  
  
!                                 ** Perform various setup operations  
      CALL SETDIS( &
          dsdh, nid, tausla, tauslau, mu2, &
          CMU, CWT, DELTAM, DTAUC, DTAUCP, EXPBEA, FBEAM, FLYR, &
          GL, HL, HLPR, IBCND, LAMBER, LAYRU, LYRCUT, &
          NCUT, NLYR, NTAU, NN, NSTR, PLANK, &
          NUMU, ONLYFL, OPRIM, PMOM, SSALB, TAUC, TAUCPR, UTAU, &
          UTAUPR, UMU, UMU0, USRTAU, USRANG )  
  
!                                 ** Print input information  
      IF ( PRNT(iONE) ) THEN  
        CALL PRTINP( NLYR, DTAUC, DTAUCP, SSALB, PMOM, &
                     NTAU, UTAU, NSTR, NUMU, UMU, &
                     PHI, UMU0, FISOT, ALBEDO, HL,& 
                     FLYR, LYRCUT, &
                     OPRIM, TAUC, TAUCPR, PRNT(7_ik) )  
      ENDIF  
  
!                              ** Handle special case for getting albedo  
!                                 and transmissivity of medium for many  
!                                 beam angles at once  
!                                   ** Calculate Planck functions  
  
         BPLANK = rZERO  
         TPLANK = rZERO  
         PKAG   = rZERO  
  
! ========  BEGIN LOOP TO SUM AZIMUTHAL COMPONENTS OF INTENSITY  =======  
!           (EQ STWJ 5)  
  
      KCONV = iZERO  
!                                    ** Azimuth-independent case  
      IF( FBEAM == rZERO .OR. (rONE - UMU0) < 1.E-5_dk .OR. ONLYFL &
                         .OR. (NUMU == iONE .AND. (rONE - UMU(iONE)) < 1.E-5_dk) ) THEN  
        NAZ = iZERO
      ELSE  
        NAZ  = NSTR - iONE
      ENDIF  
  
      AZIMUTH_LOOP: DO MAZIM = iZERO, NAZ  
         IF( MAZIM == iZERO ) THEN  
           DELM0  = rONE  
         ELSE  
           DELM0  = rZERO  
         ENDIF  
!                             ** Get normalized associated Legendre  
!                                polynomials for  
!                                (a) incident beam angle cosine  
!                                (b) computational and user polar angle  
!                                    cosines  
         IF( FBEAM > rZERO ) THEN  
            NCOS   = iONE
            ANGCOS = -UMU0  
            CALL LEPOLY( NCOS, MAZIM, NSTR - iONE, ANGCOS, YLM0 )  
         END IF  
  
  
         IF( .NOT. ONLYFL .AND. USRANG ) THEN  
            CALL LEPOLY( NUMU, MAZIM, NSTR-iONE, UMU, YLMU )  
         ENDIF  
  
         CALL LEPOLY( NN, MAZIM, NSTR-iONE, CMU, YLMC )  
  
!                       ** Get normalized associated Legendre polys.  
!                          with negative arguments from those with  
!                          positive arguments; Dave/Armstrong Eq. (15)  
         SGN  = -rONE  
  
         DO L = MAZIM, NSTR - iONE  
            SGN  = - SGN  
            DO IQ = NN + iONE, NSTR  
               YLMC(L,IQ) = SGN*YLMC(L,IQ - NN)  
            ENDDO  
         ENDDO  
!     ** Specify users bottom reflectivity  
!        and emissivity properties  
         IF ( .NOT. LYRCUT ) THEN  
           CALL  SURFAC( &
               ALBEDO, DELM0, FBEAM, HLPR, LAMBER, &
               MAZIM, NN, NUMU, NSTR, ONLYFL, &
               UMU, USRANG, YLM0(:,1), YLMC, YLMU, BDR, EMU, BEM, RMU )  
         ENDIF  
  
  
! ===================  BEGIN LOOP ON COMPUTATIONAL LAYERS  =============  
  
         DO LC = iONE, NCUT  
            CALL SOLVEC( AMB, APB, ARRAY, CMU, CWT, GL(iZERO:,LC), &
                 MAZIM, NN, NSTR, YLM0(:,iONE), YLMC, CC, &
                 EVECC, EVAL, KK(:,LC ), GC(:,:,LC), AAD, EVECCD, &
                 EVALD, WK, WKD, DELM0, FBEAM, IPVT, PI, UMU0, &
                 ZJ, ZZ(:,LC), OPRIM(LC), LC, mu2(lc), glsave, dgl)  
         ENDDO  
  
  
! ===================  END LOOP ON COMPUTATIONAL LAYERS  ===============  
  
!                      ** Set coefficient matrix of equations combining  
!                         boundary and layer interface conditions  
  
         CALL SETMTX( BDR, CBAND, CMU, CWT, DELM0, DTAUCP, GC, KK, &
                      LAMBER, LYRCUT, NCOL, NCUT, &
                      NN, NSTR, TAUCPR, WK )  
  
!                      ** Solve for constants of integration in homo-  
!                         geneous solution (general boundary conditions)  
         CALL SOLVE0( &
             B, BDR, BEM, BPLANK, CBAND, CMU, CWT, EXPBEA, &
             FBEAM, FISOT, IPVT, LAMBER, LL, LYRCUT, MAZIM, &
             NCOL, NCUT, NN, NSTR, PI, &
             TPLANK, TAUCPR, UMU0, Z, ZZ, ZPLK0, ZPLK1 )  
  
!                                  ** Compute upward and downward fluxes  
         IF ( MAZIM == iZERO ) THEN  
            CALL FLUXES( tausla, tauslau, &
                         CMU, CWT, FBEAM, GC, KK, LAYRU, LL, LYRCUT, &
                         NCUT, NN, NSTR, NTAU, &
                         PI, PRNT, SSALB, TAUCPR, UMU0, UTAU, UTAUPR, &
                         XR0, XR1, ZZ, ZPLK0, ZPLK1, DFDT, FLUP, &
                         FLDN, FLDIR, RFLDIR, RFLDN, UAVG, U0C, &
                         uavgso, uavgup, uavgdn, sindir, sinup, sindn )
         ENDIF  
  
         IF( ONLYFL ) THEN  
            IF( MXUMU >= NSTR ) THEN  
!     ** Save azimuthal-avg intensities  
!        at quadrature angles  
               DO LU = iONE, NTAU  
                  DO IQ = iONE, NSTR  
                     U0U(IQ,LU) = U0C(IQ,LU)  
                  ENDDO  
               ENDDO  
            ENDIF  
            EXIT AZIMUTH_LOOP  
         ENDIF  
  
         UUM = rZERO  
  
         IF( MAZIM == iZERO ) THEN  
!     ** Save azimuthally averaged intensities  
            DO LU = iONE, NTAU  
               DO IU = iONE, NUMU  
                  U0U(IU,LU) = UUM(IU,LU)  
                  DO J = iONE, NPHI  
                     UU(IU,LU,J) = UUM(IU,LU)  
                  ENDDO  
               ENDDO  
            ENDDO  
!                              ** Print azimuthally averaged intensities  
!                                 at user angles  
  
            IF( PRNT( 4 ) ) THEN  
               CALL PRAVIN( UMU, NUMU, UTAU, NTAU, U0U )  
            ENDIF  
            IF( NAZ > iZERO ) THEN  
               PHIRAD = rZERO  
               DO J = iONE, NPHI  
                  PHIRAD(J) = RPD*(PHI(J) - PHI0)  
               ENDDO  
            END IF  
         ELSE  
!                                ** Increment intensity by current  
!                                   azimuthal component (Fourier  
!                                   cosine series);  Eq SD(2)  
            AZERR  = rZERO  
            DO J = iONE, NPHI  
               COSPHI = COS( MAZIM*PHIRAD(J) )  
               DO LU = iONE, NTAU  
                  DO IU = iONE, NUMU  
                     AZTERM = UUM(IU,LU)*COSPHI  
                     UU(IU,LU,J) = UU(IU,LU,J) + AZTERM  
                     AZERR = MAX( AZERR, RATIO( ABS(AZTERM), ABS(UU(IU,LU,J))) )  
                  ENDDO  
               ENDDO  
            ENDDO  
  
            IF( AZERR <= ACCUR ) KCONV  = KCONV + iONE
  
            IF( KCONV >= iTWO ) THEN  
               EXIT AZIMUTH_LOOP  
            ENDIF  
         ENDIF  
  
      ENDDO AZIMUTH_LOOP  
  
! ===================  END LOOP ON AZIMUTHAL COMPONENTS  ===============  
  
!                                          ** Print intensities  
      IF( PRNT( 5 ) .AND. .NOT. ONLYFL ) THEN  
        CALL PRTINT( UU, UTAU, NTAU, UMU, NUMU, PHI )  
      ENDIF  
  
      END SUBROUTINE PSNDO  
  
      SUBROUTINE ASYMTX( AA, EVEC, EVAL, IER, WKD, AAD, &
                         EVECD, EVALD )  
  
!    =======  D O U B L E    P R E C I S I O N    V E R S I O N  ======  
  
!       Solves eigenfunction problem for real asymmetric matrix  
!       for which it is known a priori that the eigenvalues are real.  
  
!       This is an adaptation of a subroutine EIGRF in the IMSL  
!       library to use real instead of complex arithmetic, accounting  
!       for the known fact that the eigenvalues and eigenvectors in  
!       the discrete ordinate solution are real.  Other changes include  
!       putting all the called subroutines in-line, deleting the  
!       performance index calculation, updating many DO-loops  
!       to Fortran77, and in calculating the machine precision  
!       TOL instead of specifying it in a data statement.  
  
!       EIGRF is based primarily on EISPACK routines.  The matrix is  
!       first balanced using the Parlett-Reinsch algorithm.  Then  
!       the Martin-Wilkinson algorithm is applied.  
  
!       References:  
!          Dongarra, J. and C. Moler, EISPACK -- A Package for Solving  
!             Matrix Eigenvalue Problems, in Cowell, ed., 1984:  
!             Sources and Development of Mathematical Software,  
!             Prentice-Hall, Englewood Cliffs, NJ  
!         Parlett and Reinsch, 1969: Balancing a Matrix for Calculation  
!             of Eigenvalues and Eigenvectors, Num. Math. 13, 293-304  
!         Wilkinson, J., 1965: The Algebraic Eigenvalue Problem,  
!             Clarendon Press, Oxford  
  
!   I N P U T    V A R I A B L E S:  
  
!       AA    :  input asymmetric matrix, destroyed after solved  
!        M    :  order of  AA  
!       IA    :  first dimension of  AA  
!    IEVEC    :  first dimension of  EVEC  
  
!   O U T P U T    V A R I A B L E S:  
  
!       EVEC  :  (unnormalized) eigenvectors of  AA  
!                   ( column J corresponds to EVAL(J) )  
  
!       EVAL  :  (unordered) eigenvalues of AA ( dimension at least M )  
  
!       IER   :  if .NE. 0, signals that EVAL(IER) failed to converge;  
!                   in that case eigenvalues IER+1,IER+2,...,M  are  
!                   correct but eigenvalues 1,...,IER are set to zero.  
  
!   S C R A T C H   V A R I A B L E S:  
  
!       WKD   :  work area ( dimension at least 2*M )  
!       AAD   :  double precision stand-in for AA  
!       EVECD :  double precision stand-in for EVEC  
!       EVALD :  double precision stand-in for EVAL  
  
!   Called by- SOLEIG  
!   Calls- D1MACH, ERRMSG  
! +-------------------------------------------------------------------+  
  
!     .. Scalar Arguments ..  
  
      integer(ik), intent(out) :: IER  
!     ..  
!     .. Array Arguments ..  
  
      real(dk), intent(in)  :: AA(:,:)  
      real(dk), intent(out) :: EVAL(:), EVEC(:,:)  
      real(dk), intent(inout) :: WKD(:)  
      real(dk), intent(out)   :: AAD(:,:)  
      real(dk), intent(out)   :: EVALD(:), EVECD(:,:)  
!     ..  
!     .. Local Scalars ..  
  
      integer(ik) :: I, II, IN, J, K, KA, KKK, L, LB, LLL, N, N1, N2  
      integer(ik) :: M, IA, IEVEC  
      logical(lk) :: Converged, NOTLAS, Reset
  
      real(dk) :: COL, DISCRI, F, G, H,  
      real(dk) :: P, Q, R, REPL, RNORM, ROW, S, SCALE, SGN, T
      real(dk) :: TOL, UU, VV, W, X, Y, Z
  
      real(dk), parameter :: C1 = 0.4375_8_dk
      real(dk), parameter :: C2 = 0.5_8_dk
      real(dk), parameter :: C3 = 0.75_8_dk
      real(dk), parameter :: C4 = 0.95_8_dk
      real(dk), parameter :: C5 = 16._8_dk
      real(dk), parameter :: C6 = 256._8_dk
      real(dk), parameter :: ZERO = rZERO
      real(dk), parameter :: ONE  = rONE
  
      IER  = iZERO
!     TOL  = D1MACH( 4 )  
      TOL  = EPSILON( rONE )
  
      M     = size(EVAL)  
      IA    = size(AA,dim=1)  
      IEVEC = size(EVECD,dim=1)  
  
      IF( M < iONE .OR. IA < M .OR. IEVEC < M ) &
        CALL ERRMSG( 'ASYMTX--bad input variable(s)', .TRUE. )  
  
!                           ** Handle 1x1 and 2x2 special cases  
      IF( M == iONE ) THEN  
         EVAL(iONE)      = AA(iONE,iONE)  
         EVEC(iONE,iONE) = rONE  
         RETURN  
      ELSE IF( M == iTWO ) THEN  
         DISCRI = ( AA(iONE,iONE) - AA(iTWO,iTWO) )**2  + 4._dk*AA(iONE,iTWO)*AA(iTWO,iONE)  
         IF( DISCRI < rZERO ) &
             CALL ERRMSG( 'ASYMTX--complex evals in 2x2 case',.TRUE. )  
  
         SGN  = rONE  
         IF( AA(iONE,iONE) < AA(iTWO,iTWO) ) SGN  = -rONE  
  
         EVAL(iONE) = 0.5_dk*(AA(iONE,iONE) + AA(iTWO,iTWO) + SGN*SQRT( DISCRI ) )  
         EVAL(iTWO) = 0.5_dk*(AA(iONE,iONE) + AA(iTWO,iTWO) - SGN*SQRT( DISCRI ) )  
         EVEC(iONE,iONE) = rONE  
         EVEC(iTWO,iTWO) = rONE  
  
         IF( AA(iONE,iONE) == AA(iTWO,iTWO) &
             .AND. (AA(iTWO,iONE) == rZERO .OR. AA(iONE,iTWO) == rZERO) ) THEN  
  
            RNORM  = ABS( AA(iONE,iONE) ) + ABS( AA(iONE,iTWO) ) &
                   + ABS( AA(iTWO,iONE) ) + ABS( AA(iTWO,iTWO) )  
            W  = TOL*RNORM  
            EVEC(iTWO,iONE) =   AA(iTWO,iONE) / W  
            EVEC(iONE,iTWO) = - AA(iONE,iTWO) / W  
         ELSE  
            EVEC(iTWO,iONE) = AA(iTWO,iONE) / (EVAL(iONE) - AA(iTWO,iTWO))  
            EVEC(iONE,iTWO) = AA(iONE,iTWO) / (EVAL(iTWO) - AA(iONE,iONE))  
         END IF  

         RETURN  
      END IF  
!                               ** Put s.p. matrix into d.p. matrix  
      DO K = iONE, M  
        AAD(:M,K) = AA(:M,K)
      ENDDO  
!                                ** Initialize output variables  
      IER  = iZERO
  
      DO I = iONE, M  
        EVALD(I) = ZERO  
        DO J = iONE, M  
          EVECD(I,J) = ZERO  
        ENDDO
        EVECD(I,I) = ONE  
      ENDDO
  
!                  ** Balance the input matrix and reduce its norm by  
!                     diagonal similarity transformation stored in WK;  
!                     then search for rows isolating an eigenvalue  
!                     and push them down  
      RNORM  = ZERO  
      L  = iONE
      K  = M  
  
Balance_loop: &
      DO
        KKK  = K  
        DO J = KKK, iONE, -iONE  
          ROW  = ZERO  
          DO I = iONE, K  
            IF( I /= J ) ROW  = ROW + ABS( AAD(J,I) )  
          ENDDO
  
          IF( ROW == ZERO ) THEN  
            WKD(K) = J  
  
            IF( J /= K ) THEN  
               DO I = iONE, K  
                  REPL     = AAD(I,J)  
                  AAD(I,J) = AAD(I,K)  
                  AAD(I,K) = REPL  
               ENDDO
  
               DO I = L, M  
                  REPL     = AAD(J,I)  
                  AAD(J,I) = AAD(K,I)  
                  AAD(K,I) = REPL  
               ENDDO
            END IF  
  
            K  = K - iONE  
            CYCLE Balance_loop
          ENDIF  
        ENDDO
        EXIT Balance_loop
      ENDDO Balance_loop
   
!                                ** Search for columns isolating an  
!                                   eigenvalue and push them left  
Column_search_loop: &
      DO
        LLL  = L  
        DO J = LLL, K  
          COL  = ZERO  
          DO I = L, K
            IF( I /= J ) COL = COL + ABS( AAD(I,J) )  
          ENDDO
  
          IF( COL == ZERO ) THEN  
            WKD( L ) = J  
            IF( J /= L ) THEN  
               DO I = iONE, K  
                  REPL     = AAD(I,J)  
                  AAD(I,J) = AAD(I,L)  
                  AAD(I,L) = REPL  
               ENDDO
  
               DO I = L, M  
                  REPL     = AAD(J,I)  
                  AAD(J,I) = AAD(L,I)  
                  AAD(L,I) = REPL  
               ENDDO
            END IF  
            L  = L + iONE
            CYCLE Column_search_loop
          ENDIF  
        ENDDO
        EXIT Column_search_loop
      ENDDO Column_search_loop
!                           ** Balance the submatrix in rows L through K  
      WKD(L:K) = rONE  
  
      Converged = .TRUE._lk
Convergence_loop: &
      DO 
        DO I = L, K  
          COL  = ZERO  
          ROW  = ZERO  
  
          DO J = L, K  
            IF( J /= I ) THEN  
              COL  = COL + ABS( AAD(J,I) )  
              ROW  = ROW + ABS( AAD(I,J) )  
            END IF  
          ENDDO
  
          F  = ONE  
          G  = ROW / C5  
          H  = COL + ROW  
          DO WHILE( COL < G )
            F    = F*C5  
            COL  = COL*C6  
          ENDDO
  
          G  = ROW*C5  
          DO WHILE( COL >= G )
            F    = F / C5  
            COL  = COL / C6  
          ENDDO
!                                                ** Now balance  
          IF( (COL + ROW) / F < C4*H ) THEN  
            WKD( I ) = WKD( I )*F  
            Converged = .FALSE._lk
  
            DO J = L, M  
              AAD(I,J) = AAD(I,J) / F  
            ENDDO
  
            DO J = iONE, K  
              AAD(J,I) = AAD(J,I)*F  
            ENDDO
          END IF  
        ENDDO
        IF( Converged ) THEN
          EXIT Convergence_loop
        ENDIF
        Converged = .TRUE._lk
      ENDDO Convergence_loop
  
!                                   ** Is A already in Hessenberg form?  
Is_Hessenberg: &
      IF( K-iONE >= L+iONE ) THEN
!                                   ** Transfer A to a Hessenberg form  
        DO N = L + iONE, K - iONE  
          H  = rZERO  
          WKD(N + M) = rZERO  
          SCALE  = rZERO  
!                                                 ** Scale column  
          DO I = N, K  
            SCALE  = SCALE + ABS( AAD(I,N - iONE) )  
          ENDDO
  
          IF( SCALE /= ZERO ) THEN  
            DO I = K, N, -iONE  
              WKD(I + M) = AAD(I,N - iONE) / SCALE  
              H  = H + WKD(I + M)**2  
            ENDDO
  
            G    = - SIGN( SQRT(H), WKD(N + M) )  
            H    = H - WKD(N + M)*G  
            WKD(N + M) = WKD(N + M) - G  
!                                            ** Form (I-(U*UT)/H)*A  
            DO J = N, M  
               F  = ZERO  
  
               DO I = K, N, -iONE
                 F  = F + WKD(I + M)*AAD(I,J)  
               ENDDO
  
               DO I = N, K  
                 AAD(I,J) = AAD(I,J) - WKD(I + M)*F / H  
               ENDDO
            ENDDO
!                                    ** Form (I-(U*UT)/H)*A*(I-(U*UT)/H)  
            DO I = iONE, K  
               F  = ZERO  
  
               DO J = K, N, -iONE
                  F  = F + WKD(J + M)*AAD(I,J)  
               ENDDO
  
               DO J = N, K  
                  AAD(I,J) = AAD(I,J) - WKD(J + M)*F / H  
               ENDDO
            ENDDO
  
            WKD(N + M) = SCALE*WKD( N + M )  
            AAD(N,N - iONE) = SCALE*G  
          ENDIF  
        ENDDO
  
        DO N = K - iTWO, L, -iONE
          N1   = N + iONE
          N2   = N + iTWO  
          F  = AAD(N + iONE,N)  
  
          IF( F /= ZERO ) THEN  
            F  = F*WKD( N + iONE + M )  
  
            DO I = N + iTWO, K  
               WKD( I + M ) = AAD( I, N )  
            ENDDO
  
            IF( (N + iONE) <= K ) THEN  
               DO J = iONE, M  
                  G  = ZERO  
                  DO I = N + iONE, K  
                     G  = G + WKD(I + M)*EVECD(I,J)  
                  ENDDO
                  G  = G / F  
                  DO I = N + iONE, K  
                     EVECD(I,J) = EVECD(I,J) + G*WKD(I + M)  
                  ENDDO
               ENDDO
            END IF  
          END IF  
        ENDDO
      ENDIF Is_Hessenberg
  
      N  = iONE
  
      DO I = iONE, M
         DO J = N, M  
           RNORM  = RNORM + ABS( AAD(I,J) )  
         ENDDO
  
         N  = I  
         IF( I < L .OR. I > K ) EVALD(I) = AAD(I,I)  
      ENDDO
  
      N  = K  
      T  = ZERO  
  
!                                      ** Search for next eigenvalues  
      Reset = .true._lk
Eigen_loop: &
      DO
        IF( N >= L ) THEN
          IF( Reset ) THEN
            IN  = iZERO
            N1  = N - iONE  
            N2  = N - iTWO  
          ENDIF
!                          ** Look for single small sub-diagonal element  
        DO I = L, N  
          LB  = N + L - I  
          IF( LB == L ) EXIT
          S  = ABS( AAD(LB - iONE,LB - iONE) ) + ABS( AAD(LB,LB) )  
          IF( S == ZERO ) S  = RNORM  
          IF( ABS( AAD(LB,LB-iONE) ) <= TOL*S ) EXIT
        ENDDO
  
        X  = AAD(N,N)  
  
        IF( LB == N ) THEN  
!                                        ** One eigenvalue found  
          AAD(N,N) = X + T  
          EVALD(N) = AAD(N,N)  
          N  = N1  
          Reset = .true._lk
          CYCLE Eigen_loop
        END IF  
  
! next line has been included to avoid run time error caused by xlf  
        IF( N1 <= iZERO .OR. N <= iZERO ) THEN  
          WRITE(0,*) 'Subscript out of bounds in ASYMTX'  
          STOP 9999  
        ENDIF  
  
        Y = AAD(N1,N1)  
        W  = AAD(N,N1)*AAD(N1,N)  
        IF( LB == N1 ) THEN  
!                                        ** Two eigenvalues found  
          P  = ( Y - X )*C2  
          Q  = P**2 + W  
          Z  = SQRT( ABS( Q ) )  
          AAD(N,N) = X + T  
          X  = AAD(N,N)  
          AAD(N1,N1) = Y + T  
!                                        ** Real pair  
          Z  = P + SIGN( Z, P )  
          EVALD(N1) = X + Z  
          EVALD(N)  = EVALD(N1)  
          IF( Z /= ZERO ) EVALD(N) = X - W / Z  
          X  = AAD(N,N1)  
!                                  ** Employ scale factor in case  
!                                     X and Z are very small  
          R  = SQRT( X*X + Z*Z )  
          P  = X / R  
          Q  = Z / R  
!                                             ** Row modification  
          DO J = N1, M  
            Z  = AAD(N1,J)  
            AAD(N1,J) = Q*Z + P*AAD(N,J)  
            AAD(N,J)  = Q*AAD(N,J) - P*Z  
          ENDDO
!                                             ** Column modification  
          DO I = iONE, N  
            Z  = AAD(I,N1)  
            AAD(I,N1) = Q*Z + P*AAD(I,N)  
            AAD(I,N)  = Q*AAD(I,N) - P*Z  
          ENDDO
!                                          ** Accumulate transformations  
          DO I = L, K  
            Z  = EVECD(I,N1)  
            EVECD(I,N1) = Q*Z + P*EVECD(I,N)  
            EVECD(I,N)  = Q*EVECD(I,N) - P*Z  
          ENDDO
          N  = N2  
          Reset = .true._lk
          CYCLE Eigen_loop
        END IF  
  
        IF( IN == 30_ik ) THEN  
!                    ** No convergence after 30 iterations; set error  
!                       indicator to the index of the current eigenvalue  
          IER  = N  
          GO TO  700  
        END IF  
!                                                  ** Form shift  
        IF( IN == 10_ik .OR. IN == 20_ik ) THEN  
          T  = T + X  
  
          DO I = L, N  
            AAD(I,I) = AAD(I,I) - X  
          ENDDO
  
          S  = ABS( AAD(N,N1) ) + ABS( AAD(N1,N2) )  
          X  = C3*S  
          Y  = X  
          W  = -C1*S**2  
        END IF  
  
        IN  = IN + iONE
  
!                ** Look for two consecutive small sub-diagonal elements  
        DO J = LB, N2  
          I  = N2 + LB - J  
          Z  = AAD( I, I )  
          R  = X - Z  
          S  = Y - Z  
          P  = ( R*S - W ) / AAD( I + iONE, I ) + AAD( I, I + iONE )  
          Q  = AAD( I + iONE, I + iONE ) - Z - R - S  
          R  = AAD( I + iTWO, I + iONE )  
          S  = ABS( P ) + ABS( Q ) + ABS( R )  
          P  = P / S  
          Q  = Q / S  
          R  = R / S  
  
          IF( I.EQ.LB ) EXIT
  
          UU   = ABS( AAD( I, I-iONE ) )*( ABS( Q ) + ABS( R ) )  
          VV   = ABS( P ) * ( ABS( AAD( I-iONE, I-iONE ) ) + ABS( Z ) + &
                             ABS( AAD( I+iONE, I+iONE ) ) )  
  
          IF( UU .LE. TOL*VV ) EXIT
        ENDDO
  
        AAD( I+iTWO, I ) = rZERO  
  
!                      ** fpp vectorization of this loop triggers  
!                         array bounds errors, so inhibit  
        DO J = I + 3_ik, N  
          AAD( J, J - iTWO ) = rZERO  
          AAD( J, J - 3_ik ) = rZERO  
        ENDDO
  
!             ** Double QR step involving rows K to N and columns M to N  
        DO KA = I, N1  
          NOTLAS = KA /= N1  
          IF( KA == I ) THEN  
            S  = SIGN( SQRT( P*P + Q*Q + R*R ), P )  
            IF( LB /= I ) AAD(KA,KA - iONE) = -AAD(KA,KA - iONE )  
          ELSE  
            P  = AAD(KA,KA - iONE)  
            Q  = AAD(KA + iONE,KA - iONE)  
            R  = rZERO  
            IF( NOTLAS ) R  = AAD(KA + iTWO,KA - iONE)  
  
            X  = ABS( P ) + ABS( Q ) + ABS( R )  
  
            IF( X == rZERO ) CYCLE
  
            P  = P / X  
            Q  = Q / X  
            R  = R / X  
            S  = SIGN( SQRT( P*P + Q*Q + R*R ), P )  
            AAD(KA,KA - iONE) = -S*X  
          END IF  
  
          P  = P + S  
          X  = P / S  
          Y  = Q / S  
          Z  = R / S  
          Q  = Q / P  
          R  = R / P  
!                                              ** Row modification  
          DO J = KA, M  
            P  = AAD(KA,J) + Q*AAD(KA + iONE,J)  
            IF( NOTLAS ) THEN  
               P  = P + R*AAD(KA + iTWO,J)  
               AAD(KA + iTWO,J) = AAD(KA + iTWO,J) - P*Z  
            END IF  
            AAD(KA + iONE,J ) = AAD(KA + iONE,J) - P*Y  
            AAD(KA,J) = AAD(KA,J) - P*X  
          ENDDO
!                                                 ** Column modification  
          DO II = iONE, MIN( N, KA + 3_ik )  
            P  = X*AAD(II,KA) + Y*AAD(II,KA + iONE)  
            IF( NOTLAS ) THEN  
               P  = P + Z*AAD(II,KA + iTWO)  
               AAD(II,KA + iTWO) = AAD(II,KA + iTWO) - P*R  
            END IF  
            AAD(II,KA + iONE) = AAD(II,KA + iONE) - P*Q  
            AAD(II,KA) = AAD(II,KA) - P  
          ENDDO
!                                          ** Accumulate transformations  
          DO II = L, K  
            P  = X*EVECD(II,KA) + Y*EVECD(II,KA + iONE)  
            IF( NOTLAS ) THEN  
               P  = P + Z*EVECD(II,KA + 2_ik)  
               EVECD(II,KA + 2_ik) = EVECD(II,KA + 2_ik) - P*R  
            END IF  
            EVECD(II,KA + iONE) = EVECD(II,KA + iONE) - P*Q  
            EVECD(II,KA) = EVECD(II,KA) - P  
          ENDDO
        ENDDO
          Reset = .false._ik
        ELSE
          EXIT Eigen_loop
        ENDIF
      ENDDO Eigen_loop
  
!                     ** All evals found, now backsubstitute real vector  
      IF( RNORM /= rZERO ) THEN  
         DO N = M, iONE, -iONE  
            N2   = N  
            AAD(N,N) = ONE  
  
            DO I = N - iONE, iONE, -iONE  
               W  = AAD(I,I) - EVALD(N)  
               IF( W == ZERO ) W  = TOL*RNORM  
               R  = AAD(I,N)  
               DO J = N2, N - iONE
                  R  = R + AAD(I,J)*AAD(J,N)  
               ENDDO
               AAD(I,N) = -R / W  
               N2   = I  
            ENDDO
         ENDDO
!                      ** End backsubstitution vectors of isolated evals  
         DO I = iONE, M  
            IF( I < L .OR. I > K ) THEN  
               DO J = I, M  
                  EVECD(I,J) = AAD(I,J)  
               ENDDO
            END IF  
         ENDDO
!                                   ** Multiply by transformation matrix  
         IF( K /= iZERO ) THEN  
            DO J = M, L, -iONE
               DO I = L, K  
                  Z  = rZERO  
                  DO N = L, MIN( J, K )  
                     Z  = Z + EVECD(I,N)*AAD(N,J)  
                  ENDDO
                  EVECD(I,J) = Z  
               ENDDO
            ENDDO
         ENDIF  
      ENDIF  
  
  
      DO I = L, K  
         DO J = iONE, M  
            EVECD(I,J) = EVECD(I,J)*WKD(I)  
         ENDDO
      ENDDO
  
!                           ** Interchange rows if permutations occurred  
      DO I = L-iONE, iONE, -iONE  
         J  = WKD(I)  
         IF( I /= J ) THEN  
            DO N = iONE, M  
               REPL   = EVECD(I,N)  
               EVECD(I,N) = EVECD(J,N)  
               EVECD(J,N) = REPL  
            ENDDO
         ENDIF  
      ENDDO
  
      DO I = K + iONE, M  
         J  = WKD(I)  
         IF( I /= J ) THEN  
            DO N = iONE, M  
               REPL   = EVECD(I,N)  
               EVECD(I,N) = EVECD(J,N)  
               EVECD(J,N) = REPL  
            ENDDO
         END IF  
      ENDDO
  
!                         ** Put results into output arrays  
  700 CONTINUE  
  
      DO J = iONE, M  
         EVAL(J) = EVALD(J)  
         DO K = iONE, M  
            EVEC(J,K) = EVECD(J,K)  
         ENDDO
      ENDDO
  
      END SUBROUTINE ASYMTX  
  
      SUBROUTINE CHEKIN( NLYR, DTAUC, SSALB, PMOM, &
                         NTAU, UTAU, NSTR, NUMU, &
                         UMU, NPHI, PHI, UMU0, &
                         FISOT, ALBEDO, HL, TAUC )  
  
!           Checks the input dimensions and variables  
  
!   Calls- WRTBAD, WRTDIM, DREF, ERRMSG  
!   Called by- DISORT  
! --------------------------------------------------------------------  
  
!     .. Scalar Arguments ..  
  
      integer(ik), intent(in) ::   NLYR, NPHI, NSTR, NTAU, NUMU  
      real(dk), intent(in)    ::   ALBEDO, FISOT, UMU0  
!     ..  
!     .. Array Arguments ..  
  
      real(dk), intent(in) ::  
                DTAUC(:), HL(0:), PHI(:),  
                PMOM(0:,:), SSALB(:),  
                TAUC(0:), UMU(:)  
      real(dk), intent(inout) :: UTAU(:)  
!     ..  
!     .. Local Scalars ..  
  
      logical(lk)  :: INPERR  
      integer(ik)  :: IRMU, IU, J, K, LC, LU  
      real(dk)     :: FLXALB, RMU  
  
!     ..  
  
      INPERR = .FALSE.  
  
      IF( NLYR < 1 ) INPERR = WRTBAD( 'NLYR' )  
  
      IF( NLYR > MXCLY ) INPERR = WRTBAD( 'MAXCLY' )  
  
      IF( ANY( DTAUC(:) < rZERO ) ) THEN  
         INPERR = WRTBAD( 'DTAUC' )  
      ENDIF  
      IF( ANY( SSALB(:) < rZERO ) .or. ANY( SSALB(:) > rONE ) ) THEN  
         INPERR = WRTBAD( 'SSALB' )  
      ENDIF  
      DO LC = 1, NLYR  
         IF( ANY( PMOM(:,LC) < -rONE ) .or.   
             ANY( PMOM(:,LC) > rONE ) ) THEN  
            INPERR = WRTBAD( 'PMOM' )  
         ENDIF  
      ENDDO  
  
      IF( MXULV < NLYR + iONE ) INPERR = WRTBAD( 'MAXULV' )  
  
      IF( NSTR < 2_ik .OR. MOD(NSTR,2_ik) /= 0 ) INPERR = WRTBAD( 'NSTR' )  
  
  
      IF( NSTR > MXCMU ) INPERR = WRTBAD( 'MAXCMU' )  
  
      IF( USRANG ) THEN  
         IF( NUMU < 0 ) INPERR = WRTBAD( 'NUMU' )  
  
         IF( .NOT. ONLYFL .AND. NUMU == iZERO ) INPERR = WRTBAD( 'NUMU' )  
  
         IF( NUMU > MXUMU ) INPERR = WRTBAD( 'MXUMU' )  
  
         IF( IBCND == iONE .AND. 2_ik*NUMU > MXUMU )  
             INPERR = WRTBAD( 'MXUMU' )  
  
         DO IU = iONE, NUMU  
            IF( UMU(IU) < -rONE .OR. UMU(IU) > rONE .OR. &
                UMU(IU) == rZERO ) INPERR = WRTBAD( 'UMU' )  
  
            IF( IBCND == iONE .AND. UMU(IU) < rZERO )  INPERR = WRTBAD( 'UMU' )  
  
            IF( IU > iONE ) THEN  
               IF( UMU(IU) < UMU(IU-1) ) INPERR = WRTBAD( 'UMU' )  
            END IF  
         ENDDO  
      ELSE  
         IF( MXUMU < NSTR ) INPERR = WRTBAD( 'MAXUMU' )  
      END IF  
  
  
      IF( .NOT. ONLYFL .AND. IBCND /= iONE ) THEN  
         IF( NPHI <= 0 ) INPERR = WRTBAD( 'NPHI' )  
  
         IF( NPHI > MXPHI ) INPERR = WRTBAD( 'MAXPHI' )  
  
         IF( ANY( PHI(:) < rZERO ) .OR. ANY( PHI(:) > 360.0_dk ) )  INPERR = WRTBAD( 'PHI' )  
      END IF  
  
  
      IF( IBCND.LT.iZERO .OR. IBCND.GT.iONE ) INPERR = WRTBAD( 'IBCND' )  
  
      IF( IBCND == iZERO ) THEN  
         IF( FBEAM.LT.rZERO ) INPERR = WRTBAD( 'FBEAM' )  
  
         IF( FBEAM.GT.rZERO .AND. abs(UMU0).GT.rONE ) INPERR = WRTBAD( 'UMU0' )  
  
         IF( FBEAM.GT.rZERO .AND. (PHI0.LT.rZERO .OR.PHI0.GT.360.0_dk) ) INPERR = WRTBAD( 'PHI0' )  
  
         IF( FISOT.LT.rZERO ) INPERR = WRTBAD( 'FISOT' )  
!                    ** Make sure flux albedo at dense mesh of incident  
!                       angles does not assume unphysical values  
         IF( (.NOT. ONLYFL .AND. USRANG) .OR. .NOT. LAMBER ) THEN  
           DO IRMU = iZERO, 100_ik
             RMU  = real(IRMU,dk)*0.01_dk
             FLXALB = DREF(RMU,HL,NSTR)  
             IF( FLXALB < rZERO .OR. FLXALB > rONE )  
                INPERR = WRTBAD( 'HL' )  
           ENDDO  
         ENDIF  
      ELSE IF( IBCND == 1 ) THEN  
         IF( ALBEDO < rZERO .OR. ALBEDO > rONE ) INPERR = WRTBAD( 'ALBEDO' )  
      END IF  
  
      IF( ACCUR < rZERO .OR. ACCUR > 1.E-2_dk ) INPERR = WRTBAD( 'ACCUR' )  
  
      IF( MXCLY.LT.NLYR ) INPERR = WRTDIM( 'MXCLY', NLYR )  
  
      IF( IBCND /= 1 ) THEN  
         IF( USRTAU .AND. MXULV.LT.NTAU ) INPERR = WRTDIM( 'MXULV',NTAU )  
         IF( .NOT.USRTAU .AND. MXULV .LT. NLYR + iONE ) THEN
             INPERR = WRTDIM( 'MXULV', NLYR + iONE )  
         ENDIF
      ELSE  
         IF( MXULV.LT.2_ik ) INPERR = WRTDIM( 'MXULV', 2_ik )  
      END IF  
  
      IF( MXCMU.LT.NSTR ) INPERR = WRTDIM( 'MXCMU', NSTR )  
  
      IF( USRANG .AND. MXUMU.LT.NUMU ) INPERR = WRTDIM( 'MXUMU', NUMU )  
  
      IF( USRANG .AND. IBCND.EQ.iONE .AND.MXUMU.LT.2_ik*NUMU ) THEN
        INPERR = WRTDIM( 'MXUMU', NUMU )  
      ENDIF
  
      IF( .NOT.USRANG .AND. MXUMU.LT.NSTR ) INPERR = WRTDIM( 'MXUMU', NSTR )  
  
      IF( .NOT.ONLYFL .AND. IBCND.NE.iONE .AND. MXPHI.LT.NPHI ) THEN
          INPERR = WRTDIM( 'MXPHI', NPHI )  
      ENDIF
  
      IF( INPERR )  
          CALL ERRMSG( 'DISORT--input and/or dimension errors',.True.)  
  
      END SUBROUTINE CHEKIN  
  
      SUBROUTINE FLUXES( tausla, tauslau, &
                         CMU, CWT, FBEAM, GC, KK, LAYRU, LL, LYRCUT, &
                         NCUT, NN, NSTR, NTAU, PI, &
                         PRNT, SSALB, TAUCPR, UMU0, UTAU, UTAUPR, XR0, &
                         XR1, ZZ, ZPLK0, ZPLK1, DFDT, FLUP, FLDN, FLDIR, &
                         RFLDIR, RFLDN, UAVG, U0C, &
                         uavgso, uavgup, uavgdn, sindir, sinup, sindn)  
  
!       Calculates the radiative fluxes, mean intensity, and flux  
!       derivative with respect to optical depth from the m=0 intensity  
!       components (the azimuthally-averaged intensity)  
  
!    I N P U T     V A R I A B L E S:  
  
!       CMU      :  Abscissae for Gauss quadrature over angle cosine  
!       CWT      :  Weights for Gauss quadrature over angle cosine  
!       GC       :  Eigenvectors at polar quadrature angles, SC(1)  
!       KK       :  Eigenvalues of coeff. matrix in Eq. SS(7)  
!       LAYRU    :  Layer number of user level UTAU  
!       LL       :  Constants of integration in Eq. SC(1), obtained  
!                     by solving scaled version of Eq. SC(5);  
!                     exponential term of Eq. SC(12) not included  
!       LYRCUT   :  Logical flag for truncation of comput. layer  
!       NN       :  Order of double-Gauss quadrature (NSTR/2)  
!       NCUT     :  Number of computational layer where absorption  
!                     optical depth exceeds ABSCUT  
!       TAUCPR   :  Cumulative optical depth (delta-M-scaled)  
!       UTAUPR   :  Optical depths of user output levels in delta-M  
!                     coordinates;  equal to UTAU if no delta-M  
!       XR0      :  Expansion of thermal source function in Eq. SS(14)  
!       XR1      :  Expansion of thermal source function Eqs. SS(16)  
!       ZZ       :  Beam source vectors in Eq. SS(19)  
!       ZPLK0    :  Thermal source vectors Z0, by solving Eq. SS(16)  
!       ZPLK1    :  Thermal source vectors Z1, by solving Eq. SS(16)  
!       (remainder are DISORT input variables)  
  
  
!                   O U T P U T     V A R I A B L E S:  
  
!       U0C      :  Azimuthally averaged intensities  
!                   ( at polar quadrature angles )  
!       (RFLDIR, RFLDN, FLUP, DFDT, UAVG are DISORT output variables)  
  
  
!                   I N T E R N A L       V A R I A B L E S:  
  
!       DIRINT   :  Direct intensity attenuated  
!       FDNTOT   :  Total downward flux (direct + diffuse)  
!       FLDIR    :  Direct-beam flux (delta-M scaled)  
!       FLDN     :  Diffuse down-flux (delta-M scaled)  
!       FNET     :  Net flux (total-down - diffuse-up)  
!       FACT     :  EXP( - UTAUPR / UMU0 )  
!       PLSORC   :  Planck source function (thermal)  
!       ZINT     :  Intensity of m = 0 case, in Eq. SC(1)  
  
!   Called by- DISORT  
! +-------------------------------------------------------------------+  
  
!     .. Scalar Arguments ..  
  
      integer(ik), intent(in) :: NCUT, NN, NSTR, NTAU  
      real(dk), intent(in)    :: FBEAM, PI, UMU0  
      logical(lk), intent(in) :: LYRCUT  
!     ..  
!     .. Array Arguments ..  
  
      logical(lk), intent(in) :: PRNT(:)  
      integer(ik), intent(in) :: LAYRU(:)  
      real(dk),    intent(in) :: CMU(:), CWT(:),  
     &          GC(:,:,:), KK(:,:), LL(:,:),  
     &          SSALB(:),  TAUCPR(0:),  
     &          UTAU(:), UTAUPR(:), XR0(:), XR1(:),  
     &          ZPLK0(:,:), ZPLK1(:,:), ZZ(:,:)  
      real(dk), intent(in)  :: tausla(0:), tauslau(0:)  
      real(dk), intent(out) :: U0C(:,:)  
      real(dk), intent(out) :: RFLDIR(:), RFLDN(:), FLUP(:)  
      real(dk), intent(out) :: DFDT(:), UAVG(:)  
      real(dk), intent(out) :: uavgso(:), uavgup(:), uavgdn(:)  
      real(dk), intent(out) :: sindir(:), sinup(:), sindn(:)  
      real(dk), intent(out) :: FLDIR(:), FLDN(:)  
!     ..  
!     .. Local Scalars ..  
  
      integer(ik) :: IQ, JQ, LU, LYU  
      real(dk)    :: ANG1, ANG2, DIRINT, FACT, FDNTOT, FNET, PLSORC, ZINT  
!     ..  
  
      IF( PRNT( 2 ) ) WRITE( *, 9000 )  
!                                          ** Zero DISORT output arrays  
      U0C   = rZERO
      FLDIR = rZERO
      FLDN  = rZERO
      uavgso = rZERO
      uavgup = rZERO
      uavgdn = rZERO
      sindir = rZERO
      sinup  = rZERO
      sindn  = rZERO
  
!    ** Loop over user levels  
      LEVEL_LOOP: DO LU = iONE, NTAU  
         LYU  = LAYRU( LU )  
         IF( LYRCUT .AND. LYU > NCUT ) THEN  
!                                                ** No radiation reaches  
!                                                ** this level  
            FDNTOT = rZERO  
            FNET   = rZERO  
            PLSORC = rZERO  
            IF( PRNT( 2 ) ) WRITE( *, FMT = 9010 ) UTAU( LU ), LYU, &
              RFLDIR( LU ), RFLDN( LU ), FDNTOT, FLUP( LU ), FNET, & 
              UAVG( LU ), PLSORC, DFDT( LU )  
            CYCLE LEVEL_LOOP
         END IF  
  
         IF( FBEAM > rZERO ) THEN  
            FACT  = EXP( - tausla(LU-1) )  
            DIRINT       = FBEAM*FACT  
            FLDIR( LU )  = UMU0*( FBEAM*FACT )  
            RFLDIR( LU ) = UMU0*FBEAM * EXP( -tauslau(lu-1) )  
            sindir( LU ) = SQRT(rONE - UMU0*UMU0)*FBEAM * EXP( -tauslau(lu-1) )
         ELSE  
            DIRINT       = rZERO  
            FLDIR( LU )  = rZERO  
            RFLDIR( LU ) = rZERO  
            sindir( LU ) = rZERO  
         END IF  
  
  
         DO IQ = iONE, NN  
            ZINT   = rZERO  
            DO JQ = iONE, NN  
               ZINT   = ZINT + GC(IQ,JQ,LYU)*LL(JQ,LYU) &
                             * EXP( -KK(JQ,LYU)*(UTAUPR(LU) - TAUCPR(LYU)) )
            ENDDO  
  
            DO JQ = NN + iONE, NSTR  
               ZINT   = ZINT + GC(IQ,JQ,LYU)*LL(JQ,LYU) &
                        * EXP( -KK(JQ,LYU)*(UTAUPR(LU) - TAUCPR(LYU - iONE)) )  
            ENDDO  
  
            U0C(IQ,LU) = ZINT  
  
            IF( FBEAM > rZERO ) THEN  
              U0C(IQ,LU) = ZINT + ZZ(IQ,LYU)*FACT  
            ENDIF  
  
            U0C(IQ,LU) = U0C(IQ,LU) + ZPLK0(IQ,LYU) + ZPLK1(IQ,LYU)*UTAUPR(LU)  
            UAVG(LU)   = UAVG(LU) + CWT(NN + iONE - IQ)*U0C(IQ,LU)  
            uavgdn(lu) = uavgdn(lu) + cwt(nn+iONE-iq) * u0c(iq,lu)  
            sindn(lu)  = sindn(lu)  + cwt(nn+iONE-iq) * SQRT(rONE - CMU(NN+iONE-IQ)*CMU(NN+iONE-IQ)) * U0C(IQ,LU)  
            FLDN(LU)   = FLDN(LU) + CWT(NN + iONE - IQ) * CMU(NN + iONE - IQ)*U0C(IQ,LU)  
         ENDDO  
  
  
         DO IQ = NN + iONE, NSTR  
            ZINT   = rZERO  
            DO JQ = 1, NN  
               ZINT   = ZINT + GC(IQ,JQ,LYU)*LL(JQ,LYU) &
                        * EXP( -KK(JQ,LYU)*(UTAUPR(LU) - TAUCPR(LYU)) )  
            ENDDO  
  
            DO JQ = NN + 1, NSTR  
               ZINT   = ZINT + GC(IQ,JQ,LYU)*LL(JQ,LYU ) &
                        * EXP( -KK(JQ,LYU)*(UTAUPR(LU) - TAUCPR(LYU - iONE)) )
            ENDDO  
  
            U0C( IQ, LU ) = ZINT  
  
            IF( FBEAM > rZERO ) THEN  
              U0C(IQ,LU) = ZINT + ZZ(IQ,LYU)*FACT  
            ENDIF  
  
            U0C(IQ,LU) = U0C(IQ,LU) + ZPLK0(IQ,LYU) + ZPLK1(IQ,LYU)*UTAUPR(LU)  
            UAVG(LU)   = UAVG(LU) + CWT(IQ - NN)*U0C(IQ,LU)  
            uavgup(lu) = uavgup(lu) + cwt(iq-nn) * u0c(iq,lu)  
            sinup (lu) = sinup(lu)  + cwt(iq-nn) * SQRT(rONE-CMU(IQ-NN)*CMU(IQ-NN))*U0C(IQ,LU)  
            FLUP(LU)   = FLUP(LU) + CWT(IQ - NN)*CMU(IQ - NN) * U0C(IQ,LU)  
         ENDDO  
  
  
         FLUP(LU)  = TWOPI*FLUP(LU)  
         FLDN(LU)  = TWOPI*FLDN(LU)  
         FDNTOT    = FLDN(LU) + FLDIR(LU)  
         FNET      = FDNTOT - FLUP(LU)  
         RFLDN(LU) = FDNTOT - RFLDIR(LU)  
         UAVG(LU)  = (TWOPI*UAVG(LU) + DIRINT) / FOURPI
         uavgso( lu ) = dirint / FOURPI
         uavgup( lu ) = (TWOPI * uavgup(lu) ) / FOURPI
         uavgdn( lu)  = (TWOPI * uavgdn(lu) ) / FOURPI
         sindn ( lu ) = TWOPI*sindn(LU)  
         sinup ( lu ) = TWOPI*sinup(LU)  
  
         PLSORC    = XR0(LYU) + XR1(LYU)*UTAUPR(LU)  
         DFDT(LU)  = (rONE - SSALB(LYU)) * FOURPI * (UAVG(LU) - PLSORC)
  
         IF( PRNT( 2 ) ) WRITE( *, FMT = 9010 ) UTAU( LU ), LYU, &
             RFLDIR( LU ), RFLDN( LU ), FDNTOT, FLUP( LU ), FNET, & 
             UAVG( LU ), PLSORC, DFDT( LU )  
      ENDDO LEVEL_LOOP  
  
  
      IF( PRNT( 3 ) ) THEN  
         WRITE( *, FMT = 9020 )  
         DO LU = iONE, NTAU  
            WRITE( *, FMT = 9030 ) UTAU( LU )  
            DO IQ = iONE, NN  
               ANG1   = 180._dk/ PI* ACOS( CMU(2_ik*NN - IQ + iONE) )  
               ANG2   = 180._dk/ PI* ACOS( CMU(IQ) )  
               WRITE( *, 9040 ) ANG1, CMU(2_ik*NN-IQ+iONE), U0C(IQ,LU), &
                                ANG2, CMU(IQ),        U0C(IQ+NN,LU)  
            ENDDO  
         ENDDO  
      END IF  
  
  
 9000 FORMAT( //, 21X,  
       '<----------------------- FLUXES ----------------------->', /,  
       '   Optical  Compu    Downward    Downward    Downward     ',  
       ' Upward                    Mean      Planck   d(Net Flux)', /,  
       '     Depth  Layer      Direct     Diffuse       Total     ',  
       'Diffuse         Net   Intensity      Source   / d(Op Dep)', / )  
 9010 FORMAT( F10.4, I7, 1P, 7E12.3, E14.3 )  
 9020 FORMAT( / , / , ' ******** AZIMUTHALLY AVERAGED INTENSITIES',  
     &      ' ( at polar quadrature angles ) *******' )  
 9030 FORMAT( /, ' Optical depth =', F10.4, //,  
        '     Angle (deg)   cos(Angle)     Intensity',  
        '     Angle (deg)   cos(Angle)     Intensity' )  
 9040 FORMAT( 2( 0P,F16.4,F13.5,1P,E14.3 ) )  
  
      END SUBROUTINE FLUXES  
  
      SUBROUTINE LEPOLY( NMU, M, TWONM1, MU, YLM )  
  
!       Computes the normalized associated Legendre polynomial,  
!       defined in terms of the associated Legendre polynomial  
!       Plm = P-sub-l-super-m as  
  
!             Ylm(MU) = sqrt( (l-m)!/(l+m)! ) * Plm(MU)  
  
!       for fixed order m and all degrees from l = m to TWONM1.  
!       When m.GT.0, assumes that Y-sub(m-1)-super(m-1) is available  
!       from a prior call to the routine.  
  
!       REFERENCE: Dave, J.V. and B.H. Armstrong, Computations of  
!                  High-Order Associated Legendre Polynomials,  
!                  J. Quant. Spectrosc. Radiat. Transfer 10,  
!                  557-562, 1970.  (hereafter D/A)  
  
!       METHOD: Varying degree recurrence relationship.  
  
!       NOTE 1: The D/A formulas are transformed by  
!               setting  M = n-1; L = k-1.  
!       NOTE 2: Assumes that routine is called first with  M = 0,  
!               then with  M = 1, etc. up to  M = TWONM1.  
!       NOTE 3: Loops are written in such a way as to vectorize.  
  
!  I N P U T     V A R I A B L E S:  
  
!       NMU    :  Number of arguments of YLM  
!       M      :  Order of YLM  
!       MAXMU  :  First dimension of YLM  
!       TWONM1 :  Max degree of YLM  
!       MU(i)  :  Arguments of YLM (i = 1 to NMU)  
  
!       If M.GT.0, YLM(M-1,i) for i = 1 to NMU is assumed to exist  
!       from a prior call.  
  
!  O U T P U T     V A R I A B L E:  
  
!       YLM(l,i) :  l = M to TWONM1, normalized associated Legendre  
!                   polynomials evaluated at argument MU(i)  
  
!   Called by- DISORT, ALBTRN, SURFAC  
!   Calls- ERRMSG  
! +-------------------------------------------------------------------+  
  
!     .. Parameters ..  
  
      integer(ik), PARAMETER ::  MAXSQT = 1000_ik
!     ..  
!     .. Scalar Arguments ..  
  
      integer(ik), intent(in) :: M, NMU, TWONM1  
!     ..  
!     .. Array Arguments ..  
  
      real(dk), intent(in)  :: MU(:)  
      real(dk), intent(out) :: YLM(0:,:)  
!     ..  
!     .. Local Scalars ..  
  
      logical(lk)   :: PASS1  
      integer(ik)   :: I, L, NS  
      real(dk)      :: TMP1, TMP2  
!     ..  
!     .. Local Arrays ..  
  
      real(dk) ::  SQT( MAXSQT )  
!     ..  
!     .. External Subroutines ..  
  
      SAVE      SQT, PASS1  
      DATA      PASS1 / .TRUE._lk /  
  
  
      IF( PASS1 ) THEN  
         PASS1  = .FALSE._lk  
         DO NS = 1, MAXSQT  
            SQT( NS ) = SQRT( real(dk)( NS ) )  
         ENDDO  
      END IF  
  
      IF( iTWO*TWONM1 > MAXSQT )  
     &    CALL ERRMSG('LEPOLY--need to increase param MAXSQT',.True.)  
  
  
      IF( M == iZERO ) THEN  
!                             ** Upward recurrence for ordinary  
!                                Legendre polynomials  
         DO I = iONE, NMU  
            YLM(iZERO,I) = rONE  
            YLM(iONE,I)  = MU(I)  
         ENDDO  
  
         DO L = iTWO, TWONM1  
            DO I = iONE, NMU  
               YLM(L,I) = ((iTWO*L - iONE)*MU(I)*YLM(L - iONE,I) - (L - iONE)*YLM(L - iTWO,I)) / L  
            ENDDO  
         ENDDO  
      ELSE  
         DO I = iONE, NMU  
!                               ** Y-sub-m-super-m; derived from  
!                               ** D/A Eqs. (11,12)  
            YLM(M,I) = -SQT(real(iTWO*M - iONE),dk) / SQT(real(iTWO*M,dk)) &
                       * SQRT(rONE - MU(I)**2)*YLM(M - iONE,I)  
  
!                              ** Y-sub-(m+1)-super-m; derived from  
!                              ** D/A Eqs.(13,14) using Eqs.(11,12)  
            YLM(M + iONE,I) = SQT(real(iTWO*M + iONE,dk))*MU(I)*YLM(M,I)  
         ENDDO  
!                                   ** Upward recurrence; D/A EQ.(10)  
         DO L = M + iTWO, TWONM1  
            TMP1 = SQT(L - M )*SQT(L + M)  
            TMP2 = SQT(L - M - iONE)*SQT(L + M - iONE)  
            DO I = iONE, NMU  
               YLM(L,I) = ((real(iTWO*L - 1,dk))*MU(I)*YLM(L-iONE,I) - TMP2*YLM(L-iTWO,I)) / TMP1  
            ENDDO  
         ENDDO  
      END IF  
  
      END SUBROUTINE LEPOLY  
  
      SUBROUTINE PRAVIN( UMU, NUMU, UTAU, NTAU, U0U )  
  
!        Print azimuthally averaged intensities at user angles  
  
!   Called by- DISORT  
  
!     LENFMT   Max number of polar angle cosines UMU that can be  
!                printed on one line, as set in FORMAT statement  
! --------------------------------------------------------------------  
  
!     .. Scalar Arguments ..  
  
      integer(ik)   NTAU, NUMU  
!     ..  
!     .. Array Arguments ..  
  
      real(dk)      U0U( :, : ), UMU( : ), UTAU( : )  
!     ..  
!     .. Local Scalars ..  
  
      integer(ik)   IU, IUMAX, IUMIN, LENFMT, LU, NP, NPASS  
!     ..  
!     .. Intrinsic Functions ..  
  
      INTRINSIC MIN  
!     ..  
  
  
      IF( NUMU.LT.iONE )  RETURN  
  
      WRITE( *, '(//,A)' )  
     &   ' *******  AZIMUTHALLY AVERAGED INTENSITIES ' //  
     &   '(at user polar angles)  ********'  
  
      LENFMT = 8_ik  
      NPASS  = iONE + (NUMU-iONE) / LENFMT  
  
      WRITE( *,'(/,A,/,A)') '   Optical   Polar Angle Cosines Depth'  
  
      DO NP = iONE, NPASS  
         IUMIN  = iONE + LENFMT * ( NP - iONE )  
         IUMAX  = MIN( LENFMT*NP, NUMU )  
         WRITE( *,'(/,10X,8F14.5)') (UMU(IU), IU = IUMIN, IUMAX)
  
         DO LU = iONE, NTAU  
            WRITE( *, '(0P,F10.4,1P,8E14.4)' ) UTAU( LU ), (U0U(IU,LU), IU = IUMIN, IUMAX)  
         ENDDO
      ENDDO
  
      END SUBROUTINE PRAVIN 
  
      SUBROUTINE PRTINP( NLYR, DTAUC, DTAUCP, SSALB, PMOM, &
                         NTAU, UTAU, NSTR, NUMU, UMU, &
                         PHI, UMU0, FISOT, &
                         ALBEDO, HL, FLYR, LYRCUT, &
                         OPRIM, TAUC, TAUCPR, PRTMOM )  
  
!        Print values of input variables  
  
!   Called by- DISORT  
! --------------------------------------------------------------------  
  
!     .. Scalar Arguments ..  
  
      logical(lk), intent(in) :: LYRCUT, PRTMOM  
      integer(ik), intent(in) :: NLYR, NSTR, NTAU, NUMU  
      real(dk), intent(in)    :: ALBEDO, FISOT, UMU0  
!     ..  
!     .. Array Arguments ..  
  
      real(dk), intent(in) ::  DTAUC(:), DTAUCP(:), FLYR(:), HL(0:),  
                OPRIM(:), PHI(:), PMOM(0:,:), SSALB(:),  
                TAUC(0:), TAUCPR(0:), UMU(:), UTAU(:)  
!     ..  
!     .. Local Scalars ..  
  
      integer(ik) :: IU, J, K, LC, LU  
      real(dk)    :: YESSCT  
!     ..  
  
  
      WRITE( *, '(/,A,I4,A,I4)' ) ' No. streams =', NSTR,'     No. computational layers =', NLYR  
  
      IF( IBCND /= iONE ) WRITE( *, '(I4,A,10F10.4,/,(26X,10F10.4))' ) &
          NTAU,' User optical depths :', ( UTAU(LU), LU = 1, NTAU )  
  
      IF( .NOT. ONLYFL ) WRITE( *, '(I4,A,10F9.5,/,(31X,10F9.5))' ) &
          NUMU,' User polar angle cosines :',( UMU(IU), IU = 1, NUMU )  
  
      IF( .NOT. ONLYFL .AND. IBCND /= iONE ) &
          WRITE( *, '(I4,A,10F9.2,/,(28X,10F9.2))' ) &
                 NPHI,' User azimuthal angles :',( PHI(J), J = iONE, NPHI )  
  
      IF( .NOT. PLANK .OR. IBCND == iONE ) WRITE( *, '(A)' ) ' No thermal emission'  
  
  
      WRITE( *, '(A,I2)' ) ' Boundary condition flag: IBCND =', IBCND  
  
      IF( IBCND == iZERO ) THEN  
         WRITE( *, '(A,1P,E11.3,A,0P,F8.5,A,F7.2,/,A,1P,E11.3)' ) &
                '    Incident beam with intensity =', FBEAM, &
                ' and polar angle cosine = ', UMU0, &
                '  and azimuth angle =', PHI0, &
                '    plus isotropic incident intensity =', FISOT  
  
         IF( LAMBER ) WRITE( *, '(A,0P,F8.4)' ) '    Bottom albedo (Lambertian) =', ALBEDO  
  
         IF( .NOT. LAMBER ) WRITE( *, '(A,/,(10X,10F9.5))' ) &
           '    Legendre coeffs of bottom bidirectional reflectivity :', &
               (HL(K), K = iZERO, NSTR)  
      ELSE IF( IBCND == iONE ) THEN  
         WRITE(*,'(A)') '    Isotropic illumination from top and bottom'  
         WRITE( *, '(A,0P,F8.4)' ) '    Bottom albedo (Lambertian) =', ALBEDO  
      END IF  
  
  
      IF( DELTAM ) WRITE( *, '(A)' ) ' Uses delta-M method'  
      IF( .NOT.DELTAM ) WRITE( *, '(A)' ) ' Does not use delta-M method'  
  
  
      IF( IBCND == iONE ) THEN  
         WRITE( *, '(A)' ) ' Calculate albedo and transmissivity of'// &
                           ' medium vs. incident beam angle'  
      ELSE IF( ONLYFL ) THEN  
         WRITE( *, '(A)' ) &
                ' Calculate fluxes and azim-averaged intensities only'  
      ELSE  
         WRITE( *, '(A)' ) ' Calculate fluxes and intensities'  
      END IF  
  
      WRITE( *, '(A,1P,E11.2)' ) &
             ' Relative convergence criterion for azimuth series =', &
             ACCUR  
  
      IF( LYRCUT ) WRITE( *, '(A)' ) &
          ' Sets radiation = 0 below absorption optical depth 10'  
  
  
!                                        ** Print layer variables  
      IF( PLANK ) WRITE( *, FMT = 9180 )  
      IF( .NOT. PLANK ) WRITE( *, FMT = 9190 )  
  
      YESSCT = rZERO  
  
      DO LC = iONE, NLYR  
         YESSCT = YESSCT + SSALB(LC)  
         IF( PLANK ) THEN
             WRITE(*,'(I4,2F10.4,F10.5,F12.5,2F10.4,F10.5,F9.4,F14.3)') &
                   LC, DTAUC( LC ), TAUC( LC ), SSALB( LC ), FLYR( LC ),& 
                   DTAUCP( LC ), TAUCPR( LC ), OPRIM( LC ), PMOM(1,LC)  
         ELSE
             WRITE(*,'(I4,2F10.4,F10.5,F12.5,2F10.4,F10.5,F9.4)') &
                   LC, DTAUC( LC ), TAUC( LC ), SSALB( LC ), FLYR( LC ), &
                   DTAUCP( LC ), TAUCPR( LC ), OPRIM( LC ), PMOM( 1,LC )  
         ENDIF
      ENDDO  
  
  
      IF( PRTMOM .AND. YESSCT > rZERO ) THEN  
         WRITE( *, '(/,A)' ) ' Layer   Phase Function Moments'  
         DO LC = iONE, NLYR  
            IF( SSALB(LC).GT.rZERO ) &
                WRITE( *, '(I6,10F11.6,/,(6X,10F11.6))' ) &
                       LC, (PMOM(K,LC), K = iZERO,NSTR)  
         ENDDO  
      END IF  
  
!                ** (Read every other line in these formats)  
  
 9180 FORMAT( /, 37X, '<------------- Delta-M --------------->', /,  
     &'                   Total    Single                           ',  
     &               'Total    Single', /,  
     &'       Optical   Optical   Scatter   Truncated   ',  
     &   'Optical   Optical   Scatter    Asymm', /,  
     &'         Depth     Depth    Albedo    Fraction     ',  
     &     'Depth     Depth    Albedo   Factor   Temperature' )  
 9190 FORMAT( /, 37X, '<------------- Delta-M --------------->', /,  
     &'                   Total    Single                           ',  
     &               'Total    Single', /,  
     &'       Optical   Optical   Scatter   Truncated   ',  
     &   'Optical   Optical   Scatter    Asymm', /,  
     &'         Depth     Depth    Albedo    Fraction     ',  
     &     'Depth     Depth    Albedo   Factor' )  
  
      END SUBROUTINE PRTINP
  
      SUBROUTINE PRTINT( UU, UTAU, NTAU, UMU, NUMU, PHI )  
  
!         Prints the intensity at user polar and azimuthal angles  
  
!     All arguments are DISORT input or output variables  
  
!   Called by- DISORT  
  
!     LENFMT   Max number of azimuth angles PHI that can be printed  
!                on one line, as set in FORMAT statement  
! +-------------------------------------------------------------------+  
  
  
!     .. Scalar Arguments ..  
  
      integer(ik)   NPHI, NTAU, NUMU  
!     ..  
!     .. Array Arguments ..  
  
      real(dk)      PHI(:), UMU(:), UTAU(:), UU(:,:,:)  
!     ..  
!     .. Local Scalars ..  
  
      integer(ik)   IU, J, JMAX, JMIN, LENFMT, LU, NP, NPASS  
!     ..  
!     .. Intrinsic Functions ..  
  
      INTRINSIC MIN  
!     ..  
  
  
      IF( NPHI.LT.iONE )  RETURN  
  
      WRITE( *, '(//,A)' )  
     &   ' *********  I N T E N S I T I E S  *********'  
  
      LENFMT = 10_ik  
      NPASS  = iONE + (NPHI-iONE) / LENFMT  
  
      WRITE( *, '(/,A,/,A,/,A)' ) &
         '             Polar   Azimuth angles (degrees)', &
         '   Optical   Angle', &
         '    Depth   Cosine'  
  
      DO LU = iONE, NTAU  
         DO NP = iONE, NPASS  
            JMIN   = iONE + LENFMT * (NP - iONE)  
            JMAX   = MIN( LENFMT*NP, NPHI )  
  
            WRITE( *, '(/,18X,10F11.2)' ) ( PHI(J), J = JMIN, JMAX )  
  
            IF( NP.EQ.iONE ) WRITE( *, '(F10.4,F8.4,1P,10E11.3)' ) &
                   UTAU(LU), UMU(1), (UU(1, LU, J), J = JMIN, JMAX)  
            IF( NP.GT.iONE ) WRITE( *, '(10X,F8.4,1P,10E11.3)' ) &
                             UMU(1), (UU(1, LU, J), J = JMIN, JMAX)  
  
            DO IU = iTWO, NUMU  
               WRITE( *, '(10X,F8.4,1P,10E11.3)' ) &
                       UMU(IU), (UU(IU,LU,J), J = JMIN, JMAX)
            ENDDO
         ENDDO
      ENDDO
  
      END SUBROUTINE PRTINT  
  
      SUBROUTINE QGAUSN( M, GMU, GWT )  
  
!       Compute weights and abscissae for ordinary Gaussian quadrature  
!       on the interval (0,1);  that is, such that  
  
!           sum(i=1 to M) ( GWT(i) f(GMU(i)) )  
  
!       is a good approximation to  
  
!           integral(0 to 1) ( f(x) dx )  
  
!   INPUT :    M       order of quadrature rule  
  
!   OUTPUT :  GMU(I)   array of abscissae (I = 1 TO M)  
!             GWT(I)   array of weights (I = 1 TO M)  
  
!   REFERENCE:  Davis, P.J. and P. Rabinowitz, Methods of Numerical  
!                   Integration, Academic Press, New York, pp. 87, 1975  
  
!   METHOD:  Compute the abscissae as roots of the Legendre  
!            polynomial P-sub-M using a cubically convergent  
!            refinement of Newton's method.  Compute the  
!            weights from EQ. 2.7.3.8 of Davis/Rabinowitz.  Note  
!            that Newton's method can very easily diverge; only a  
!            very good initial guess can guarantee convergence.  
!            The initial guess used here has never led to divergence  
!            even for M up to 1000.  
  
!   ACCURACY:  relative error no better than TOL or computer  
!              precision (machine epsilon), whichever is larger  
  
!   INTERNAL VARIABLES:  
  
!    ITER      : number of Newton Method iterations  
!    MAXIT     : maximum allowed iterations of Newton Method  
!    PM2,PM1,P : 3 successive Legendre polynomials  
!    PPR       : derivative of Legendre polynomial  
!    P2PRI     : 2nd derivative of Legendre polynomial  
!    TOL       : convergence criterion for Legendre poly root iteration  
!    X,XI      : successive iterates in cubically-convergent version  
!                of Newtons Method (seeking roots of Legendre poly.)  
  
!   Called by- SETDIS, SURFAC  
!   Calls- D1MACH, ERRMSG  
! +-------------------------------------------------------------------+  
  
!     .. Scalar Arguments ..  
  
      integer(ik), intent(in) :: M  
!     ..  
!     .. Array Arguments ..  
  
      real(dk), intent(out)  :: GMU(:), GWT(:)  
!     ..  
!     .. Local Scalars ..  
  
      integer(ik), PARAMETER ::  MAXIT = 1000_ik
      real(dk)(8), PARAMETER ::  ONE = 1.D0, TWO = 2.D0  
  
      integer(ik) ::  ITER, K, LIM, NN, NP1  
      real(dk)    ::  CONA, PI, T  
      real(dk)    ::  EN, NNP1, P, P2PRI, PM1, PM2, PPR, PROD, &
                      TMP, TOL, X, XI  
  
      SAVE      PI, TOL  
      DATA      PI / rZERO /   
  
      IF( PI == rZERO ) THEN  
         PI   = rTWO*ASIN( rONE )  
!        TOL  = 10._dk*D1MACH( 4 )  
         TOL  = 10._dk*EPSILON( rONE )  
      END IF  
  
      IF( M < iONE ) CALL ERRMSG( 'QGAUSN--Bad value of M',.True.)  
  
      IF( M == iONE ) THEN  
         GMU(iONE) = 0.5_dk
         GWT(iONE) = rONE  
         RETURN  
      END IF  
  
      EN   = M  
      NP1  = M + iONE  
      NNP1 = M*NP1  
      CONA = real(M - 1,dk) / real(8_ik*M**3,dk)  
  
      LIM  = M / iTWO
  
      DO K = iONE, LIM  
!                                        ** Initial guess for k-th root  
!                                           of Legendre polynomial, from  
!                                           Davis/Rabinowitz (2.7.3.3a)  
         T  = (4_ik*K - iONE)*PI / real( 4_ik*M + iTWO,dk )  
         X  = COS( T + CONA / TAN( T ) )  
         ITER = iZERO
!                                        ** Upward recurrence for  
!                                           Legendre polynomials  
   10    CONTINUE  
         ITER   = ITER + iONE  
         PM2    = rONE  
         PM1    = X  
  
         DO NN = iTWO, M  
            P    = (real((iTWO*NN - iONE),dk)*X*PM1 - real(NN - iONE,dk)*PM2) / real(NN,dk)
            PM2  = PM1  
            PM1  = P  
         ENDDO
!                                              ** Newton Method  
         TMP    = rONE / (rONE - X**2)  
         PPR    = EN*(PM2 - X*P)*TMP  
         P2PRI  = (rTWO*X*PPR - NNP1*P)*TMP  
         XI     = X - (P/PPR)*(rONE + (P / PPR )*P2PRI / (rTWO*PPR))  
  
!                                              ** Check for convergence  
         IF( ABS( XI - X ) > TOL ) THEN  
            IF( ITER.GT.MAXIT ) CALL ERRMSG( 'QGAUSN--max iteration count',.True.)  
            X  = XI  
            GO TO  10  
         END IF  
!                             ** Iteration finished--calculate weights,  
!                                abscissae for (-1,1)  
         GMU(K) = -X  
         GWT(K) = rTWO / (TMP*(EN*PM2)**2)  
         GMU(NP1 - K) = -GMU(K)  
         GWT(NP1 - K) = GWT(K)  
      ENDDO
!                                    ** Set middle abscissa and weight  
!                                       for rules of odd order  
      IF( MOD( M,iTWO ).NE.iZERO ) THEN  
         GMU(LIM + iONE) = rZERO
         PROD   = rONE  
         DO K = 3_ik, M, iTWO  
           PROD = PROD * real(K,dk) / real(K - iONE,dk)  
         ENDDO
         GWT(LIM + iONE) = rTWO / PROD**2  
      END IF  
  
!                                        ** Convert from (-1,1) to (0,1)  
      DO K = iONE, M  
         GMU(K) = ONEHALF*GMU(K) + ONEHALF
         GWT(K) = ONEHALF*GWT(K)  
      ENDDO
  
      END SUBROUTINE QGAUSN  
  
      SUBROUTINE SETDIS( dsdh, nid, tausla, tauslau, mu2, &
                         CMU, CWT, DELTAM, DTAUC, DTAUCP, EXPBEA, FBEAM, &
                         FLYR, GL, HL, HLPR, IBCND, LAMBER, LAYRU, &
                         LYRCUT, NCUT, NLYR, &
                         NTAU, NN, NSTR, PLANK, NUMU, ONLYFL, OPRIM, &
                         PMOM, SSALB, TAUC, TAUCPR, UTAU, UTAUPR, UMU,& 
                         UMU0, USRTAU, USRANG )  
  
!          Perform miscellaneous setting-up operations  
  
!       INPUT :  all are DISORT input variables (see DOC file)  
  
!       OUTPUT:  NTAU,UTAU   if USRTAU = FALSE  
!                NUMU,UMU    if USRANG = FALSE  
!                CMU,CWT     computational polar angles and  
!                               corresponding quadrature weights  
!                EXPBEA      transmission of direct beam  
!                FLYR        truncated fraction in delta-M method  
!                GL          phase function Legendre coefficients multi-  
!                              plied by (2L+1) and single-scatter albedo  
!                HLPR        Legendre moments of surface bidirectional  
!                              reflectivity, times 2K+1  
!                LAYRU       Computational layer in which UTAU falls  
!                LYRCUT      flag as to whether radiation will be zeroed  
!                              below layer NCUT  
!                NCUT        computational layer where absorption  
!                              optical depth first exceeds  ABSCUT  
!                NN          NSTR / 2  
!                OPRIM       delta-M-scaled single-scatter albedo  
!                TAUCPR      delta-M-scaled optical depth  
!                UTAUPR      delta-M-scaled version of  UTAU  
  
!   Called by- DISORT  
!   Calls- QGAUSN, ERRMSG  
! ----------------------------------------------------------------------  
  
      use tuv_params, only : largest  
  
!     .. Scalar Arguments ..  
  
      integer(ik), intent(in)  ::  IBCND, NLYR, NSTR  
      integer(ik), intent(out) ::  NCUT, NN, NTAU, NUMU  
      logical(lk), intent(in)  ::  DELTAM, LAMBER, ONLYFL  
      logical(lk), intent(in)  ::  PLANK, USRANG, USRTAU  
      logical(lk), intent(out) ::  LYRCUT  
      real(dk), intent(in)     ::  FBEAM, UMU0  
  
! geometry  
      integer(ik), intent(in) :: nid(0:)  
      real(dk), intent(in)    :: dsdh(0:,:)  
      real(dk), intent(out)   :: tausla(0:), tauslau(0:), mu2(0:)  
  
      real(dk) :: sum, sumu  
!     ..  
!     .. Array Arguments ..  
  
      integer(ik), intent(out) :: LAYRU(:)  
      real(dk), intent(in)   :: DTAUC(:), HL(0:), SSALB(:), TAUC(0:)  
      real(dk), intent(out)  :: UTAU(:), UTAUPR(:)  
      real(dk), intent(out)  :: CMU(:), CWT(:), DTAUCP(:), EXPBEA(0:)  
      real(dk), intent(out)  :: FLYR(:), GL(0:,:), HLPR(0:), UMU(:)  
      real(dk), intent(out)  :: OPRIM(:), PMOM(0:,:), TAUCPR(0:)   
  
!     ..  
!     .. Local Scalars ..  
  
      real(dk), PARAMETER :: ABSCUT = 10000._dk
  
      integer(ik)   :: IQ, IU, K, LC, LU, I  
      real(dk)      :: ABSTAU, F  
  
      IF( .NOT. USRTAU ) THEN  
!  ** Set output levels at computational layer boundaries  
         NTAU  = NLYR + iONE  
         DO LC = iZERO, NTAU - iONE  
            UTAU(LC + 1) = TAUC(LC)  
         ENDDO  
      END IF  
!                        ** Apply delta-M scaling and move description  
!                           of computational layers to local variables  
      EXPBEA(iZERO) = rONE  
      TAUCPR(iZERO) = rZERO  
      ABSTAU        = rZERO  
  
      tausla  = rZERO  
      tauslau = rZERO  
      mu2 = rONE/largest  
  
      DO LC = iONE, NLYR  
         PMOM( iZERO, LC ) = rONE  
         IF( ABSTAU < ABSCUT ) NCUT  = LC  
  
         ABSTAU = ABSTAU + (rONE - SSALB(LC))*DTAUC(LC)  
  
         IF( .NOT. DELTAM ) THEN  
            OPRIM(LC)  = SSALB(LC)  
            DTAUCP(LC) = DTAUC(LC)  
            TAUCPR(LC) = TAUC(LC)  
  
            DO K = iZERO, NSTR - iONE
               GL(K,LC) = real(iTWO*K + iONE,dk)*OPRIM(LC)*PMOM(K,LC)  
            ENDDO  
  
            F  = rZERO  
         ELSE  
!                                    ** Do delta-M transformation  
            F  = PMOM(NSTR,LC)  
            OPRIM(LC)  = SSALB(LC) * (rONE - F) / (rONE - F*SSALB(LC))  
            DTAUCP(LC) = (rONE - F*SSALB(LC))*DTAUC(LC)  
            TAUCPR(LC) = TAUCPR(LC - iONE) + DTAUCP(LC)  
  
            DO K = iZERO, NSTR - iONE
               GL(K,LC) = real((iTWO*K + iONE,dk) * OPRIM(LC) * (PMOM(K,LC) - F) / (rONE - F)  
            ENDDO  
         END IF  
         FLYR(LC) = F  
         EXPBEA(LC) = rZERO  
      ENDDO  
!   
! calculate slant optical depth  
!                
         IF(umu0 <  rZERO) THEN  
           IF(nid(iZERO) < iZERO) THEN  
             tausla(iZERO)  = largest  
             tauslau(iZERO) = largest  
           ELSE  
             sum  = rZERO  
             sumu = rZERO  
             DO lc = iONE, nid(iZERO)  
               sum  = sum + rTWO*dtaucp(lc)*dsdh(iZERO,lc)  
               sumu = sumu + rTWO*dtauc(lc)*dsdh(iZERO,lc)  
             END DO  
             tausla(iZERO)  = sum   
             tauslau(iZERO) = sumu   
           END IF  
         END IF  
  
         expbea(iZERO) = EXP( -tausla(iZERO) )  
  
!  
         DO lc = iONE, nlyr  
           IF(nid(lc) < iZERO) THEN  
             tausla(lc)  = largest  
             tauslau(lc) = largest  
           ELSE  
             sum  = rZERO  
             sumu = rZERO  
             DO lu = iONE, MIN(nid(lc),lc)  
               sum  = sum + dtaucp(lu)*dsdh(lc,lu)  
               sumu = sumu + dtauc(lu)*dsdh(lc,lu)  
             ENDDO  
             DO lu = MIN(nid(lc),lc)+1,nid(lc)  
               sum  = sum + rTWO*dtaucp(lu)*dsdh(lc,lu)  
               sumu = sumu + rTWO*dtauc(lu)*dsdh(lc,lu)  
             ENDDO  
             tausla(lc) = sum   
             tauslau(lc) = sumu   
             IF(tausla(lc) == tausla(lc-iONE)) THEN  
               mu2(lc) = largest  
             ELSE  
               mu2(lc) = (taucpr(lc) - taucpr(lc-iONE))/(tausla(lc) - tausla(lc-iONE))  
               mu2(lc) = SIGN( MAX(ABS(mu2(lc)),rONE/largest),mu2(lc) )  
             END IF  
           END IF  
           expbea(lc) = EXP( -tausla(lc) )  
         ENDDO  
  
!                      ** If no thermal emission, cut off medium below  
!                         absorption optical depth = ABSCUT ( note that  
!                         delta-M transformation leaves absorption  
!                         optical depth invariant ).  Not worth the  
!                         trouble for one-layer problems, though.  
  
      LYRCUT = ABSTAU >= ABSCUT .AND. .NOT. PLANK .AND. IBCND /= iONE .AND. NLYR > iONE  
  
      IF( .NOT.LYRCUT ) NCUT = NLYR  
  
!                             ** Set arrays defining location of user  
!                             ** output levels within delta-M-scaled  
!                             ** computational mesh  
      DO LU = iONE, NTAU  
         DO LC = iONE, NLYR  
            IF( UTAU(LU) >= TAUC(LC - iONE ) .AND. UTAU(LU) <= TAUC(LC) ) EXIT  
         ENDDO  
         LC = MIN( NLYR,LC )  
         UTAUPR(LU) = UTAU(LU)  
         IF( DELTAM ) THEN  
           UTAUPR( LU ) = TAUCPR(LC - iONE) + (rONE - SSALB(LC)*FLYR(LC)) &
                          * (UTAU(LU) - TAUC(LC-iONE))  
         ENDIF  
         LAYRU(LU) = LC  
      ENDDO  
!                      ** Calculate computational polar angle cosines  
!                         and associated quadrature weights for Gaussian  
!                         quadrature on the interval (0,1) (upward)  
      NN   = NSTR / iTWO
  
      CALL QGAUSN( NN, CMU, CWT )  
!                                  ** Downward (neg) angles and weights  
      DO IQ = iONE, NN  
         CMU(IQ + NN) = -CMU(IQ)  
         CWT(IQ + NN) = CWT(IQ)  
      ENDDO  
  
  
      DO IQ = iONE, NN  
!                      ** Dither mu2 if it is close to one of the   
!                         quadrature angles.  
        DO  lc = 1, nlyr  
          IF (  ABS(mu2(lc)) < 1.E5_dk ) THEN  
            IF( ABS(rONE - ABS(mu2(lc))/CMU(IQ)) < 0.05_dk ) mu2(lc) = mu2(lc)*0.999_dk
          ENDIF  
        END DO  
      END DO  
  
      IF( .NOT. USRANG .OR. ( ONLYFL .AND. MXUMU >= NSTR ) ) THEN  
!                                   ** Set output polar angles to  
!                                      computational polar angles  
        NUMU = NSTR  
        DO IU = iONE, NN  
          UMU(IU) = -CMU(NN + iONE - IU)  
        ENDDO  
  
        DO IU = NN + iONE, NSTR  
          UMU(IU) = CMU(IU - NN)  
        ENDDO  
      END IF  
  
  
      IF( USRANG .AND. IBCND == iONE ) THEN  
!                               ** Shift positive user angle cosines to  
!                                  upper locations and put negatives  
!                                  in lower locations  
         DO IU = iONE, NUMU  
            UMU(IU+NUMU) = UMU(IU)  
         ENDDO  
  
         DO IU = iONE, NUMU  
            UMU(IU) = -UMU(iTWO*NUMU + iONE - IU)  
         ENDDO  
         NUMU   = iTWO*NUMU  
      END IF  
  
  
      IF( .NOT. LYRCUT .AND. .NOT.LAMBER ) THEN  
         DO K = iZERO, NSTR  
            HLPR(K) = real(iTWO*K + iONE,dk)*HL(K)  
         ENDDO  
      END IF  
  
      END SUBROUTINE SETDIS  
  
      SUBROUTINE SETMTX( BDR, CBAND, CMU, CWT, DELM0, DTAUCP, GC, KK, &
                         LAMBER, LYRCUT, NCOL, NCUT, &
                         NN, NSTR, TAUCPR, WK )  
  
!        Calculate coefficient matrix for the set of equations  
!        obtained from the boundary conditions and the continuity-  
!        of-intensity-at-layer-interface equations;  store in the  
!        special banded-matrix format required by LINPACK routines  
  
!     I N P U T      V A R I A B L E S:  
  
!       BDR      :  Surface bidirectional reflectivity  
!       CMU      :  Abscissae for Gauss quadrature over angle cosine  
!       CWT      :  Weights for Gauss quadrature over angle cosine  
!       DELM0    :  Kronecker delta, delta-sub-m0  
!       GC       :  Eigenvectors at polar quadrature angles, SC(1)  
!       KK       :  Eigenvalues of coeff. matrix in Eq. SS(7)  
!       LYRCUT   :  Logical flag for truncation of comput. layer  
!       NN       :  Number of streams in a hemisphere (NSTR/2)  
!       NCUT     :  Total number of computational layers considered  
!       TAUCPR   :  Cumulative optical depth (delta-M-scaled)  
!       (remainder are DISORT input variables)  
  
!   O U T P U T     V A R I A B L E S:  
  
!       CBAND    :  Left-hand side matrix of linear system Eq. SC(5),  
!                      scaled by Eq. SC(12); in banded form required  
!                      by LINPACK solution routines  
!       NCOL     :  Counts of columns in CBAND  
  
!   I N T E R N A L    V A R I A B L E S:  
  
!       IROW     :  Points to row in CBAND  
!       JCOL     :  Points to position in layer block  
!       LDA      :  Row dimension of CBAND  
!       NCD      :  Number of diagonals below or above main diagonal  
!       NSHIFT   :  For positioning number of rows in band storage  
!       WK       :  Temporary storage for EXP evaluations  
  
!   Called by- DISORT, ALBTRN  
! +--------------------------------------------------------------------+  
  
  
!     .. Scalar Arguments ..  
  
      logical(lk), intent(in)  :: LAMBER, LYRCUT  
      integer(ik), intent(in)  :: NCUT, NN, NSTR  
      integer(ik), intent(out) :: NCOL  
      real(dk), intent(in)     :: DELM0  
!     ..  
!     .. Array Arguments ..  
  
      real(dk), intent(in)  ::  BDR(:,0:), CMU(:),  
                CWT(:), DTAUCP(:), GC(:,:,:),  
                KK(:,:), TAUCPR(0:)  
      real(dk), intent(inout) :: WK(:)  
      real(dk), intent(out)   :: CBAND(:,:)  
!     ..  
!     .. Local Scalars ..  
  
      integer(ik)   :: IQ, IROW, JCOL, JQ, K, LC, LDA, NCD, NNCOL, NSHIFT  
      real(dk)      :: EXPA, SUM  
!     ..  
  
      CBAND = rZERO  
  
      NCD    = iTHREE*NN - iONE  
      LDA    = iTHREE*NCD + iONE  
      NSHIFT = LDA - iTWO*NSTR + iONE  
      NCOL   = iZERO
!                         ** Use continuity conditions of Eq. STWJ(17)  
!                            to form coefficient matrix in STWJ(20);  
!                            employ scaling transformation STWJ(22)  
      DO LC = iONE, NCUT  
         DO IQ = iONE, NN  
            WK(IQ) = EXP( KK(IQ,LC)*DTAUCP(LC) )  
         ENDDO
  
         JCOL  = iZERO
         DO IQ = iONE, NN  
            NCOL  = NCOL + iONE  
            IROW  = NSHIFT - JCOL  
            DO JQ = iONE, NSTR  
               CBAND(IROW + NSTR,NCOL) =   GC(JQ,IQ,LC)  
               CBAND(IROW,NCOL)        = - GC(JQ,IQ,LC)*WK(IQ)  
               IROW  = IROW + iONE  
            ENDDO
            JCOL  = JCOL + iONE  
         ENDDO
  
         DO IQ = NN + iONE, NSTR  
            NCOL  = NCOL + iONE  
            IROW  = NSHIFT - JCOL  
            DO JQ = iONE, NSTR  
               CBAND(IROW + NSTR,NCOL) =   GC(JQ,IQ,LC) * WK(NSTR + iONE - IQ)  
               CBAND(IROW,NCOL)        = - GC(JQ,IQ,LC)  
               IROW  = IROW + iONE  
            ENDDO
            JCOL  = JCOL + iONE  
         ENDDO
      ENDDO
!                  ** Use top boundary condition of STWJ(20a) for  
!                     first layer  
  
      JCOL  = iZERO
  
      DO IQ = iONE, NN  
         EXPA  = EXP( KK(IQ,iONE)*TAUCPR(iONE) )  
         IROW  = NSHIFT - JCOL + NN  
         DO JQ = NN, iONE, -iONE  
            CBAND(IROW,JCOL + iONE) = GC(JQ,IQ,iONE)*EXPA  
            IROW  = IROW + iONE  
         ENDDO
         JCOL  = JCOL + iONE  
      ENDDO
  
  
      DO IQ = NN + iONE, NSTR  
         IROW  = NSHIFT - JCOL + NN  
         DO JQ = NN, iONE, -iONE  
            CBAND(IROW,JCOL + iONE) = GC(JQ,IQ,iONE)  
            IROW  = IROW + iONE  
         ENDDO
         JCOL  = JCOL + iONE  
      ENDDO
!                           ** Use bottom boundary condition of  
!                              STWJ(20c) for last layer  
  
      NNCOL = NCOL - NSTR  
      JCOL  = iZERO
  
      DO IQ = iONE, NN  
         NNCOL  = NNCOL + iONE  
         IROW   = NSHIFT - JCOL + NSTR  
  
         DO JQ = NN + iONE, NSTR  
            IF( LYRCUT .OR. (LAMBER .AND. DELM0 == rZERO) ) THEN  
!                          ** No azimuthal-dependent intensity if Lam-  
!                             bert surface; no intensity component if  
!                             truncated bottom layer  
  
               CBAND(IROW,NNCOL) = GC(JQ,IQ,NCUT)  
            ELSE  
               SUM  = rZERO
               DO K = iONE, NN  
                  SUM  = SUM + CWT(K)*CMU(K)*BDR(JQ - NN,K) * GC(NN + iONE - K,IQ,NCUT)  
               ENDDO
               CBAND(IROW,NNCOL) = GC(JQ,IQ,NCUT) -  (rONE + DELM0)*SUM  
            END IF  
            IROW  = IROW + iONE  
         ENDDO
         JCOL  = JCOL + iONE  
      ENDDO
  
  
      DO IQ = NN + iONE, NSTR  
         NNCOL  = NNCOL + iONE  
         IROW   = NSHIFT - JCOL + NSTR  
         EXPA   = WK(NSTR + iONE - IQ)  
         DO JQ = NN + iONE, NSTR  
            IF( LYRCUT .OR. (LAMBER .AND. DELM0 == rZERO) ) THEN  
               CBAND(IROW,NNCOL) = GC(JQ,IQ,NCUT)*EXPA  
            ELSE  
               SUM  = rZERO
               DO K = iONE, NN  
                  SUM  = SUM + CWT(K)*CMU(K)*BDR(JQ - NN,K) * GC(NN + iONE - K,IQ,NCUT)  
               ENDDO
               CBAND(IROW,NNCOL) = (GC(JQ,IQ,NCUT) -  (rONE + DELM0)*SUM)*EXPA  
            END IF  
            IROW  = IROW + iONE  
         ENDDO
         JCOL  = JCOL + iONE  
      ENDDO
  
      END SUBROUTINE SETMTX  
  
      SUBROUTINE SOLEIG( AMB, APB, ARRAY, CMU, CWT, GL, MAZIM, &
                         NN, NSTR, YLMC, CC, EVECC, EVAL, KK, GC, &
                         AAD, EVECCD, EVALD, WKD )  
  
!         Solves eigenvalue/vector problem necessary to construct  
!         homogeneous part of discrete ordinate solution; STWJ(8b)  
!         ** NOTE ** Eigenvalue problem is degenerate when single  
!                    scattering albedo = 1;  present way of doing it  
!                    seems numerically more stable than alternative  
!                    methods that we tried  
  
!   I N P U T     V A R I A B L E S:  
  
!       GL     :  Delta-M scaled Legendre coefficients of phase function  
!                    (including factors 2l+1 and single-scatter albedo)  
!       CMU    :  Computational polar angle cosines  
!       CWT    :  Weights for quadrature over polar angle cosine  
!       MAZIM  :  Order of azimuthal component  
!       NN     :  Half the total number of streams  
!       YLMC   :  Normalized associated Legendre polynomial  
!                    at the quadrature angles CMU  
!       (remainder are DISORT input variables)  
  
!   O U T P U T    V A R I A B L E S:  
  
!       CC     :  C-sub-ij in Eq. SS(5); needed in SS(15&18)  
!       EVAL   :  NN eigenvalues of Eq. SS(12) on return from ASYMTX  
!                    but then square roots taken  
!       EVECC  :  NN eigenvectors  (G+) - (G-)  on return  
!                    from ASYMTX ( column j corresponds to EVAL(j) )  
!                    but then  (G+) + (G-)  is calculated from SS(10),  
!                    G+  and  G-  are separated, and  G+  is stacked on  
!                    top of  G-  to form NSTR eigenvectors of SS(7)  
!       GC     :  Permanent storage for all NSTR eigenvectors, but  
!                    in an order corresponding to KK  
!       KK     :  Permanent storage for all NSTR eigenvalues of SS(7),  
!                    but re-ordered with negative values first ( square  
!                    roots of EVAL taken and negatives added )  
  
!   I N T E R N A L   V A R I A B L E S:  
  
!       AMB,APB :  Matrices (alpha-beta), (alpha+beta) in reduced  
!                    eigenvalue problem  
!       ARRAY   :  Complete coefficient matrix of reduced eigenvalue  
!                    problem: (alfa+beta)*(alfa-beta)  
!       GPPLGM  :  (G+) + (G-) (cf. Eqs. SS(10-11))  
!       GPMIGM  :  (G+) - (G-) (cf. Eqs. SS(10-11))  
!       WKD     :  Scratch array required by ASYMTX  
  
!   Called by- DISORT, ALBTRN  
!   Calls- ASYMTX, ERRMSG  
! +-------------------------------------------------------------------+  
  
  
!     .. Scalar Arguments ..  
  
      integer(ik), intent(in) ::  MAZIM, NN, NSTR  
!     ..  
!     .. Array Arguments ..  
  
      real(dk), intent(out) ::  EVAL(:), KK(:)  
      real(dk), intent(out) ::  CC(:,:), EVECC(:,:), GC(:,:)  
      real(dk), intent(out) ::  AMB(:,:), APB(:,:), ARRAY(:,:)  
      real(dk), intent(in)  ::  CMU(:), CWT(:), GL(0:), YLMC(0:,:)  
      real(dk)(8), intent(inout)  :: AAD(:,:), EVALD(:), EVECCD(:,:), WKD(:)  
  
!     ..  
!     .. Local Scalars ..  
  
      integer(ik) :: IER, IQ, JQ, KQ, L  
      real(dk)    :: ALPHA, BETA, GPMIGM, GPPLGM, SUM  
  
!                             ** Calculate quantities in Eqs. SS(5-6)  
      DO IQ = iONE, NN  
         DO JQ = iONE, NSTR  
            SUM  = rZERO
            DO L = MAZIM, NSTR - iONE  
               SUM  = SUM + GL(L)*YLMC(L,IQ)*YLMC(L,JQ)  
            ENDDO 
            CC(IQ,JQ) = ONEHALF*SUM*CWT(JQ)  
         ENDDO 
  
         DO JQ = iONE, NN  
!                             ** Fill remainder of array using symmetry  
!                                relations  C(-mui,muj) = C(mui,-muj)  
!                                and        C(-mui,-muj) = C(mui,muj)  
            CC(IQ + NN,JQ) = CC(IQ,JQ + NN)  
            CC(IQ + NN,JQ + NN) = CC(IQ,JQ)  
!                                       ** Get factors of coeff. matrix  
!                                          of reduced eigenvalue problem  
            ALPHA  = CC(IQ,JQ) / CMU(IQ)  
            BETA   = CC(IQ,JQ + NN) / CMU(IQ)  
            AMB(IQ,JQ) = ALPHA - BETA  
            APB(IQ,JQ) = ALPHA + BETA  
         ENDDO 
         AMB(IQ,IQ) = AMB(IQ,IQ) - rONE / CMU(IQ)  
         APB(IQ,IQ) = APB(IQ,IQ) - rONE / CMU(IQ)  
      ENDDO 
!                      ** Finish calculation of coefficient matrix of  
!                         reduced eigenvalue problem:  get matrix  
!                         product (alfa+beta)*(alfa-beta); SS(12)  
      DO IQ = iONE, NN  
         DO JQ = iONE, NN  
            SUM  = rZERO
            DO KQ = iONE, NN  
               SUM  = SUM + APB(IQ,KQ)*AMB(KQ,JQ)  
            ENDDO 
            ARRAY(IQ,JQ) = SUM  
         ENDDO 
      ENDDO 
!                      ** Find (real) eigenvalues and eigenvectors  
  
      CALL ASYMTX( ARRAY, EVECC, EVAL, IER, WKD, AAD, EVECCD, EVALD )  
  
      IF( IER.GT.iZERO ) THEN  
         WRITE( *, FMT = '(//,A,I4,A)' ) ' ASYMTX--eigenvalue no. ', &
            IER, '  didnt converge.  Lower-numbered eigenvalues wrong.'  
         CALL ERRMSG( 'ASYMTX--convergence problems',.True.)  
      END IF  
  
      DO IQ = iONE, NN  
         EVAL(IQ)    = SQRT( ABS( EVAL(IQ) ) )  
         KK(IQ + NN) = EVAL(IQ)  
!                                      ** Add negative eigenvalue  
         KK(NN + iONE - IQ) = -EVAL(IQ)  
      ENDDO
  
!                          ** Find eigenvectors (G+) + (G-) from SS(10)  
!                             and store temporarily in APB array  
      DO JQ = iONE, NN  
         DO IQ = iONE, NN  
            SUM  = rZERO
            DO KQ = iONE, NN  
               SUM  = SUM + AMB(IQ,KQ)*EVECC(KQ,JQ)  
            ENDDO
            APB(IQ,JQ) = SUM / EVAL(JQ)  
         ENDDO 
      ENDDO 
  
  
      DO JQ = iONE, NN  
         DO IQ = iONE, NN  
            GPPLGM = APB(IQ,JQ)  
            GPMIGM = EVECC(IQ,JQ)  
!                                ** Recover eigenvectors G+,G- from  
!                                   their sum and difference; stack them  
!                                   to get eigenvectors of full system  
!                                   SS(7) (JQ = eigenvector number)  
            EVECC(IQ,      JQ ) = ONEHALF*(GPPLGM + GPMIGM)  
            EVECC(IQ + NN, JQ ) = ONEHALF*(GPPLGM - GPMIGM)  
!                                ** Eigenvectors corresponding to  
!                                   negative eigenvalues (corresp. to  
!                                   reversing sign of 'k' in SS(10) )  
            GPPLGM = - GPPLGM  
            EVECC(IQ,   JQ+NN) = ONEHALF * (GPPLGM + GPMIGM)  
            EVECC(IQ+NN,JQ+NN) = ONEHALF * (GPPLGM - GPMIGM)  
            GC(IQ+NN,   JQ+NN) = EVECC(IQ,JQ)  
            GC(NN+iONE-IQ, JQ+NN)      = EVECC(IQ+NN, JQ)  
            GC(IQ+NN,   NN+iONE-JQ)    = EVECC(IQ,    JQ+NN)  
            GC(NN+iONE-IQ, NN+iONE-JQ) = EVECC(IQ+NN, JQ+NN)  
         ENDDO 
      ENDDO 
  
      END SUBROUTINE SOLEIG  
  
      SUBROUTINE SOLVE0( B, BDR, BEM, BPLANK, CBAND, CMU, CWT, EXPBEA, &
                         FBEAM, FISOT, IPVT, LAMBER, LL, LYRCUT, MAZIM, &
                         NCOL, NCUT, NN, NSTR, &
                         PI, TPLANK, TAUCPR, UMU0, Z, ZZ, ZPLK0, ZPLK1 )  
  
      use linalgebra, only : linalgebra_t  
  
!        Construct right-hand side vector B for general boundary  
!        conditions STWJ(17) and solve system of equations obtained  
!        from the boundary conditions and the continuity-of-  
!        intensity-at-layer-interface equations.  
!        Thermal emission contributes only in azimuthal independence.  
  
!     I N P U T      V A R I A B L E S:  
  
!       BDR      :  Surface bidirectional reflectivity  
!       BEM      :  Surface bidirectional emissivity  
!       BPLANK   :  Bottom boundary thermal emission  
!       CBAND    :  Left-hand side matrix of linear system Eq. SC(5),  
!                   scaled by Eq. SC(12); in banded form required  
!                   by LINPACK solution routines  
!       CMU      :  Abscissae for Gauss quadrature over angle cosine  
!       CWT      :  Weights for Gauss quadrature over angle cosine  
!       EXPBEA   :  Transmission of incident beam, EXP(-TAUCPR/UMU0)  
!       LYRCUT   :  Logical flag for truncation of comput. layer  
!       MAZIM    :  Order of azimuthal component  
!       ncol     :  Counts of columns in CBAND  
!       NN       :  Order of double-Gauss quadrature (NSTR/2)  
!       NCUT     :  Total number of computational layers considered  
!       TPLANK   :  Top boundary thermal emission  
!       TAUCPR   :  Cumulative optical depth (delta-M-scaled)  
!       ZZ       :  Beam source vectors in Eq. SS(19)  
!       ZPLK0    :  Thermal source vectors Z0, by solving Eq. SS(16)  
!       ZPLK1    :  Thermal source vectors Z1, by solving Eq. SS(16)  
!       (remainder are DISORT input variables)  
  
!   O U T P U T     V A R I A B L E S:  
  
!       B        :  Right-hand side vector of Eq. SC(5) going into  
!                   SGBSL; returns as solution vector of Eq. SC(12),  
!                   constants of integration without exponential term  
!  
!      LL        :  Permanent storage for B, but re-ordered  
  
!   I N T E R N A L    V A R I A B L E S:  
  
!       IPVT     :  Integer vector of pivot indices  
!       IT       :  Pointer for position in  B  
!       NCD      :  Number of diagonals below or above main diagonal  
!       RCOND    :  Indicator of singularity for CBAND  
!       Z        :  Scratch array required by SGBCO  
  
!   Called by- DISORT  
!   Calls- SGBCO, ERRMSG, SGBSL  
! +-------------------------------------------------------------------+  
  
  
!     .. Scalar Arguments ..  
  
      logical(lk), intent(in) :: LAMBER, LYRCUT  
      integer(ik), intent(in) :: MAZIM, NCOL, NCUT, NN, NSTR  
      real(dk),    intent(in) :: BPLANK, FBEAM, FISOT, PI, TPLANK, UMU0  
!     ..  
!     .. Array Arguments ..  
  
      integer(ik), intent(inout) :: IPVT(:)  
      real(dk), intent(out)      :: B(:), LL(:,:)  
      real(dk), intent(inout)    :: Z(:)  
      real(dk), intent(inout)    :: CBAND(:,:)  
      real(dk), intent(in)       :: BDR(:,0:), BEM(:), &
                                    CMU(:), CWT(:), &
                                    EXPBEA(0:), TAUCPR(0:), &
                                    ZPLK0(:,:), ZPLK1(:,:), ZZ(:,:)  
!     ..  
!     .. Local Scalars ..  
  
      integer(ik) ::   IPNT, IQ, IT, JQ, LC, NCD  
      real(dk)    :: RCOND, SUM  
      TYPE(linalgebra_t) :: linpack  
  
      B = rZERO  
!                              ** Construct B,  STWJ(20a,c) for  
!                                 parallel beam + bottom reflection +  
!                                 thermal emission at top and/or bottom  
  
      IF( MAZIM.GT.iZERO .AND. FBEAM.GT.rZERO ) THEN  
!                                         ** Azimuth-dependent case  
!                                            (never called if FBEAM = 0)  
         IF( LYRCUT .OR. LAMBER ) THEN  
!               ** No azimuthal-dependent intensity for Lambert surface;  
!                  no intensity component for truncated bottom layer  
            DO IQ = iONE, NN  
!                                                  ** Top boundary  
               B(IQ) = - ZZ(NN + iONE - IQ,iONE)*EXPBEA(iZERO)  
!                                                  ** Bottom boundary  
               B(NCOL - NN + IQ) = -ZZ(IQ + NN,NCUT)*EXPBEA(NCUT)  
            ENDDO
         ELSE  
            DO IQ = iONE, NN  
               B(IQ) = - ZZ(NN + iONE - IQ,iONE)*EXPBEA(iZERO)  
               SUM  = rZERO
               DO JQ = iONE, NN  
                  SUM  = SUM + CWT(JQ)*CMU(JQ)*BDR(IQ,JQ) * ZZ(NN + iONE - JQ,NCUT)*EXPBEA(NCUT)  
               ENDDO
               B(NCOL - NN + IQ) = SUM  
               IF( FBEAM.GT.rZERO ) THEN
                  B(NCOL - NN + IQ) = SUM + (BDR(IQ,iZERO)*UMU0*FBEAM / PI - ZZ(IQ + NN,NCUT)) &
                                            * EXPBEA(NCUT)  
               ENDIF
            ENDDO
         END IF  
!                             ** Continuity condition for layer  
!                                interfaces of Eq. STWJ(20b)  
         IT = NN  
         DO LC = iONE, NCUT - iONE  
            DO IQ = iONE, NSTR  
               IT   = IT + iONE  
               B(IT) = (ZZ(IQ,LC+iONE) - ZZ(IQ,LC))*EXPBEA(LC)  
            ENDDO
         ENDDO
      ELSE  
!                                   ** Azimuth-independent case  
         IF( FBEAM.EQ.rZERO ) THEN  
            DO IQ = iONE, NN  
!                                      ** Top boundary  
               B(IQ) = -ZPLK0(NN + iONE - IQ, iONE) + FISOT + TPLANK  
            ENDDO
            IF( LYRCUT ) THEN  
!                               ** No intensity component for truncated  
!                                  bottom layer  
               DO IQ = iONE, NN  
!                                      ** Bottom boundary  
                  B(NCOL - NN + IQ) = - ZPLK0(IQ + NN,NCUT) - ZPLK1(IQ + NN,NCUT) * TAUCPR(NCUT)  
               ENDDO
            ELSE  
               DO IQ = iONE, NN  
                  SUM  = 0.  
                  DO JQ = iONE, NN  
                     SUM  = SUM + CWT(JQ)*CMU(JQ)*BDR(IQ,JQ) &
                                * (ZPLK0(NN + iONE - JQ,NCUT) + ZPLK1(NN + iONE - JQ,NCUT)*TAUCPR(NCUT))  
                  ENDDO
                  B(NCOL - NN + IQ) = rTWO*SUM + BEM(IQ)*BPLANK - ZPLK0(IQ + NN,NCUT) &
                                               - ZPLK1(IQ + NN,NCUT) * TAUCPR(NCUT)  
               ENDDO
            END IF  
!                             ** Continuity condition for layer  
!                                interfaces, STWJ(20b)  
            IT   = NN  
            DO LC = iONE, NCUT - iONE  
               DO IQ = iONE, NSTR  
                  IT   = IT + iONE  
                  B(IT) = ZPLK0(IQ,LC + iONE) - ZPLK0(IQ,LC) &
                        + (ZPLK1(IQ,LC + iONE) - ZPLK1(IQ,LC)) * TAUCPR(LC)  
               ENDDO 
            ENDDO 
         ELSE  
            DO IQ = iONE, NN  
               B(IQ) = - ZZ(NN + iONE - IQ,iONE)*EXPBEA(iZERO) &
                       - ZPLK0(NN + iONE - IQ,iONE) + FISOT + TPLANK  
            ENDDO
            IF( LYRCUT ) THEN  
               DO IQ = iONE, NN  
                  B(NCOL-NN+IQ) = - ZZ(IQ+NN, NCUT) * EXPBEA(NCUT) &
                                  - ZPLK0(IQ+NN,NCUT) &
                                  - ZPLK1(IQ+NN,NCUT) * TAUCPR(NCUT)  
               ENDDO
            ELSE  
               DO IQ = iONE, NN  
                  SUM  = rZERO
                  DO JQ = iONE, NN  
                     SUM = SUM + CWT(JQ) * CMU(JQ) * BDR(IQ,JQ) &
                                *(ZZ(NN+iONE-JQ,NCUT) * EXPBEA(NCUT) &
                                  + ZPLK0(NN+iONE-JQ,NCUT) &
                                  + ZPLK1(NN+iONE-JQ,NCUT) * TAUCPR(NCUT))  
                  ENDDO
                  B(NCOL-NN+IQ) = rTWO*SUM + (BDR(IQ,iZERO) * UMU0*FBEAM/PI - ZZ(IQ+NN,NCUT)) * EXPBEA(NCUT) &
                                + BEM(IQ) * BPLANK  - ZPLK0(IQ+NN,NCUT) - ZPLK1(IQ+NN,NCUT) * TAUCPR(NCUT)  
               ENDDO
            END IF  
  
            IT   = NN  
            DO LC = iONE, NCUT - iONE  
               DO IQ = iONE, NSTR  
                  IT   = IT + iONE  
                  B(IT) = (ZZ(IQ,LC+iONE) - ZZ(IQ,LC)) * EXPBEA(LC) &
                        + ZPLK0(IQ,LC+iONE) - ZPLK0(IQ,LC) &
                        + (ZPLK1(IQ,LC+iONE) - ZPLK1(IQ,LC)) * TAUCPR(LC)  
               ENDDO 
            ENDDO 
         END IF  
      END IF  
!                     ** Find L-U (lower/upper triangular) decomposition  
!                        of band matrix CBAND and test if it is nearly  
!                        singular (note: CBAND is destroyed)  
!                        (CBAND is in LINPACK packed format)  
      RCOND  = rZERO  
      NCD    = iTHREE*NN - iONE
  
      CALL linpack%SGBCO( CBAND, NCOL, NCD, NCD, IPVT, RCOND, Z )  
  
      IF( rONE + RCOND == rONE )  
          CALL ERRMSG('SOLVE0--SGBCO says matrix near singular',.FALSE.)  
  
!                   ** Solve linear system with coeff matrix CBAND  
!                      and R.H. side(s) B after CBAND has been L-U  
!                      decomposed.  Solution is returned in B.  
  
      CALL linpack%SGBSL( CBAND, NCOL, NCD, NCD, IPVT, B, 0 )  
  
!                   ** Zero CBAND (it may contain 'foreign'  
!                      elements upon returning from LINPACK);  
!                      necessary to prevent errors  
  
      CBAND = rZERO  
  
      DO LC = iONE, NCUT  
         IPNT  = LC*NSTR - NN  
         DO IQ = iONE, NN  
            LL(NN + iONE - IQ,LC ) = B(IPNT + iONE - IQ)  
            LL(IQ + NN,       LC ) = B(IQ + IPNT)  
         ENDDO 
      ENDDO 
  
      END SUBROUTINE SOLVE0  
  
      SUBROUTINE SURFAC( ALBEDO, DELM0, FBEAM, HLPR, LAMBER, MAZIM, &
                         NN, NUMU, NSTR, ONLYFL, UMU, &
                         USRANG, YLM0, YLMC, YLMU, BDR, EMU, BEM, RMU )  
  
!       Specifies user's surface bidirectional properties, STWJ(21)  
  
!   I N P U T     V A R I A B L E S:  
  
!       DELM0  :  Kronecker delta, delta-sub-m0  
!       HLPR   :  Legendre moments of surface bidirectional reflectivity  
!                    (with 2K+1 factor included)  
!       MAZIM  :  Order of azimuthal component  
!       NN     :  Order of double-Gauss quadrature (NSTR/2)  
!       YLM0   :  Normalized associated Legendre polynomial  
!                 at the beam angle  
!       YLMC   :  Normalized associated Legendre polynomials  
!                 at the quadrature angles  
!       YLMU   :  Normalized associated Legendre polynomials  
!                 at the user angles  
!       (remainder are DISORT input variables)  
  
!    O U T P U T     V A R I A B L E S:  
  
!       BDR :  Surface bidirectional reflectivity (computational angles)  
!       RMU :  Surface bidirectional reflectivity (user angles)  
!       BEM :  Surface directional emissivity (computational angles)  
!       EMU :  Surface directional emissivity (user angles)  
  
!    I N T E R N A L     V A R I A B L E S:  
  
!       DREF      Directional reflectivity  
!       NMUG   :  Number of angle cosine quadrature points on (0,1) for  
!                   integrating bidirectional reflectivity to get  
!                   directional emissivity (it is necessary to use a  
!                   quadrature set distinct from the computational  
!                   angles, because the computational angles may not be  
!                   dense enough--NSTR may be too small--to give an  
!                   accurate approximation for the integration).  
!       GMU    :  The NMUG angle cosine quadrature points on (0,1)  
!       GWT    :  The NMUG angle cosine quadrature weights on (0,1)  
!       YLMG   :  Normalized associated Legendre polynomials  
!                   at the NMUG quadrature angles  
  
!   Called by- DISORT  
!   Calls- QGAUSN, LEPOLY, ERRMSG  
! +-------------------------------------------------------------------+  
  
!     .. Parameters ..  
  
      integer(ik), PARAMETER ::   NMUG = 10_ik, MAXSTR = 100_ik 
!     ..  
!     .. Scalar Arguments ..  
  
      logical(lk), intent(in) :: LAMBER, ONLYFL, USRANG  
      integer(ik), intent(in) :: MAZIM, NN, NSTR, NUMU  
      real(dk), intent(in)    :: ALBEDO, DELM0, FBEAM  
!     ..  
!     .. Array Arguments ..  
  
      real(dk), intent(in)  ::  HLPR(0:), UMU(:), &
                                YLM0(0:), YLMC(0:,:), YLMU(0:,:)  
      real(dk), intent(out) ::  BDR(:,0:), BEM(:), EMU(:), RMU(:,0:)  
!     ..  
!     .. Local Scalars ..  
  
      logical(lk) :: PASS1  
      integer(ik) :: IQ, IU, JG, JQ, K  
      real(dk)    :: DREF, SGN, SUM  
!     ..  
!     .. Local Arrays ..  
  
      real(dk)    :: GMU(NMUG), GWT(NMUG), YLMG(0:MAXSTR,NMUG)  
  
      SAVE      PASS1, GMU, GWT, YLMG  
      DATA      PASS1 / .TRUE._lk /  
  
  
      IF( PASS1 ) THEN  
         PASS1  = .FALSE._lk  
         CALL QGAUSN( NMUG, GMU, GWT )  
         CALL LEPOLY( NMUG, 0, MAXSTR, GMU, YLMG )  
!                       ** Convert Legendre polys. to negative GMU  
         SGN  = -rONE
         DO K = iZERO, MAXSTR  
            SGN  = -SGN  
            DO JG = iONE, NMUG  
               YLMG(K,JG) = SGN*YLMG(K,JG)  
            ENDDO  
         ENDDO  
      END IF  
  
      BDR = rZERO  
      BEM = rZERO  
  
      IF( LAMBER .AND. MAZIM == iZERO ) THEN  
         DO IQ = iONE, NN  
            BEM(IQ) = rONE - ALBEDO  
            DO JQ = iZERO, NN  
               BDR(IQ,JQ) = ALBEDO  
            ENDDO  
         ENDDO  
      ELSE IF( .NOT. LAMBER ) THEN  
!                                  ** Compute surface bidirectional  
!                                     properties at computational angles  
         DO IQ = iONE, NN  
            DO JQ = iONE, NN  
               SUM  = rZERO  
               DO K = MAZIM, NSTR - iONE  
                  SUM  = SUM + HLPR(K)*YLMC(K,IQ)*YLMC(K,JQ+NN)  
               ENDDO  
               BDR(IQ,JQ) = (rTWO - DELM0)*SUM  
            ENDDO  
  
            IF( FBEAM > rZERO ) THEN  
               SUM  = rZERO  
               DO K = MAZIM, NSTR - 1  
                  SUM  = SUM + HLPR(K)*YLMC(K,IQ)*YLM0(K)  
               ENDDO  
               BDR(IQ,iZERO) = (rTWO - DELM0)*SUM  
            END IF  
         ENDDO  
  
         IF( MAZIM == 0 ) THEN  
            IF( NSTR > MAXSTR )  
                CALL ERRMSG('SURFAC--parameter MAXSTR too small',.True.)  
!                              ** Integrate bidirectional reflectivity  
!                                 at reflection polar angles CMU and  
!                                 incident angles GMU to get  
!                                 directional emissivity at  
!                                 computational angles CMU.  
            DO IQ = iONE, NN  
               DREF  = rZERO  
               DO JG = iONE, NMUG  
                  SUM  = rZERO  
                  DO K = iZERO, NSTR - iONE  
                     SUM  = SUM + HLPR(K)*YLMC(K,IQ)*YLMG(K,JG)  
                  ENDDO  
                  DREF  = DREF + rTWO*GWT(JG)*GMU(JG)*SUM  
               ENDDO  
               BEM(IQ) = rONE - DREF  
            ENDDO  
         END IF  
      END IF  
!                                       ** Compute surface bidirectional  
!                                          properties at user angles  
  
      IF( .NOT. ONLYFL .AND. USRANG ) THEN  
         EMU = rZERO  
         RMU = rZERO  
         DO IU = iONE, NUMU  
            IF( UMU(IU) > rZERO ) THEN  
               IF( LAMBER .AND. MAZIM == iZERO ) THEN  
                  DO IQ = iZERO, NN  
                     RMU( IU, IQ ) = ALBEDO  
                  ENDDO  
                  EMU( IU ) = rONE - ALBEDO  
               ELSE IF( .NOT.LAMBER ) THEN  
                  DO IQ = iONE, NN  
                     SUM  = rZERO  
                     DO K = MAZIM, NSTR - iONE  
                        SUM = SUM + HLPR(K)*YLMU(K,IU)*YLMC(K,IQ + NN)  
                     ENDDO  
                     RMU(IU,IQ) = (rTWO - DELM0)*SUM  
                  ENDDO  
  
                  IF( FBEAM > rZERO ) THEN  
                     SUM  = rZERO  
                     DO K = MAZIM, NSTR - iONE  
                        SUM  = SUM + HLPR(K)*YLMU(K,IU)*YLM0(K)  
                     ENDDO  
                     RMU(IU,iZERO) = (rTWO - DELM0)*SUM  
                  END IF  
  
                  IF( MAZIM == iZERO ) THEN  
!                               ** Integrate bidirectional reflectivity  
!                                  at reflection angles UMU and  
!                                  incident angles GMU to get  
!                                  directional emissivity at  
!                                  user angles UMU.  
                     DREF  = rZERO  
                     DO JG = iONE, NMUG  
                        SUM  = rZERO  
                        DO K = iZERO, NSTR - iONE  
                           SUM = SUM + HLPR(K)*YLMU(K,IU)*YLMG(K,JG)  
                        ENDDO  
                        DREF  = DREF + rTWO*GWT(JG)*GMU(JG)*SUM  
                     ENDDO  
                     EMU(IU) = rONE - DREF  
                  END IF  
               END IF  
            END IF  
         ENDDO  
      END IF  
  
      END SUBROUTINE SURFAC  
  
      SUBROUTINE SOLVEC( AMB, APB, ARRAY, CMU, CWT, GL, &
           MAZIM, NN, NSTR, YLM0, YLMC, CC, &
           EVECC, EVAL, KK, GC, AAD, EVECCD, EVALD, &
           WK, WKD, DELM0, FBEAM, IPVT, PI, UMU0, ZJ, ZZ, &
           OPRIM, LC, mu2, glsave, dgl)  
  
!bm  SOLVEC calls SOLEIG and UPBEAM; if UPBEAM reports a potenially   
!bm  unstable solution, the calculation is repeated with a slightly   
!bm  changed single scattering albedo; this process is iterates   
!bm  until a stable solution is found; as stable solutions may be   
!bm  reached either by increasing or by decreasing the single   
!bm  scattering albedo, both directions are explored ('upward' and  
!bm  'downward' iteration); the solution which required the smaller   
!bm  change in the single scattering albedo is finally returned   
!bm  by SOLVEC.  
  
!gy added glsave and dgl to call to allow adjustable dimensioning  
  
!     .. Scalar Arguments ..  
  
      integer(ik), intent(in) :: MAZIM, NN, NSTR, LC  
      real(dk), intent(in)    :: DELM0, FBEAM, PI, UMU0, OPRIM  
      real(dk), intent(in)    :: mu2  
  
!     ..  
!     .. Array Arguments ..  
  
      integer(ik), intent(inout) ::   IPVT(:)  
        
      real(dk), intent(inout) :: WK(:)  
      real(dk), intent(in)    :: CMU(:), CWT(:), YLM0(0:), YLMC(0:,:)  
  
      real(dk), intent(inout) :: AMB(:,:), APB(:,:), ARRAY(:,:)  
      real(dk), intent(inout) :: EVAL(:), EVECC(:,:), KK(:), GC(:,:)  
      real(dk), intent(inout) :: GLSAVE(0:), DGL(0:), CC(:,:), ZJ(:), ZZ(:)  
      real(dk), intent(out)   :: GL(0:)  
  
      integer(ik) :: K  
      real(dk)    :: AAD(:,:), EVALD(:), EVECCD(:,:), WKD(:)  
  
!bm   Variables for instability fix  
        
      integer(ik) :: UAGAIN, DAGAIN  
      real(dk)    :: MINRCOND, ADD, UADD, DADD, SSA, DSSA, FACTOR  
        
      logical(lk) ::  DONE, NOUP, NODN, DEBUG, INSTAB  
        
!bm   reset parameters  
  
      DONE = .FALSE._lk
      NOUP = .FALSE._lk  
      NODN = .FALSE._lk  
  
  
!bm   flag for printing debugging output        
!      DEBUG  = .TRUE.  
      DEBUG  = .FALSE._lk  
  
!bm   instability parameter; the solution is considered   
!bm   unstable, if the RCOND reported by SGECO is smaller   
!bm   than MINRCOND  
!     MINRCOND = 5000._dk * R1MACH(4)  
      MINRCOND = 5000._dk * EPSILON( rONE )
  
!bm   if an instability is detected, the single scattering albedo  
!bm   is iterated downwards in steps of DADD and upwards in steps   
!bm   of UADD; in practice, MINRCOND and -MINRCOND should   
!bm   be reasonable choices for these parameters  
      DADD    = -MINRCOND  
      UADD    = MINRCOND  
  
      UAGAIN = iZERO
      DAGAIN = iZERO
      ADD    = DADD  
  
!bm   save array GL( ) because it will be   
!bm   changed if an iteration should be neccessary  
      DO K = MAZIM, NSTR - iONE
         GLSAVE(K) =  GL(K)  
      ENDDO  
        
      SSA = OPRIM  
!bm   in case of an instability reported by UPBEAM (INSTAB)  
!bm   the single scattering albedo will be changed by a small   
!bm   amount (ADD); this is indicated by DAGAIN or UAGAIN   
!bm   being larger than 0; a change in the single scattering   
!bm   albedo is equivalent to scaling the array GL( )  
  
 666  IF ( DAGAIN > iZERO .OR. UAGAIN > iZERO)  THEN  
         FACTOR = (SSA + ADD) / SSA  
         DO K = MAZIM, NSTR - iONE
            GL(K) =  GL(K) * FACTOR  
         ENDDO  
  
         SSA = SSA + ADD  
!bm   if the single scattering albedo is now smaller than 0  
!bm   the downward iteration is stopped and upward iteration   
!bm   is forced instead  
  
         IF( SSA < DITHER) THEN  
            NODN = .TRUE._lk
            DAGAIN = -iONE
            goto 778  
         ENDIF  
  
!bm   if the single scattering albedo is now larger than its maximum   
!bm   allowed value (1.0 - DITHER), the upward iteration is   
!bm   stopped and downward iteration is forced instead  
  
         IF( SSA > rONE - DITHER) THEN  
            NOUP = .TRUE._lk
            UAGAIN = -iONE
            goto 888  
         ENDIF  
      ENDIF  
  
!     ** Solve eigenfunction problem in Eq. STWJ(8B);  
!        return eigenvalues and eigenvectors  
 777     CALL SOLEIG( AMB, APB, ARRAY, CMU, CWT, GL, &
           MAZIM, NN, NSTR, YLMC, CC, EVECC, EVAL, &
           KK, GC, AAD, EVECCD, EVALD, WKD )  
  
!     ** Calculate particular solutions of  
!        q.SS(18) for incident beam source  
      IF ( FBEAM > rZERO ) THEN  
         CALL  UPBEAM( mu2, &
              ARRAY, CC, CMU, DELM0, FBEAM, GL, &
              IPVT, MAZIM, NN, NSTR, PI, UMU0, WK, &
              YLM0, YLMC, ZJ, ZZ, MINRCOND, INSTAB)  
      ENDIF  
        
!     ** Calculate particular solutions of  
!        Eq. SS(15) for thermal emission source  
!        (not available in psndo.f)  
!bm   finished if the result is stable on the first try  
      IF ( (.NOT. INSTAB) .AND. (UAGAIN == iZERO) .AND. (DAGAIN == iZERO)) THEN  
         goto 999  
      ENDIF  
  
!bm   downward iteration  
      IF( INSTAB .AND. UAGAIN == iZERO )  THEN  
         DAGAIN = DAGAIN + iONE
         GOTO 666  
      ENDIF  
        
!bm   upward iteration  
      IF( INSTAB .AND. UAGAIN > iZERO )  THEN  
         UAGAIN = UAGAIN + iONE
         GOTO 666  
      ENDIF  
  
  
!bm   ( DAGAIN .NE. 0 ) at this place means that the downward  
!bm   iteration is finished   
  
 778  IF (DAGAIN /= iZERO .AND. UAGAIN == iZERO) THEN  
!bm   save downward iteration data for later use and   
!bm   restore original input data  
         DO K = MAZIM, NSTR - iONE
            DGL(K) =  GL(K)  
            GL(K)  =  GLSAVE(K)  
         ENDDO  
  
         DSSA = SSA  
         SSA = OPRIM  
!bm   start upward iteration  
         ADD = UADD  
         UAGAIN = UAGAIN + iONE
         GOTO 666  
      ENDIF  
  
!bm   both iterations finished  
 888  IF (DONE) THEN  
         goto 998  
      ENDIF  
  
!bm  if neither upward nor downward iteration converged, the   
!bm  original conditions are restored and SOLEIG/UPBEAM   
!bm  is called for the last time   
      IF (NOUP .AND. NODN) THEN  
         DO K = MAZIM, NSTR - iONE
            GL(K) =  GLSAVE(K)  
         ENDDO  
           
         SSA = OPRIM  
           
         IF (DEBUG) THEN  
            write (*,*) '! *** Neither upward nor downward iteration'  
            write (*,*) '! *** converged; using original result.'  
         ENDIF  
  
         DONE = .TRUE._lk
         GOTO 777  
      ENDIF  
  
!bm  if upward iteration did not converge, the stable downward conditions  
!bm  are restored and SOLEIG/UPBEAM is called for the last time  
      IF (NOUP) THEN  
         DO K = MAZIM, NSTR - iONE
            GL(K) =  DGL(K)  
         ENDDO  
           
         SSA = DSSA  
           
         IF (DEBUG) THEN  
            write (*,*) '! *** The upward iteration did not converge.'  
            write (*,*) '! *** Had to iterate ', DAGAIN,' times in layer LC =', LC,';'  
            write (*,*) '! *** changed SSA from ',OPRIM,' to ', SSA,','  
            write (*,*) '! *** by a factor of ', SSA/OPRIM  
         ENDIF  
  
         DONE = .TRUE._lk
         GOTO 777  
      ENDIF  
  
!bm  if downward iteration did not converge, we are done   
!bm  (the result of the upward iteration will be used)  
      IF (NODN) THEN  
         IF (DEBUG) THEN  
            write (*,*) '! *** The downward iteration did not converge.'  
            write (*,*) '! *** Had to iterate ', UAGAIN,' times in layer LC =', LC,';'  
            write (*,*) '! *** changed SSA from ',OPRIM, ' to ', SSA,','  
            write (*,*) '! *** by a factor of ', SSA/OPRIM  
         ENDIF  
           
         DONE = .TRUE._lk
         GOTO 998  
      ENDIF  
        
!bm   if both iterations converged, and if the upward iteration   
!bm   required more steps than the downward iteration, the stable   
!bm   downward conditions are restored and SOLEIG/UPBEAM is   
!bm   called for the last time   
      IF (UAGAIN > DAGAIN) THEN  
         DO K = MAZIM, NSTR - iONE
            GL(K) =  DGL(K)  
         ENDDO  
           
         SSA = DSSA  
           
         IF (DEBUG) THEN  
            write (*,*) '! *** Both iterations converged; using downward.'  
            write (*,*) '! *** Had to iterate ',DAGAIN,' times in layer LC =', LC,';'  
            write (*,*) '! *** changed SSA from ',OPRIM, ' to ', SSA,','  
            write (*,*) '! *** by a factor of ', SSA/OPRIM  
         ENDIF  
  
         DONE = .TRUE._lk
         GOTO 777  
      ELSE  
         IF (DEBUG) THEN  
            write (*,*) '! *** Both iterations converged; using upward.'  
            write (*,*) '! *** Had to iterate ', UAGAIN,' times in layer LC =', LC,';'  
            write (*,*) '! *** changed SSA from ',OPRIM, ' to ', SSA,','  
            write (*,*) '! *** by a factor of ', SSA/OPRIM  
         ENDIF  
  
         DONE = .TRUE._lk
         goto 998  
      ENDIF  
        
!bm   finally restore original input data  
 998  DO K = MAZIM, NSTR - iONE
         GL(K) =  GLSAVE(K)  
      ENDDO  
        
 999  CONTINUE  
  
      END SUBROUTINE SOLVEC  
  
      SUBROUTINE UPBEAM( mu2, &
                         ARRAY, CC, CMU, DELM0, FBEAM, GL, IPVT, MAZIM, &
                         NN, NSTR, PI, UMU0, WK, YLM0, YLMC, ZJ, &
                         ZZ, MINRCOND, INSTAB )  
  
      use linalgebra, only : linalgebra_t  
  
!         Finds the incident-beam particular solution of SS(18)  
  
!   I N P U T    V A R I A B L E S:  
  
!       CC     :  C-sub-ij in Eq. SS(5)  
!       CMU    :  Abscissae for Gauss quadrature over angle cosine  
!       DELM0  :  Kronecker delta, delta-sub-m0  
!       GL     :  Delta-M scaled Legendre coefficients of phase function  
!                    (including factors 2L+1 and single-scatter albedo)  
!       MAZIM  :  Order of azimuthal component  
!       YLM0   :  Normalized associated Legendre polynomial  
!                    at the beam angle  
!       YLMC   :  Normalized associated Legendre polynomial  
!                    at the quadrature angles  
!       (remainder are DISORT input variables)  
  
!   O U T P U T    V A R I A B L E S:  
  
!       ZJ     :  Right-hand side vector X-sub-zero in SS(19); also the  
!                 solution vector Z-sub-zero after solving that system  
  
!       ZZ     :  Permanent storage for ZJ, but re-ordered  
  
!   I N T E R N A L    V A R I A B L E S:  
  
!       ARRAY  :  Coefficient matrix in left-hand side of Eq. SS(19)  
!       IPVT   :  Integer vector of pivot indices required by LINPACK  
!       WK     :  Scratch array required by LINPACK  
  
!   Called by- DISORT  
!   Calls- SGECO, ERRMSG, SGESL  
! +-------------------------------------------------------------------+  
  
  
!     .. Scalar Arguments ..  
  
      integer(ik), intent(in)  :: MAZIM, NN, NSTR  
      logical(lk), intent(out) :: INSTAB  
      real(dk), intent(in)     :: MINRCOND  
      real(dk), intent(in)     :: DELM0, FBEAM, PI, UMU0  
      real(dk), intent(in)     :: mu2  
!     ..  
!     .. Array Arguments ..  
  
      integer(ik), intent(inout) :: IPVT(:)  
  
      real(dk), intent(in) :: CC(:,:), CMU(:), &
                              GL(0:), YLM0(0:), YLMC(0:,:)  
      real(dk), intent(inout) :: WK(:)  
      real(dk), intent(out) :: ARRAY(:,:)  
      real(dk), intent(out) :: ZJ(:), ZZ(:)  
!     ..  
!     .. Local Scalars ..  
  
      integer(ik)   :: IQ, JOB, JQ, K  
      real(dk)      :: RCOND, SUM  
  
      TYPE(linalgebra_t) :: linpack  
  
      DO IQ = iONE, NSTR  
         DO JQ = iONE, NSTR  
            ARRAY(IQ,JQ) = -CC(IQ,JQ)  
         ENDDO  
  
         ARRAY(IQ,IQ) = rONE + CMU(IQ) / mu2 + ARRAY(IQ,IQ)  
         SUM  = rZERO  
         DO K = MAZIM, NSTR - iONE
            SUM  = SUM + GL( K )*YLMC(K,IQ)*YLM0(K)  
         ENDDO  
         ZJ(IQ) = (rTWO - DELM0)*FBEAM*SUM / FOURPI
      ENDDO  
  
!                  ** Find L-U (lower/upper triangular) decomposition  
!                     of ARRAY and see if it is nearly singular  
!                     (NOTE:  ARRAY is destroyed)  
      RCOND  = rZERO  
  
      CALL linpack%SGECO( ARRAY, NSTR, IPVT, RCOND, WK )  
  
!bm      IF( 1.0 + RCOND.EQ.1.0 )  
!bm     &    CALL ERRMSG('UPBEAM--SGECO says matrix near singular',.FALSE.)  
!bm  
!bm   replaced original check of RCOND by the following:  
      INSTAB = .FALSE._lk
      IF( ABS(RCOND) .LT. MINRCOND )  THEN  
         INSTAB = .TRUE._lk
         RETURN  
      ENDIF  
  
!                ** Solve linear system with coeff matrix ARRAY  
!                   (assumed already L-U decomposed) and R.H. side(s)  
!                   ZJ;  return solution(s) in ZJ  
      JOB  = iZERO
  
      CALL linpack%SGESL( ARRAY, NSTR, IPVT, ZJ, JOB )  
  
      DO IQ = iONE, NN  
         ZZ(IQ + NN)        = ZJ(IQ)  
         ZZ(NN + iONE - IQ) = ZJ(IQ + NN)  
      ENDDO  
  
      END SUBROUTINE UPBEAM  
  
      SUBROUTINE ZEROAL( EXPBEA, FLYR, OPRIM, TAUCPR, XR0, XR1, &
                         CMU, CWT, PSI, WK, Z0, Z1, ZJ, &
                         HLPR, YLM0, &
                         ARRAY, CC, EVECC, &
                         GL, YLMC, YLMU, &
                         KK, LL, ZZ, ZPLK0, ZPLK1, &
                         GC, LAYRU, UTAUPR, &
                         GU, Z0U, Z1U, ZBEAM, &
                         EVAL, AMB, APB, IPVT, Z, &
                         RFLDIR, RFLDN, FLUP, UAVG, DFDT, &
                         TRNMED, U0U, UU )  
  
!         ZERO ARRAYS  
  
!   Called by- DISORT  
! --------------------------------------------------------------------  
  
!     .. Array Arguments ..  
  
      integer(ik), intent(out) ::   IPVT(:), LAYRU(:)  
      real(dk), intent(out)    :: &
                AMB(:,:), APB(:,:), ARRAY(:,:), CC(:,:), &
                CMU(:), CWT(:), DFDT(:), EVAL(:), EVECC(:,:), &
                EXPBEA(:), FLUP(:), FLYR(:), GC(:,:,:), GL(:,:), &
                GU(:,:,:), HLPR(:), KK(:,:), LL(:,:), OPRIM(:), &
                PSI(:), RFLDIR(:), RFLDN(:), TAUCPR(:), &
                TRNMED(:), U0U(:,:), UAVG(:), UTAUPR(:),& 
                UU(:,:,:), &
                WK(:), XR0(:), XR1(:), YLM0(:), YLMC(:,:), &
                YLMU(:,:), Z(:), Z0(:), Z0U(:,:), Z1(:), Z1U(:,:), &
                ZBEAM(:,:), ZJ(:), ZPLK0(:,:), ZPLK1(:,:), ZZ(:,:)  
  
!     ..  
  
      EXPBEA = rZERO  
      FLYR   = rZERO  
      OPRIM  = rZERO  
      TAUCPR = rZERO  
      XR0    = rZERO  
      XR1    = rZERO  
  
         CMU = rZERO  
         CWT = rZERO  
         PSI = rZERO  
         WK  = rZERO  
         Z0  = rZERO  
         Z1  = rZERO  
         ZJ  = rZERO  
  
         HLPR = rZERO  
         YLM0 = rZERO  
  
         ARRAY = rZERO  
         CC    = rZERO  
         EVECC = rZERO  
  
         GL = rZERO  
  
         YLMC = rZERO  
  
         YLMU = rZERO  
  
         KK    = rZERO  
         LL    = rZERO  
         ZZ    = rZERO  
         ZPLK0 = rZERO  
         ZPLK1 = rZERO  
  
         GC = rZERO  
  
         LAYRU  = iZERO
         UTAUPR = rZERO  
  
         GU = rZERO  
  
         Z0U   = rZERO  
         Z1U   = rZERO  
         ZBEAM = rZERO  
  
         EVAL = rZERO  
  
         AMB = rZERO  
         APB = rZERO  
  
         IPVT = iZERO
         Z    = rZERO  
  
         RFLDIR = rZERO  
         RFLDN  = rZERO  
         FLUP   = rZERO  
         UAVG   = rZERO  
         DFDT   = rZERO  
  
         TRNMED = rZERO  
  
         U0U = rZERO  
  
         UU = rZERO  
  
      END SUBROUTINE ZEROAL  
  
      real(dk) FUNCTION DREF( MU, HL, NSTR )  
  
!        Exact flux albedo for given angle of incidence, given  
!        a bidirectional reflectivity characterized by its  
!        Legendre coefficients ( NOTE** these will only agree  
!        with bottom-boundary albedos calculated by DISORT in  
!        the limit as number of streams go to infinity, because  
!        DISORT evaluates the integral 'CL' only approximately,  
!        by quadrature, while this routine calculates it exactly.)  
  
!  INPUT :   MU     Cosine of incidence angle  
!            HL     Legendre coefficients of bidirectional reflectivity  
!          NSTR     Number of elements of HL to consider  
  
!  INTERNAL VARIABLES (P-sub-L is the L-th Legendre polynomial) :  
  
!       CL      Integral from 0 to 1 of  MU * P-sub-L(MU)  
!                   (vanishes for  L = 3, 5, 7, ... )  
!       PL      P-sub-L  
!       PLM1    P-sub-(L-1)  
!       PLM2    P-sub-(L-2)  
  
!   Called by- CHEKIN  
!   Calls- ERRMSG  
! +-------------------------------------------------------------------+  
  
!     .. Parameters ..  
  
      integer(ik), PARAMETER ::  MAXTRM = 100_ik
!     ..  
!     .. Scalar Arguments ..  
  
      integer(ik), intent(in) :: NSTR  
      real(dk), intent(in)    :: MU  
!     ..  
!     .. Array Arguments ..  
  
      real(dk), intent(in)   :: HL(0:NSTR)  
!     ..  
!     .. Local Scalars ..  
  
      logical(lk)   :: PASS1  
      integer(ik)   :: L  
      real(dk)      :: CL, PL, PLM1, PLM2  
!     ..  
!     .. Local Arrays ..  
  
      real(dk)      :: C(MAXTRM)  
!     ..  
  
      SAVE      PASS1, C  
      DATA      PASS1 / .TRUE._lk /  
!     ..  
  
  
      IF( PASS1 ) THEN  
         PASS1  = .FALSE._lk  
         CL     = 0.125_dk
         C(2)   = 10._dk*CL  
  
         DO L = 4_ik, MAXTRM, iTWO
            CL   = -CL*real(L - iTHREE,dk) / real(L + iTWO,dk)
            C(L) = rTWO*real(iTWO*L + iONE,dk)*CL  
         ENDDO  
      END IF  
  
      IF( NSTR < iTWO .OR. ABS(MU) > rONE ) &
          CALL ERRMSG( 'DREF--input argument error(s)',.True. )  
  
      IF( NSTR > MAXTRM ) &
          CALL ERRMSG( 'DREF--parameter MAXTRM too small',.True. )  
  
  
      DREF  = HL(iZERO) - rTWO*HL(iONE)*MU  
      PLM2  = rONE  
      PLM1  = -MU  
  
      DO L = iTWO, NSTR - iONE
!                                ** Legendre polynomial recurrence  
         PL = (real(iTWO*L - iONE,dk)*(-MU)*PLM1 - real(L-iONE,dk)*PLM2) / real(L,dk)
         IF( MOD( L,iTWO ) == iZERO ) DREF = DREF + C(L)*HL(L)*PL  
         PLM2  = PLM1  
         PLM1  = PL  
      ENDDO  
  
      IF( DREF < rZERO .OR. DREF > rONE ) &
          CALL ERRMSG( 'DREF--albedo value not in (0,1)',.False. )  
  
      END FUNCTION DREF  
  
      real(dk) FUNCTION RATIO( A, B )  
  
!        Calculate ratio  A/B  with over- and under-flow protection  
!        (thanks to Prof. Jeff Dozier for some suggestions here).  
!        Since this routine takes two logs, it is no speed demon,  
!        but it is invaluable for comparing results from two runs  
!        of a program under development.  
  
!        NOTE:  In Fortran90, built-in functions TINY and HUGE  
!               can replace the R1MACH calls.  
! ---------------------------------------------------------------  
  
!     .. Scalar Arguments ..  
  
      real(dk), intent(in) :: A, B  
!     ..  
!     .. Local Scalars ..  
  
      logical(lk)   PASS1  
      real(dk)      ABSA, ABSB, POWA, POWB, POWMAX, POWMIN
!     ..  
!     .. External Functions ..  
  
!     EXTERNAL  R1MACH  
!     ..  
!     .. Intrinsic Functions ..  
  
      INTRINSIC ABS, LOG10, SIGN  
!     ..  
      SAVE      PASS1, POWMAX, POWMIN  
      DATA      PASS1 / .TRUE._lk /  
!     ..  
  
  
      IF( PASS1 ) THEN  
!        TINY   = R1MACH( 1 )  
!        HUGE   = R1MACH( 2 )  
         POWMAX = LOG10( HUGE(rONE) )  
         POWMIN = LOG10( TINY(rONE) )  
         PASS1  = .FALSE._lk
      END IF  
  
  
      IF( A == rZERO ) THEN  
         IF( B == rZERO ) THEN  
            RATIO  = rONE  
         ELSE  
            RATIO  = rZERO  
         END IF  
      ELSE IF( B == rZERO ) THEN  
         RATIO  = SIGN( HUGE(rONE), A )  
      ELSE  
         ABSA   = ABS( A )  
         ABSB   = ABS( B )  
         POWA   = LOG10( ABSA )  
         POWB   = LOG10( ABSB )  
  
         IF( ABSA < TINY(rONE) .AND. ABSB < TINY(rONE) ) THEN  
            RATIO  = rONE  
         ELSE IF( POWA - POWB >= POWMAX ) THEN  
            RATIO  = HUGE( rONE )  
         ELSE IF( POWA - POWB <= POWMIN ) THEN  
            RATIO  = TINY( rONE )
         ELSE  
            RATIO  = ABSA / ABSB  
         END IF  
!                      ** DONT use old trick of determining sign  
!                      ** from A*B because A*B may (over/under)flow  
  
         IF( ( A > rZERO .AND. B < rZERO ) .OR. &
             ( A < rZERO .AND. B > rZERO ) ) RATIO = -RATIO  
      END IF  
  
      END FUNCTION RATIO  
  
      SUBROUTINE  ErrMsg( MESSAG, FATAL )  
  
!        Print out a warning or error message;  abort if error  
!        after making symbolic dump (machine-specific)  
  
      logical(lk)               :: FATAL, MsgLim, Cray  
      CHARACTER*(*), intent(in) :: MESSAG  

      integer(ik)   :: MaxMsg, NumMsg  
      SAVE          MaxMsg, NumMsg, MsgLim  
      DATA NumMsg / iZERO /,  MaxMsg / 100_ik /,  MsgLim / .FALSE._lk /  
  
      IF ( FATAL )  THEN  
         WRITE ( *, '(//,2A,//)' )  ' ******* ERROR >>>>>>  ', MESSAG  
         STOP  
      END IF  
  
      NumMsg = NumMsg + 1  
      IF( MsgLim )  RETURN  
  
      IF ( NumMsg.LE.MaxMsg )  THEN  
         WRITE ( *, '(/,2A,/)' )  ' ******* WARNING >>>>>>  ', MESSAG  
      ELSE  
         WRITE ( *,99 )  
         MsgLim = .True.  
      ENDIF  
  
   99 FORMAT( //,' >>>>>>  TOO MANY WARNING MESSAGES --  ',  
         'They will no longer be printed  <<<<<<<', // )  
  
      END SUBROUTINE ErrMsg  
  
      logical(lk) FUNCTION  WrtBad ( VarNam )  
  
!          Write names of erroneous variables and return 'TRUE'  
  
!      INPUT :   VarNam = Name of erroneous variable to be written  
!                         ( CHARACTER, any length )  
  
      CHARACTER*(*)  VarNam  
      integer(ik)        MaxMsg, NumMsg  
      SAVE  NumMsg, MaxMsg  
      DATA  NumMsg / 0 /,  MaxMsg / 50 /  
  
  
      WrtBad = .TRUE.  
      NumMsg = NumMsg + 1  
      WRITE ( *, '(3A)' )  ' ****  Input variable  ', VarNam,  
                           '  in error  ****'  
      IF ( NumMsg.EQ.MaxMsg )  
         CALL  ErrMsg ( 'Too many input errors.  Aborting...', .TRUE. )  
  
      END FUNCTION WrtBad  
  
      logical(lk) FUNCTION  WrtDim ( DimNam, MinVal )  
  
!          Write name of too-small symbolic dimension and  
!          the value it should be increased to;  return 'TRUE'  
  
!      INPUT :  DimNam = Name of symbolic dimension which is too small  
!                        ( CHARACTER, any length )  
!               Minval = Value to which that dimension should be  
!                        increased (at least)  
  
      CHARACTER*(*)  DimNam  
      integer(ik)        MinVal  
  
  
      WRITE ( *, '(3A,I7)' )  ' ****  Symbolic dimension  ', DimNam,  
                           '  should be increased to at least ', MinVal  
      WrtDim = .TRUE.  
  
      END FUNCTION WrtDim  
  
      logical(lk) FUNCTION  TstBad( VarNam, RelErr )  
  
!       Write name (VarNam) of variable failing self-test and its  
!       percent error from the correct value;  return  'FALSE'.  
  
      CHARACTER*(*)  VarNam  
      real(dk)           RelErr  
  
  
      TstBad = .FALSE.  
      WRITE( *, '(/,3A,1P,E11.2,A)' )  
             ' Output variable ', VarNam,' differed by ', 100.*RelErr,  
             ' per cent from correct value.  Self-test failed.'  
  
      END FUNCTION TstBad  
  
      FUNCTION D1MACH(i)  
!-----------------------------------------------------------------------------*  
!= PURPOSE:                                                                  =*  
!= D1MACH calculates various machine constants in single precision.          =*  
!-----------------------------------------------------------------------------*  
!= PARAMETERS:                                                               =*  
!=   I       -  integer(ik), identifies the machine constant (0<I<5)         (I) =*  
!=   D1MACH  -  real(dk), machine constant in single precision               (O) =*  
!=      I=1     - the smallest non-vanishing normalized floating-point       =*  
!=                power of the radix, i.e., D1MACH=FLOAT(IBETA)**MINEXP      =*  
!=      I=2     - the largest finite floating-point number.  In              =*  
!=                particular D1MACH=(1.0-EPSNEG)*FLOAT(IBETA)**MAXEXP        =*  
!=                Note - on some machines D1MACH will be only the            =*  
!=                second, or perhaps third, largest number, being            =*  
!=                too small by 1 or 2 units in the last digit of             =*  
!=                the significand.                                           =*  
!=      I=3     - A small positive floating-point number such that           =*  
!=                1.0-D1MACH .NE. 1.0. In particular, if IBETA = 2           =*  
!=                or  IRND = 0, D1MACH = FLOAT(IBETA)**NEGEPS.               =*  
!=                Otherwise,  D1MACH = (IBETA**NEGEPS)/2.  Because           =*  
!=                NEGEPS is bounded below by -(IT+3), D1MACH may not         =*  
!=                be the smallest number that can alter 1.0 by               =*  
!=                subtraction.                                               =*  
!=      I=4     - the smallest positive floating-point number such           =*  
!=                that  1.0+D1MACH .NE. 1.0. In particular, if either        =*  
!=                IBETA = 2  or  IRND = 0, D1MACH=FLOAT(IBETA)**MACHEP.      =*  
!=                Otherwise, D1MACH=(FLOAT(IBETA)**MACHEP)/2                 =*  
!=  (see routine T665D for more information on different constants)          =*  
!-----------------------------------------------------------------------------*  
  
      real(dk) :: d1mach  
      integer  :: i  
     
      logical(lk) :: doinit
      DATA doinit /.TRUE._lk/  
      SAVE doinit  
  
      real(dk) :: dmach(4)   
      SAVE dmach  
  
      IF ( i .GE. 1 .AND. i .LE. 4  ) THEN  
! compute constants at first call only  
        IF (doinit) THEN  
           CALL T665d(dmach)  
           doinit = .FALSE._lk
        ENDIF  
        d1mach = dmach(i)  
      ELSE  
        WRITE(0,*) '>>> ERROR (D1MACH) <<<  invalid argument'  
        STOP  
      ENDIF  
  
!!csm  
!!!! over-ride by sm on 5/26/03.  For some compilers than don't allow  
! calculation of d1mach(4).  Use value found on ACD server.  
  
!      if( i .eq. 4) d1mach = 2.22e-15  
  
      END FUNCTION D1MACH  
  
!      ALGORITHM 665, COLLECTED ALGORITHMS FROM ACM.  
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,  
!      VOL. 14, NO. 4, PP. 303-311.  
      SUBROUTINE T665D(DMACH)  
!-----------------------------------------------------------------------  
! This subroutine is a double precision version of subroutine T665R.  
! See code of T665R for detailed comments and explanation  
!-----------------------------------------------------------------------  
      real(dk)    :: DMACH(4)  

      integer(ik) :: I,IBETA,IEXP,IRND,IT,ITEMP,IZ,J,K,MACHEP,MAXEXP, &
                     MINEXP,MX,NEGEP,NGRD,NXRES  
      real(dk)    :: A,B,BETA,BETAIN,BETAH,CONV,EPS,EPSNEG,ONE, &
                     T,TEMP,TEMPA,TEMP1,TWO,XMAX,XMIN,Y,Z,ZERO  
!-----------------------------------------------------------------------  
      CONV(I) = real(I,dk)  
!-----------------------------------------------------------------------  
!  Determine IBETA, BETA ala Malcolm.  
!-----------------------------------------------------------------------  
      A = rONE  
      DO
        A = A + A  
        TEMP = A + rONE  
        TEMP1 = TEMP - A  
        IF (TEMP1 - rONE /= rZERO) EXIT
      ENDDO

      B = rONE  
      DO
        B = B + B  
        TEMP = A + B  
        ITEMP = INT(TEMP-A,ik)  
        IF (ITEMP /= iZERO) EXIT
      ENDDO
      IBETA = ITEMP  
      BETA = CONV(IBETA)  
!-----------------------------------------------------------------------  
!  Determine IT, IRND.  
!-----------------------------------------------------------------------  
      IT = iZERO
      B = rONE  
      DO
        IT = IT + iONE
        B = B * BETA  
        TEMP = B + rONE  
        TEMP1 = TEMP - B  
        IF (TEMP1-rONE /= rZERO) EXIT
      ENDDO

      IRND = iZERO  
      BETAH = BETA / rTWO  
      TEMP = A + BETAH  
      IF (TEMP-A /= ZERO) IRND = riONE
      TEMPA = A + BETA  
      TEMP = TEMPA + BETAH  
      IF ((IRND == iZERO) .AND. (TEMP-TEMPA /= rZERO)) IRND = iTWO
!-----------------------------------------------------------------------  
!  Determine NEGEP, EPSNEG.  
!-----------------------------------------------------------------------  
      NEGEP = IT + iTHREE  
      BETAIN = rONE / BETA  
      A = rONE  
      DO I = iONE, NEGEP  
        A = A * BETAIN  
      ENDDO

      B = A  
      DO
        TEMP = rONE - A  
        IF (TEMP-rONE /= rZERO) EXIT
        A = A * BETA  
        NEGEP = NEGEP - iONE
      ENDDO

      NEGEP = -NEGEP  
      EPSNEG = A  
      IF ( IBETA /= iTWO .AND. IRND /= iZERO ) THEN
        A = (A*(rONE+A)) / rTWO  
        TEMP = rONE - A  
        IF (TEMP-rONE /= ZERO) EPSNEG = A  
      ENDIF
!-----------------------------------------------------------------------  
!  Determine MACHEP, EPS.  
!-----------------------------------------------------------------------  
      MACHEP = -(IT + iTHREE)
      A = B  
      DO
        TEMP = rONE+A  
        IF (TEMP-rONE /= rZERO) EXIT
        A = A * BETA  
        MACHEP = MACHEP + 1  
      ENDDO

      EPS = A  
      TEMP = TEMPA+BETA*(ONE+EPS)  
      IF ( IBETA /= iTWO .AND. IRND /= iZERO ) THEN
        A = (A*(rONE + A)) / rTWO  
        TEMP = rONE + A  
        IF (TEMP-rONE /= rZERO) EPS = A  
      ENDIF
!-----------------------------------------------------------------------  
!  Determine NGRD.  
!-----------------------------------------------------------------------  
      NGRD = iZERO
      TEMP = rONE + EPS  
      IF ( IRND == 0 .AND. (TEMP*ONE-ONE /= ZERO)) NGRD = iONE
!-----------------------------------------------------------------------  
!  Determine IEXP, MINEXP, XMIN.  
!  
!  Loop to determine largest I and K = 2**I such that  
!         (1/BETA) ** (2**(I))  
!  does not underflow.  
!  Exit from loop is signaled by an underflow.  
!-----------------------------------------------------------------------  
      I = iZERO
      K = iONE
      Z = BETAIN  
      T = rONE + EPS  
      NXRES = iZERO

      DO
        Y = Z  
        Z = Y * Y  
!-----------------------------------------------------------------------  
!  Check for underflow here.  
!-----------------------------------------------------------------------  
        A = Z * rONE  
        TEMP = Z * T  
        IF ( A+A == ZERO .OR. ABS(Z) >= Y ) EXIT
        TEMP1 = TEMP * BETAIN  
        IF (TEMP1*BETA == Z) EXIT
        I = I + iONE
        K = K + K  
      ENDDO

      IF (IBETA /= 10_ik) THEN
        IEXP = I + iONE
        MX = K + K  
      ELSE
!-----------------------------------------------------------------------  
!  This segment is for decimal machines only.  
!-----------------------------------------------------------------------  
        IEXP = iTWO
        IZ = IBETA  
        DO
          IF (K < IZ) EXIT
          IZ = IZ * IBETA  
          IEXP = IEXP + iONE
        ENDDO
        MX = IZ + IZ - iONE
      ENDIF
!-----------------------------------------------------------------------  
!  Loop to determine MINEXP, XMIN.  
!  Exit from loop is signaled by an underflow.  
!-----------------------------------------------------------------------  
      DO
        XMIN = Y  
        Y = Y * BETAIN  
!-----------------------------------------------------------------------  
!  Check for underflow here.  
!-----------------------------------------------------------------------  
        A = Y * rONE  
        TEMP = Y * T  
        IF (((A+A) == ZERO) .OR. ABS(Y) >= XMIN) EXIT
        K = K + iONE
        TEMP1 = TEMP * BETAIN  
        IF (TEMP1*BETA == Y) THEN
          NXRES = iTHREE
          XMIN = Y  
          EXIT
        ENDIF
      ENDDO

      MINEXP = -K  
!-----------------------------------------------------------------------  
!  Determine MAXEXP, XMAX.  
!-----------------------------------------------------------------------  
      IF (MX <= K+K-iTHREE .AND. IBETA /= 10_ik ) THEN
        MX = MX + MX  
        IEXP = IEXP + iONE
      ENDIF
      MAXEXP = MX + MINEXP  
!-----------------------------------------------------------------  
!  Adjust IRND to reflect partial underflow.  
!-----------------------------------------------------------------  
      IRND = IRND + NXRES  
!-----------------------------------------------------------------  
!  Adjust for IEEE-style machines.  
!-----------------------------------------------------------------  
      IF ((IRND == iTWO) .OR. (IRND == 5_ik)) MAXEXP = MAXEXP - iTWO
!-----------------------------------------------------------------  
!  Adjust for non-IEEE machines with partial underflow.  
!-----------------------------------------------------------------  
      IF ((IRND == iTHREE) .OR. (IRND == 4_ik)) MAXEXP = MAXEXP - IT  
!-----------------------------------------------------------------  
!  Adjust for machines with implicit leading bit in binary  
!  significand, and machines with radix point at extreme  
!  right of significand.  
!-----------------------------------------------------------------  
      I = MAXEXP + MINEXP  
      IF ((IBETA == iTWO) .AND. (I == iZERO)) MAXEXP = MAXEXP - iONE
      IF (I > 20_ik) MAXEXP = MAXEXP - iONE
      IF (A /= Y) MAXEXP = MAXEXP - iTWO
      XMAX = rONE - EPSNEG  
      IF (XMAX*rONE /= XMAX) XMAX = rONE - BETA * EPSNEG  
      XMAX = XMAX / (BETA * BETA * BETA * XMIN)  
      I = MAXEXP + MINEXP + iTHREE

      IF (I > iZERO) THEN
        DO J = iONE, I  
          IF (IBETA == iTWO) THEN
            XMAX = XMAX + XMAX  
          ELSE
            XMAX = XMAX * BETA  
          ENDIF
        ENDDO
        DMACH(1) = XMIN  
        DMACH(2) = XMAX  
        DMACH(3) = EPSNEG  
        DMACH(4) = EPS  
      ENDIF
  
      END SUBROUTINE T665D  
  
      FUNCTION R1MACH(i)  
  
!-----------------------------------------------------------------------------*  
!= PURPOSE:                                                                  =*  
!= R1MACH calculates various machine constants in single precision.          =*  
!-----------------------------------------------------------------------------*  
!= PARAMETERS:                                                               =*  
!=   I       -  integer(ik), identifies the machine constant (0<I<5)         (I) =*  
!=   R1MACH  -  real(dk), machine constant in single precision               (O) =*  
!=      I=1     - the smallest non-vanishing normalized floating-point       =*  
!=                power of the radix, i.e., R1MACH=FLOAT(IBETA)**MINEXP      =*  
!=      I=2     - the largest finite floating-point number.  In              =*  
!=                particular R1MACH=(1.0-EPSNEG)*FLOAT(IBETA)**MAXEXP        =*  
!=                Note - on some machines R1MACH will be only the            =*  
!=                second, or perhaps third, largest number, being            =*  
!=                too small by 1 or 2 units in the last digit of             =*  
!=                the significand.                                           =*  
!=      I=3     - A small positive floating-point number such that           =*  
!=                1.0-R1MACH .NE. 1.0. In particular, if IBETA = 2           =*  
!=                or  IRND = 0, R1MACH = FLOAT(IBETA)**NEGEPS.               =*  
!=                Otherwise,  R1MACH = (IBETA**NEGEPS)/2.  Because           =*  
!=                NEGEPS is bounded below by -(IT+3), R1MACH may not         =*  
!=                be the smallest number that can alter 1.0 by               =*  
!=                subtraction.                                               =*  
!=      I=4     - the smallest positive floating-point number such           =*  
!=                that  1.0+R1MACH .NE. 1.0. In particular, if either        =*  
!=                IBETA = 2  or  IRND = 0, R1MACH=FLOAT(IBETA)**MACHEP.      =*  
!=                Otherwise, R1MACH=(FLOAT(IBETA)**MACHEP)/2                 =*  
!=  (see routine T665R for more information on different constants)          =*  
!-----------------------------------------------------------------------------*  
  
      real(dk) r1mach  
      integer(ik) i  
     
      logical(lk) doinit  
      DATA doinit/.TRUE./  
      SAVE doinit  
  
      real(dk) rmach(4)   
      SAVE rmach  
  
      IF (( i .GE. 1 ) .AND. ( i .LE. 4 )) THEN  
! compute constants at first call only  
        IF (doinit) THEN  
           CALL t665r(rmach)  
           doinit = .FALSE.  
        ENDIF  
        r1mach = rmach(i)  
      ELSE  
        WRITE(0,*) '>>> ERROR (R1MACH) <<<  invalid argument'  
        STOP  
      ENDIF  
  
      END FUNCTION  
  
  
!      ALGORITHM 665, COLLECTED ALGORITHMS FROM ACM.  
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,  
!      VOL. 14, NO. 4, PP. 303-311.  
      SUBROUTINE T665R(RMACH)  
!-----------------------------------------------------------------------  
!  This Fortran 77 subroutine is intended to determine the parameters  
!   of the floating-point arithmetic system specified below.  The  
!   determination of the first three uses an extension of an algorithm  
!   due to M. Malcolm, CACM 15 (1972), pp. 949-951, incorporating some,  
!   but not all, of the improvements suggested by M. Gentleman and S.  
!   Marovich, CACM 17 (1974), pp. 276-277.  An earlier version of this  
!   program was published in the book Software Manual for the  
!   Elementary Functions by W. J. Cody and W. Waite, Prentice-Hall,  
!   Englewood Cliffs, NJ, 1980.  
!  
!  The program as given here must be modified before compiling.  If  
!   a single (double) precision version is desired, change all  
!   occurrences of CS (CD) in columns 1 and 2 to blanks.  
!  
!  Parameter values reported are as follows:  
!  
!       IBETA   - the radix for the floating-point representation  
!       IT      - the number of base IBETA digits in the floating-point  
!                 significand  
!       IRND    - 0 if floating-point addition chops  
!                 1 if floating-point addition rounds, but not in the  
!                   IEEE style  
!                 2 if floating-point addition rounds in the IEEE style  
!                 3 if floating-point addition chops, and there is  
!                   partial underflow  
!                 4 if floating-point addition rounds, but not in the  
!                   IEEE style, and there is partial underflow  
!                 5 if floating-point addition rounds in the IEEE style,  
!                   and there is partial underflow  
!       NGRD    - the number of guard digits for multiplication with  
!                 truncating arithmetic.  It is  
!                 0 if floating-point arithmetic rounds, or if it  
!                   truncates and only  IT  base  IBETA digits  
!                   participate in the post-normalization shift of the  
!                   floating-point significand in multiplication;  
!                 1 if floating-point arithmetic truncates and more  
!                   than  IT  base  IBETA  digits participate in the  
!                   post-normalization shift of the floating-point  
!                   significand in multiplication.  
!       MACHEP  - the largest negative integer such that  
!                 1.0+FLOAT(IBETA)**MACHEP .NE. 1.0, except that  
!                 MACHEP is bounded below by  -(IT+3)  
!       NEGEPS  - the largest negative integer such that  
!                 1.0-FLOAT(IBETA)**NEGEPS .NE. 1.0, except that  
!                 NEGEPS is bounded below by  -(IT+3)  
!       IEXP    - the number of bits (decimal places if IBETA = 10)  
!                 reserved for the representation of the exponent  
!                 (including the bias or sign) of a floating-point  
!                 number  
!       MINEXP  - the largest in magnitude negative integer such that  
!                 FLOAT(IBETA)**MINEXP is positive and normalized  
!       MAXEXP  - the smallest positive power of  BETA  that overflows  
!       EPS     - the smallest positive floating-point number such  
!                 that  1.0+EPS .NE. 1.0. In particular, if either  
!                 IBETA = 2  or  IRND = 0, EPS = FLOAT(IBETA)**MACHEP.  
!                 Otherwise,  EPS = (FLOAT(IBETA)**MACHEP)/2  
!       EPSNEG  - A small positive floating-point number such that  
!                 1.0-EPSNEG .NE. 1.0. In particular, if IBETA = 2  
!                 or  IRND = 0, EPSNEG = FLOAT(IBETA)**NEGEPS.  
!                 Otherwise,  EPSNEG = (IBETA**NEGEPS)/2.  Because  
!                 NEGEPS is bounded below by -(IT+3), EPSNEG may not  
!                 be the smallest number that can alter 1.0 by  
!                 subtraction.  
!       XMIN    - the smallest non-vanishing normalized floating-point  
!                 power of the radix, i.e.,  XMIN = FLOAT(IBETA)**MINEXP  
!       XMAX    - the largest finite floating-point number.  In  
!                 particular  XMAX = (1.0-EPSNEG)*FLOAT(IBETA)**MAXEXP  
!                 Note - on some machines  XMAX  will be only the  
!                 second, or perhaps third, largest number, being  
!                 too small by 1 or 2 units in the last digit of  
!                 the significand.  
!  
!     Latest revision - April 20, 1987  
!  
!     Author - W. J. Cody  
!              Argonne National Laboratory  
!  
!-----------------------------------------------------------------------  
      real(dk) rmach(4)  
      integer(ik) I,IBETA,IEXP,IRND,IT,ITEMP,IZ,J,K,MACHEP,MAXEXP,  
     1        MINEXP,MX,NEGEP,NGRD,NXRES  
      real(dk) A,B,BETA,BETAIN,BETAH,CONV,EPS,EPSNEG,ONE,T,TEMP,TEMPA,  
     1     TEMP1,TWO,XMAX,XMIN,Y,Z,ZERO  
!D    real(dk)(8) A,B,BETA,BETAIN,BETAH,CONV,EPS,EPSNEG,ONE,  
!D   1                 T,TEMP,TEMPA,TEMP1,TWO,XMAX,XMIN,Y,Z,ZERO  
!-----------------------------------------------------------------------  
      CONV(I) = real(dk)(I)  
!D    CONV(I) = DBLE(I)  
      ONE = CONV(1)  
      TWO = ONE + ONE  
      ZERO = ONE - ONE  
!-----------------------------------------------------------------------  
!  Determine IBETA, BETA ala Malcolm.  
!-----------------------------------------------------------------------  
      A = ONE  
   10 A = A + A  
         TEMP = A+ONE  
         TEMP1 = TEMP-A  
         IF (TEMP1-ONE .EQ. ZERO) GO TO 10  
      B = ONE  
   20 B = B + B  
         TEMP = A+B  
         ITEMP = INT(TEMP-A)  
         IF (ITEMP .EQ. 0) GO TO 20  
      IBETA = ITEMP  
      BETA = CONV(IBETA)  
!-----------------------------------------------------------------------  
!  Determine IT, IRND.  
!-----------------------------------------------------------------------  
      IT = 0  
      B = ONE  
  100 IT = IT + 1  
         B = B * BETA  
         TEMP = B+ONE  
         TEMP1 = TEMP-B  
         IF (TEMP1-ONE .EQ. ZERO) GO TO 100  
      IRND = 0  
      BETAH = BETA / TWO  
      TEMP = A+BETAH  
      IF (TEMP-A .NE. ZERO) IRND = 1  
      TEMPA = A + BETA  
      TEMP = TEMPA+BETAH  
      IF ((IRND .EQ. 0) .AND. (TEMP-TEMPA .NE. ZERO)) IRND = 2  
!-----------------------------------------------------------------------  
!  Determine NEGEP, EPSNEG.  
!-----------------------------------------------------------------------  
      NEGEP = IT + 3  
      BETAIN = ONE / BETA  
      A = ONE  
      DO 200 I = 1, NEGEP  
         A = A * BETAIN  
  200 CONTINUE  
      B = A  
  210 TEMP = ONE-A  
         IF (TEMP-ONE .NE. ZERO) GO TO 220  
         A = A * BETA  
         NEGEP = NEGEP - 1  
      GO TO 210  
  220 NEGEP = -NEGEP  
      EPSNEG = A  
      IF ((IBETA .EQ. 2) .OR. (IRND .EQ. 0)) GO TO 300  
      A = (A*(ONE+A)) / TWO  
      TEMP = ONE-A  
      IF (TEMP-ONE .NE. ZERO) EPSNEG = A  
!-----------------------------------------------------------------------  
!  Determine MACHEP, EPS.  
!-----------------------------------------------------------------------  
  300 MACHEP = -IT - 3  
      A = B  
  310 TEMP = ONE+A  
         IF (TEMP-ONE .NE. ZERO) GO TO 320  
         A = A * BETA  
         MACHEP = MACHEP + 1  
      GO TO 310  
  320 EPS = A  
      TEMP = TEMPA+BETA*(ONE+EPS)  
      IF ((IBETA .EQ. 2) .OR. (IRND .EQ. 0)) GO TO 350  
      A = (A*(ONE+A)) / TWO  
      TEMP = ONE+A  
      IF (TEMP-ONE .NE. ZERO) EPS = A  
!-----------------------------------------------------------------------  
!  Determine NGRD.  
!-----------------------------------------------------------------------  
  350 NGRD = 0  
      TEMP = ONE+EPS  
      IF ((IRND .EQ. 0) .AND. (TEMP*ONE-ONE .NE. ZERO)) NGRD = 1  
!-----------------------------------------------------------------------  
!  Determine IEXP, MINEXP, XMIN.  
!  
!  Loop to determine largest I and K = 2**I such that  
!         (1/BETA) ** (2**(I))  
!  does not underflow.  
!  Exit from loop is signaled by an underflow.  
!-----------------------------------------------------------------------  
      I = 0  
      K = 1  
      Z = BETAIN  
      T = ONE + EPS  
      NXRES = 0  
  400 Y = Z  
         Z = Y * Y  
!-----------------------------------------------------------------------  
!  Check for underflow here.  
!-----------------------------------------------------------------------  
         A = Z * ONE  
         TEMP = Z * T  
         IF ((A+A .EQ. ZERO) .OR. (ABS(Z) .GE. Y)) GO TO 410  
         TEMP1 = TEMP * BETAIN  
         IF (TEMP1*BETA .EQ. Z) GO TO 410  
         I = I + 1  
         K = K + K  
      GO TO 400  
  410 IF (IBETA .EQ. 10) GO TO 420  
      IEXP = I + 1  
      MX = K + K  
      GO TO 450  
!-----------------------------------------------------------------------  
!  This segment is for decimal machines only.  
!-----------------------------------------------------------------------  
  420 IEXP = 2  
      IZ = IBETA  
  430 IF (K .LT. IZ) GO TO 440  
         IZ = IZ * IBETA  
         IEXP = IEXP + 1  
      GO TO 430  
  440 MX = IZ + IZ - 1  
!-----------------------------------------------------------------------  
!  Loop to determine MINEXP, XMIN.  
!  Exit from loop is signaled by an underflow.  
!-----------------------------------------------------------------------  
  450 XMIN = Y  
         Y = Y * BETAIN  
!-----------------------------------------------------------------------  
!  Check for underflow here.  
!-----------------------------------------------------------------------  
         A = Y * ONE  
         TEMP = Y * T  
         IF (((A+A) .EQ. ZERO) .OR. (ABS(Y) .GE. XMIN)) GO TO 460  
         K = K + 1  
         TEMP1 = TEMP * BETAIN  
         IF (TEMP1*BETA .NE. Y) GO TO 450  
      NXRES = 3  
      XMIN = Y  
  460 MINEXP = -K  
!-----------------------------------------------------------------------  
!  Determine MAXEXP, XMAX.  
!-----------------------------------------------------------------------  
      IF ((MX .GT. K+K-3) .OR. (IBETA .EQ. 10)) GO TO 500  
      MX = MX + MX  
      IEXP = IEXP + 1  
  500 MAXEXP = MX + MINEXP  
!-----------------------------------------------------------------  
!  Adjust IRND to reflect partial underflow.  
!-----------------------------------------------------------------  
      IRND = IRND + NXRES  
!-----------------------------------------------------------------  
!  Adjust for IEEE-style machines.  
!-----------------------------------------------------------------  
      IF ((IRND .EQ. 2) .OR. (IRND .EQ. 5)) MAXEXP = MAXEXP - 2  
!-----------------------------------------------------------------  
!  Adjust for non-IEEE machines with partial underflow.  
!-----------------------------------------------------------------  
      IF ((IRND .EQ. 3) .OR. (IRND .EQ. 4)) MAXEXP = MAXEXP - IT  
!-----------------------------------------------------------------  
!  Adjust for machines with implicit leading bit in binary  
!  significand, and machines with radix point at extreme  
!  right of significand.  
!-----------------------------------------------------------------  
      I = MAXEXP + MINEXP  
      IF ((IBETA .EQ. 2) .AND. (I .EQ. 0)) MAXEXP = MAXEXP - 1  
      IF (I .GT. 20) MAXEXP = MAXEXP - 1  
      IF (A .NE. Y) MAXEXP = MAXEXP - 2  
      XMAX = ONE - EPSNEG  
      IF (XMAX*ONE .NE. XMAX) XMAX = ONE - BETA * EPSNEG  
      XMAX = XMAX / (BETA * BETA * BETA * XMIN)  
      I = MAXEXP + MINEXP + 3  
      IF (I .LE. 0) GO TO 520  
      DO 510 J = 1, I  
          IF (IBETA .EQ. 2) XMAX = XMAX + XMAX  
          IF (IBETA .NE. 2) XMAX = XMAX * BETA  
  510 CONTINUE  
      RMACH(1) = XMIN  
      RMACH(2) = XMAX  
      RMACH(3) = EPSNEG  
      RMACH(4) = EPS  
  520 RETURN  
  
      END SUBROUTINE T665R  
  
      end module DISORD_SUBS  
