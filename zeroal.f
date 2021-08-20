      SUBROUTINE ZEROAL( EXPBEA, FLYR, OPRIM, TAUCPR, XR0, XR1,
     &                   CMU, CWT, PSI, WK, Z0, Z1, ZJ,
     &                   HLPR, YLM0,
     &                   ARRAY, CC, EVECC,
     &                   GL,
     &                   YLMC,
     &                   YLMU,
     &                   KK, LL, ZZ, ZPLK0, ZPLK1,
     &                   GC,
     &                   LAYRU, UTAUPR,
     &                   GU,
     &                   Z0U, Z1U, ZBEAM,
     &                   EVAL,
     &                   AMB, APB,
     &                   IPVT, Z,
     &                   RFLDIR, RFLDN, FLUP, UAVG, DFDT,
     &                   ALBMED, TRNMED,
     &                   U0U,
     &                   UU )

c         ZERO ARRAYS; NDn is dimension of all arrays following
c         it in the argument list

c   Called by- DISORT
c --------------------------------------------------------------------

c     .. Array Arguments ..

      INTEGER, intent(out) ::   IPVT( : ), LAYRU( : )
      REAL, intent(out)    ::
     $          ALBMED( : ), AMB( : ), APB( : ), ARRAY(:,:), CC(:,:),
     &          CMU( : ), CWT( : ), DFDT( : ), EVAL( : ), EVECC(:,:),
     &          EXPBEA( : ), FLUP( : ), FLYR( : ), GC( : ), GL( : ),
     &          GU( : ), HLPR( : ), KK( : ), LL( : ), OPRIM( : ),
     &          PSI( : ), RFLDIR( : ), RFLDN( : ), TAUCPR( : ),
     &          TRNMED( : ), U0U( : ), UAVG( : ), UTAUPR( : ), UU( : ),
     &          WK( : ), XR0( : ), XR1( : ), YLM0( : ), YLMC( : ),
     &          YLMU( : ), Z( : ), Z0( : ), Z0U( : ), Z1( : ), Z1U( : ),
     &          ZBEAM( : ), ZJ( : ), ZPLK0( : ), ZPLK1( : ), ZZ( : )

c     ..

      EXPBEA = 0.0
      FLYR   = 0.0
      OPRIM  = 0.0
      TAUCPR = 0.0
      XR0    = 0.0
      XR1    = 0.0

         CMU( N ) = 0.0
         CWT( N ) = 0.0
         PSI( N ) = 0.0
         WK( N )  = 0.0
         Z0( N )  = 0.0
         Z1( N )  = 0.0
         ZJ( N )  = 0.0

         HLPR( N ) = 0.0
         YLM0( N ) = 0.0

         ARRAY( N ) = 0.0
         CC( N )    = 0.0
         EVECC( N ) = 0.0

         GL( N ) = 0.0

         YLMC( N ) = 0.0

         YLMU( N ) = 0.0

         KK( N )    = 0.0
         LL( N )    = 0.0
         ZZ( N )    = 0.0
         ZPLK0( N ) = 0.0
         ZPLK1( N ) = 0.0

         GC( N ) = 0.0

         LAYRU( N )  = 0
         UTAUPR( N ) = 0.0

         GU( N ) = 0.0

         Z0U( N )   = 0.0
         Z1U( N )   = 0.0
         ZBEAM( N ) = 0.0

         EVAL( N ) = 0.0

         AMB( N ) = 0.0
         APB( N ) = 0.0

         IPVT( N ) = 0
         Z( N )    = 0.0

         RFLDIR( N ) = 0.
         RFLDN( N )  = 0.
         FLUP( N )   = 0.
         UAVG( N )   = 0.
         DFDT( N )   = 0.

         ALBMED( N ) = 0.
         TRNMED( N ) = 0.

         U0U( N ) = 0.

         UU( N ) = 0.

      END SUBROUTINE ZEROAL
