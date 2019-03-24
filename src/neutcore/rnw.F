************************************************************************
*     ---------------------------------------
      SUBROUTINE RNW(IPAR,E,X,Y,W,ITYPE,IERR)
*     ---------------------------------------
*
*     (Purpose)
*        ++ GIVE RANDOM VALUE OF X AND Y
*        THIS SUBROUTINE GIVE W VALUE WHITCH IS GREATER THAN 1.4 GEV
*
*     (Input)
*       IPAR   : NEUTRINO SPECIES
*       E      : NEUTRINO ENERGY ( GEV )
*       ITYPE  : INTERACTION TYPE
*     
*     (Output)
*       X      : Q**2/(2Mv) BJORKEN x
*       Y      : v/E        BJORKEN y
*       W      : 
*
*     (Creation Date and Author)
*       1983.11.28 ; M.NAKAHATA
*       1995.03.03 ; K. KANEYUKI  FOR S.K.
*                    RANW -> RNW
*       2016.02.28 ; C. Bronner 
*                    Separate between interactions on protons and neutrons
*                    Fix relation between x and Q2 when Bodek-Yang corrections are used.
*                    If Bodek-Yang is not used, keep the target mass corrections 
*       2016.07.20 ; C. Bronner
*                    Target mass corrections moved to PDF calculation (grv98_lo.F)
*     
************************************************************************
      DATA PMASS/0.938/
      DATA MAXTRY/5000/

      IERR=0
      ITRY=0

      DO 100 ITRY=1,MAXTRY
         CALL RNXY(IPAR,E,X,Y,ITYPE,ITARG,IERRR)
         IF(IERRR.EQ.1)GO TO 1001
         Q2=X*(2.*PMASS*Y*E)
         XNEU=Y*E
         W2=-Q2+PMASS**2+2.*PMASS*XNEU
         W=SQRT(W2)
         IF(W.GT.1.4) RETURN
 100  CONTINUE
C     
C     -- ERROR ROUTINE
C     
 1000 WRITE(6,600)E,ITYPE
  600 FORMAT(' *** WARNING *** ERROR AT RANW TRY OVER. (E,ITYPE)=',
     &     G15.7,I6)
 1001 IERR=1
      RETURN
      END
