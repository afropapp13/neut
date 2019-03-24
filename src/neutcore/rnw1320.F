************************************************************************
*     -----------------------------------------
      SUBROUTINE RNW1320(IPAR,E,X,Y,W,ITYPE,ITARG,IERR)
*     -----------------------------------------
*
*     (Purpose)
*        ++ GIVE RANDOM VALUE OF X AND Y
*        THIS SUBROUTINE GIVE W VALUE WHITCH IS 1.3-2.0 GEV
*
*     (Input)
*       IPAR   : NEUTRINO SPECIES
*       E      : NEUTRINO ENERGY ( GEV )
*       ITYPE  : INTERACTION TYPE
*       ITARG  : TARGET NUCLEON
*     
*     (Output)
*       X      : Q**2/(2Mv) BJORKEN x
*       Y      : v/E        BJORKEN y
*       W      : 
*
*     (Creation Date and Author)
*       1988.09.06 ; T.KAJITA (MOD. FROM RANW MADE BY M.NAKAHATA)
*       1995.03.03 ; K. KANEYUKI  FOR S.K.
*                    RANW13 -> RNW13
*       2016.02.28 ; C. Bronner 
*                    Separate between interactions on protons and neutrons 
*                    Fix relation between x and Q2 when Bodek-Yang corrections are used.
*                    If Bodek-Yang is not used, keep the target mass corrections  
*       2016.07.20 ; C. Bronner
*                    Target mass corrections moved to PDF calculation (grv98_lo.F)    
*     
************************************************************************

      DATA MAXTRY/10000/
      REAL *8 PMASS 
      real*8 xmp,xmn
      PARAMETER (xmp    = 0.93827231D0)
      PARAMETER (xmn    = 0.93956563D0)
C     
      IERR=0
      ITRY=0

C     -- Check target nucleon code makes sense
      IF(ITARG.NE.2112 .AND. ITARG.NE.2212)GO TO 1000

C     Assign target mass       
      IF (itarg.EQ.2112) THEN
         PMASS=xmn
      else
         PMASS=xmp
      endif

C     
C     Generate x and y using rnxy, and check it gives a W within desired range
      DO 100 ITRY=1,MAXTRY
         CALL RNXY(IPAR,E,X,Y,ITYPE,ITARG,IERRR)
         IF(IERRR.EQ.1)GO TO 1001
         Q2=X*(2.*PMASS*Y*E)         
         XNEU=Y*E
         W2=-Q2+PMASS**2+2.*PMASS*XNEU
         W=SQRT(W2)      
         IF ( W.GT.1.3 .AND. W.LT.2.0 ) RETURN
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
