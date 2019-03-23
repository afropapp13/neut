************************************************************************
*     --------------------------------------
      SUBROUTINE RNXY(IPAR,E,X,Y,ITYPE,IERR)
*     --------------------------------------
*
*     (Purpose)
*     ++ GIVE RAMDOM X AND Y VALUE
*     DISTRIBUTION IS FIT FOR GGM DATA
*
*     (Input)
*       IPAR   : NEUTRINO SPECIES
*       E      : NEUTRINO ENERGY ( GEV )
*       ITYPE  : INTERACTION TYPE
*     
*     (Output)
*       X      : Q**2/(2Mv) BJORKEN x
*       Y      : v/E        BJORKEN y
*       IERR   : ERROR CODE
*
*     (Creation Date and Author)
*       1983.11.21 ; M.NAKAHATA
*       1995.03.03 ; K. KANEYUKI  FOR S.K.
*                    RNXY -> RNXY
*     
************************************************************************
      DATA PMASS/0.938/
      REAL*8 DNDXDY
      EXTERNAL DNDXDY
      REAL*8 EE,XX,YY
      integer ierr
C     
      IERR=0
C     
C     -- FIRST FIND MAXIMUM VALUE
C     
      XYMAX=0.
      DO 100 I=1,21
         Y=FLOAT(I-1)*0.05
         DO 110 J=1,21
C            X=FLOAT(J-1)*0.05
            X=FLOAT(J-1)*0.05+0.00000001
            XMAX=4.*E**2*(1-Y)/(PMASS**2+2.*PMASS*Y*E)
            IF(X.GE.XMAX)GO TO 110
            XMAX=2.*Y*E/(2.*Y*E+PMASS)
            IF(X.GE.XMAX)GO TO 110
C            XYDIS=DNDXDY(IPAR,E,X,Y,ITYPE)
            XX = X
            YY = Y
            EE = E
            XYDIS=(DNDXDY(IPAR,2212,itype,EE,XX,YY) 
     &             + DNDXDY(ipar,2112,itype,ee,xx,yy))/2
            IF(XYMAX.GT.XYDIS)GO TO 110
            XYMAX=XYDIS
  110    CONTINUE
  100 CONTINUE
      IF(XYMAX.EQ.0.)GO TO 1000
C     
C     -- GIVE RANDOM X AND Y
C     
C  200 Y=RNDM(DUM)
C      X=RNDM(DUM)
  200 Y=RLU(IDUM)
      X=RLU(IDUM)
      XMAX=4.*E**2*(1-Y)/(PMASS**2+2.*PMASS*Y*E)
      IF(X.GE.XMAX)GO TO 200
      XMAX=2.*Y*E/(2.*Y*E+PMASS)
      IF(X.GE.XMAX)GO TO 200
C     
C      XYRAN=RNDM(DUM)*XYMAX
      XYRAN=RLU(IDUM)*XYMAX
C      XY=DNDXDY(IPAR,E,X,Y,ITYPE)
            XX = X
            YY = Y
            EE = E
      XY=(DNDXDY(IPAR,2212,itype,EE,XX,YY) 
     &     + DNDXDY(ipar,2112,itype,ee,xx,yy))/2

      IF(XYRAN.GT.XY)GO TO 200

C     
      RETURN
C     
C     -- ERROR ROUTINE
C     
 1000 WRITE(6,600)E
  600 FORMAT(' *** WARNING *** ERROR AT RNXY XMAX=0.E=',G15.7)
      IERR=1
      RETURN
      END
