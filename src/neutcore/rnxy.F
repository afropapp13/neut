************************************************************************
*     --------------------------------------
      SUBROUTINE RNXY(IPAR,E,X,Y,ITYPE,ITARG,IERR)
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
*       ITARG  : TARGET NUCLEON
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
*       2016.02.28 ; C. Bronner
*                    No longer average over neutron and proton xsec, use xsec from the right target nucleon
*                    Fix definition of x and xmax when Bodek-Yang corrections are used.
*                    If Bodek-Yang is not used, keep the target mass corrections  
*       2016.05.16 ; C. Bronner
*                    Add a 10% safety factor on the max xsec found in the grid search
*       2016.07.20 ; C. Bronner
*                    Target mass corrections moved to PDF calculation (grv98_lo.F)
*
************************************************************************
      REAL *8 PMASS     
      REAL*8 DNDXDY
      EXTERNAL DNDXDY
      REAL*8 EE,XX,YY
      integer ierr
      real*8 xmp,xmn
      PARAMETER (xmp    = 0.93827231D0)
      PARAMETER (xmn    = 0.93956563D0)
C     
      IERR=0
C     -- Check target nucleon code makes sense
      IF(ITARG.NE.2112 .AND. ITARG.NE.2212)GO TO 1000

C CB  Assign target mass       
      IF (itarg.EQ.2112) THEN
         PMASS=xmn
      else
         PMASS=xmp
      endif
C     
C     -- FIRST FIND MAXIMUM VALUE
C     Do a grid scan in x and y to find the maximum of the double differential cross section in x and y
C     If this maximum is >0, will be used as upper bound to generate x and y using acceptance-rejection sampling 
      XYMAX=0.
      DO 100 I=1,21
         Y=FLOAT(I-1)*0.05
         IF (Y.LE.0)Y=0.00000001
         DO 110 J=1,21
C            X=FLOAT(J-1)*0.05
            X=FLOAT(J-1)*0.05+0.00000001
            XMAX=4.*E**2*(1-Y)/(2.*PMASS*Y*E)           
            IF(X.GE.XMAX)GO TO 110
            XMAX=1.        
            IF(X.GE.XMAX)GO TO 110
            
            XX = X
            YY = Y
            EE = E
            XYDIS=DNDXDY(IPAR,ITARG,itype,EE,XX,YY)
            IF(XYMAX.GT.XYDIS)GO TO 110
            XYMAX=XYDIS
  110    CONTINUE
  100 CONTINUE
      IF(XYMAX.EQ.0.)GO TO 1000

C     Adding a safety factor in case real maximum is missed in the grid search      
      XYMAX=XYMAX*1.1


C     -- GIVE RANDOM X AND Y
C     Generate x and y using acceptance-rejection sampling, based on the value of the double differential cross-section in x and y
  200 Y=RLU(IDUM)
      X=RLU(IDUM)
      XMAX=4.*E**2*(1-Y)/(2.*PMASS*Y*E)
      IF(X.GE.XMAX)GO TO 200
      XMAX=1
      IF(X.GE.XMAX)GO TO 200


      XYRAN=RLU(IDUM)*XYMAX
      XX = X
      YY = Y
      EE = E
      XY=DNDXDY(IPAR,ITARG,itype,EE,XX,YY) 

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
