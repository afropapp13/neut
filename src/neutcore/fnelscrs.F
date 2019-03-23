************************************************************************
*     -------------------------
      FUNCTION FNELSCRS(E,IPAR)
*     -------------------------
*
*     (Purpose)
*       Cross section of elastic scattering
*
*     (Input)
*       E      : NEUTRINO ENERGY ( GEV )
*       IPAR   : NEUTRINO TYPE
*                  12 : NEU-E
*                 -12 : NEU-E-BAR            
*                  14 : NEU-MU
*                 -14 : NEU-MU-BAR            
*                  16 : NEU-TAU
*                 -16 : NEU-TAU-BAR
*     
*     (Output)
*       FNELSCRS : CROSS SECTION ( 10^-38cm^2 )
*
*     (Creation Date and Author)
*       19??.??.?? ; M.Nakahata
*       1987.06.26 ; N.Sato
*       1995.01.30 ; K.KANEYUKI  FOR S.K.
*                    ELSCRS -> FNELSCRS
*       1995.05.11 ; K.KANEYUKI  GPINDP -> DGAUSS ( INTEGRATION )
*     
************************************************************************
      IMPLICIT NONE

      REAL*4    FNELSCRS
      REAL*4    E
      INTEGER*4 IPAR

      REAL*8    DGAUSS 

      EXTERNAL  DGAUSS
      REAL*8    Q2MINX, Q2MAXX, EPSIN
      REAL*8    QBUF

      REAL*8    ELSCRS

      REAL*8    EE
      INTEGER*4 IPR
      COMMON /NEINT3/ EE,IPR

      REAL*8    DNELA
      EXTERNAL  DNELA

      REAL*4    FNQ2MAX, FNQ2MIN
      EXTERNAL  FNQ2MAX, FNQ2MIN

C      write(*,*) 'fnelscrs(E=,',e,',IPAR=',IPAR,')'
      EE=E
      IPR=IPAR
      FNELSCRS=0.
      IF(IABS(IPAR).EQ.14 .AND. E.LE.0.112)RETURN
      IF(IABS(IPAR).EQ.16 .AND. E.LE.3.48109)RETURN
      Q2MAXX=FNQ2MAX(E,IPAR)
      Q2MINX=FNQ2MIN(E,IPAR)
      QBUF=Q2MINX
      IF(Q2MINX.LT.0.) Q2MINX=0.
C
C****************************************
C     START DGAUSS (CERN lib)
C****************************************
C
C     DOUBLE PRECISION FUNCTION DGAUSS(FUNC,A,B,EPSIN)
C
C     PARAMETERS
C
C     FUNC    = FUNCTION ROUTINE FOR THE FUNCTION FUNC(X).TO BE DE-
C               CLARED EXTERNAL IN THE CALLING ROUTINE
C     A       = LOWER BOUNDARY
C     B       = UPPER BOUNDARY
C     EPSIN   = ACCURACY REQUIRED FOR THE APPROXINATION
C
      EPSIN = 1.0D-7

      ELSCRS = DGAUSS(DNELA,Q2MINX,Q2MAXX,EPSIN)
      
C      write(*,*) elscrs
      FNELSCRS = ELSCRS
      
      RETURN
      END
