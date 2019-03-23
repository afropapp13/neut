************************************************************************
*     -------------------------
      FUNCTION FNNCELSCRS(E,IPAR,ITARG)
*     -------------------------
*
*     (Purpose)
*       Cross section of neutral current elastic scattering
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
*       ITARG  : Target type
*                2212 : Proton
*                2112 : Neutron
*     
*     (Output)
*       FNNCELSCRS : CROSS SECTION ( 10^-38cm^2 )
*
*     (Creation Date and Author)
*       2011.03.12 ; Y.Hayato
*                    FNNCELSCRS ( modify FNELSCRS )
*     
************************************************************************
      IMPLICIT NONE

      REAL*4    FNNCELSCRS
      REAL*4    E
      INTEGER*4 IPAR,ITARG

      REAL*8    DGAUSS 

      EXTERNAL  DGAUSS
      REAL*8    Q2MINX, Q2MAXX, EPSIN
      REAL*8    QBUF

      REAL*8    ELSCRS

      REAL*8    EE
      INTEGER*4 IPR,ITRG
      COMMON /NEINT3/ EE,IPR,ITRG

      REAL*8    DNNCELA
      EXTERNAL  DNNCELA

      REAL*4    FNNCQ2MAX, FNNCQ2MIN
      EXTERNAL  FNNCQ2MAX, FNNCQ2MIN

C      write(*,*) 'fnncelscrs(E=,',e,',IPAR=',IPAR,')'
      EE=E
      IPR=IPAR
      ITRG=ITARG
      FNNCELSCRS=0.
      Q2MAXX=FNNCQ2MAX(E,ITARG)
      Q2MINX=FNNCQ2MIN(E,ITARG)
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

      ELSCRS = DGAUSS(DNNCELA,Q2MINX,Q2MAXX,EPSIN)
      
C      write(*,*) elscrs
      FNNCELSCRS = ELSCRS
      
      RETURN
      END
