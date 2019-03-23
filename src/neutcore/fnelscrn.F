************************************************************************
*     ------------------------------
      FUNCTION FNELSCRN(E,IPAR,INUC)
*     ------------------------------
*
*     (Purpose)
*       Cross section of neutral current elastic Int. on Free Nucleons
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
*       INUC   : =2212 PROTON  =2112 NEUTRON
*     
*     (Output)
*       FNELSCRN : CROSS SECTION ( 10**-38cm^2 )
*
*     (Creation Date and Author)
*       1988.10.08 ; T.KAJITA   
*       1989.06.10 ; T.KAJITA  SMALL MODIFICATION  
*         ***** AT PRESENT THE CROSS SECTIONS ARE NOT EXACT.
*         ***** CROSS SECTIONS ARE SIMPLY SCALED FROM THAT OF CC ELASTIC
*         ***** FOR THESE SCALING DATA FROM K.ABE ET AL., P.R.L.56(1986)
*               1107 ARE USED
*         ***** SIGMA(NEU    N)=1.5*SIGMA(NEU    P)           (10-JUN-89)
*               SIGMA(NEUBAR N)=1.0*SIGMA(NEUBAR P)
*               ARE USED FROM C.H.ALBRIGHT ET AL.,P.R.D14(1976)1780
*       1995.01.30 ; K. KANEYUKI  FOR S.K.
*                    ELSCRN -> FNELSCRN
*     
************************************************************************
      IMPLICIT NONE

#include "neutparams.h"      
#include "neutmodel.h"

      REAL*4 FNELSCRN
      REAL*4 E
      INTEGER*4 IPAR,INUC

      REAL*4   FNELSCRS,FNNCELSCRS
      EXTERNAL FNELSCRS,FNNCELSCRS

C ++ CROSS SECTION RATIO OF NC AND CC (NC/CC) FOR NEUTRINO AND ANTI-NEU

      REAL*4 RNEUNC,RNEBNC
      DATA RNEUNC/0.153/
      DATA RNEBNC/0.218/
C     
      INTEGER*4 IPTEMP
      REAL*4    CRS
      

      FNELSCRN=0.

      IF(INUC.NE.2212 .AND. INUC.NE.2112) THEN
         WRITE(6,*) ' FNELSCRN : INPUT NUCLEON ERROR : ',INUC
         RETURN
      ENDIF

      IPTEMP=0
      IF(IPAR.EQ.12 .OR. IPAR.EQ.14 .OR. IPAR. EQ. 16) IPTEMP=12
      IF(IPAR.EQ.-12 .OR. IPAR.EQ.-14 .OR. IPAR. EQ. -16) IPTEMP=-12
      IF(ABS(IPTEMP).NE.12) THEN
         WRITE(6,*) ' FNELSCRN : INPUT PARTICLE ERROR : ',IPAR
         RETURN
      ENDIF

      if (mod(MDLQE,100)/10.eq.0) then
         CRS=FNELSCRS(E,IPTEMP)
         
         IF(IPAR.GT.0)  THEN
            FNELSCRN=CRS*RNEUNC
            IF(INUC.EQ.2112) FNELSCRN=FNELSCRN*1.5
         END IF
         
         IF(IPAR.LT.0) FNELSCRN=CRS*RNEBNC
      else
         FNELSCRN = FNNCELSCRS(E,IPAR,INUC)
      endif
      
      RETURN
      END
