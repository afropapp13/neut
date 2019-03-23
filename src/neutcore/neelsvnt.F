************************************************************************
*     --------------------------------------
      SUBROUTINE NEELSVNT(IPAR,E,MODE,DNEUT)
*     --------------------------------------
*
*     (Purpose)
*       VECTOR GENERATION FOR N.C. ELASTIC EVENT ON FREE NUCLEONS
*
*     (Input)
*       IPAR   : NEUTRINO TYPE
*                  12 : NEU-E
*                 -12 : NEU-E-BAR            
*                  14 : NEU-MU
*                 -14 : NEU-MU-BAR            
*                  16 : NEU-TAU
*                 -16 : NEU-TAU-BAR
*       E      : NEUTRINO ENERGY ( GEV )
*       DNEUT(3) : DIRECTION OF NEUTRINO
*     
*     (Output)
*       COMMON /NEWORK/
*
*     (Creation Date and Author)
*       1988.10.08 ; T.KAJITA   
*          ***** AT PRESENT THE CROSS SECTIONS ARE NOT EXACT.
*          ***** CROSS SECTIONS ARE SIMPLY SCALED FROM THAT OF CC ELASTIC
*          ***** FOR THESE SCALING DATA FROM K.ABE ET AL., P.R.L.56(1986)1107
*          ARE USED
*          ***** SIGMA(NEU N)=SIGMA(NEU P) IS ASSUMED
*       1995.02.17 ; K. KANEYUKI  FOR S.K.
*                    ELSVNT -> FNELSVNT
*     
************************************************************************
      IMPLICIT NONE

#include "nework.h"
#include "neutmodel.h"

      INTEGER*4 IPAR
      REAL*4    E,DNEUT(3)
      INTEGER*4 MODE

      INTEGER IORG(4)/ 0, 0, 1, 2/
      INTEGER IFLG(4)/-1,-1, 2, 0/
      INTEGER ICRN(4)/ 0, 0, 1, 1/

      INTEGER*4 IPTEMP,I      
      
      IPTEMP=0
      IF(IPAR.EQ. 12 .OR. IPAR.EQ. 14 .OR. IPAR.EQ. 16) IPTEMP=12
      IF(IPAR.EQ.-12 .OR. IPAR.EQ.-14 .OR. IPAR.EQ.-16) IPTEMP=-12
      IF(IPTEMP.EQ.0) RETURN

      if ( mod(MDLQE,100)/10.eq.0 ) then
         CALL NEELSVCT(IPTEMP,E,DNEUT)

         MODENE=MODE
         IF (ABS(MODE).EQ.51) THEN
            IPNE(2)=2212
            IPNE(4)=2212
         ELSE
            IPNE(2)=2112
            IPNE(4)=2112
         ENDIF
         IPNE(1)=IPAR
         IPNE(3)=IPAR
         DO 10 I=1,4
            IORGNE(I)=IORG(I)
            IFLGNE(I)=IFLG(I)
            ICRNNE(I)=ICRN(I)
 10      CONTINUE
      else
         call NENCELSVNT(IPAR,E,MODE,DNEUT)
      endif

      RETURN
      END
