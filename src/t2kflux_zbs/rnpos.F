***********************************************************************     
*     --------------------------
      SUBROUTINE RNPOS(POS,WALL)
*     --------------------------
*     
*     ( purpose )
*       SET RANDOM POSITION
*
*     ( input ) 
*       WALL : SET RANDUM POSITION INSIDE WALL(cm) FROM PMT PLANE
*
*     ( output ) 
*       POS(3) ( cm )
*
*     ( creation date and author )
*       1981.5.28 ; First version by M.Nakahata
*       1988.7.21 ; Modified by M.Kajita
*                   add WALL to argument
*       1992.6.30 ; Modified by K.Kaneyuki
*                   Super Kamioka version
*       2001. 2. 9 ; RNDM -> RLU
*
***********************************************************************     
#include "geotnk.h"
      DIMENSION POS(3)
      INTEGER IDUM
C     
      RRR=RINTK-WALL
C      POS(3)=(HIINTK-WALL*2)*RNDM(DUM)+(ZMINTK+WALL)
      POS(3)=(HIINTK-WALL*2)*RLU(IDUM)+(ZMINTK+WALL)
C     
C  100 POS(1)=(2.*RNDM(DUM)-1)*RRR
C      POS(2)=(2.*RNDM(DUM)-1)*RRR
  100 POS(1)=(2.*RLU(IDUM)-1)*RRR
      POS(2)=(2.*RLU(IDUM)-1)*RRR
      SSS=SQRT(POS(1)**2+POS(2)**2)
      IF(SSS.GT.RRR)GO TO 100

      RETURN
      END
