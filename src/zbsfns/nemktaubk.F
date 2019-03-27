***********************************************************************
*     ------------------------
      SUBROUTINE NEMKTAUBK
*     ------------------------
*     
*     ( purpose )
*       Create TAU MC bank
*
*     ( input ) 
*       COMMON /TAUMC/
*
*     ( output ) 
*       TAU MC BANK
*
*     ( creation date and author )
*       2004.04.27 ; First version by C. Walter
*
***********************************************************************
#include "taumc.h"      

      INTEGER IDATA(100)
      REAL    RDATA(100)
      CHARACTER*32 CNAME
      EQUIVALENCE (IDATA(1),RDATA(1))
C
C     -- MAKE EMB
C
C     LEM+1 : MOMENTUM OF TAU (GEV)
C     LEM+2 : DIRECTION OF TAU (1)
C     LEM+3 : DIRECTION OF TAU (2)
C     LEM+4 : DIRECTION OF TAU (3)
C     LEM+5 : DECAY MODE (1-22)
C     LEM+6 : TAU POLARIZATION(1)
C     LEM+7 : TAU POLARIZATION(2)
C     LEM+8 : TAU POLARIZATION(3)
C
      RDATA(1)=TAUMOM
      RDATA(2)=TAUDIR(1)
      RDATA(3)=TAUDIR(2)
      RDATA(4)=TAUDIR(3)
      IDATA(5)=TAUDCY
      RDATA(6)=POLTAU(1)
      RDATA(7)=POLTAU(2)
      RDATA(8)=POLTAU(3)

      CNAME='TAUMC'

      LMAST=8
      NSEGM=0
      LSEGM=0

      CALL KZBCR1(CNAME,LMAST,'4F 1I 3F',NSEGM,LSEGM,'',IERR)      
      IF (IERR.NE.0) GOTO 1000
      CALL KZREP1(CNAME,0,LMAST,IDATA(1))      

      RETURN

C     
C     ++ ERROR MESSAGE
C     
 1000 WRITE(6,900)
  900 FORMAT(' *** ERROR *** ( BANK CREATION IN NEMKTAUBK )')
      RETURN
      END
