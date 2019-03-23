***********************************************************************
*     ------------------------
      SUBROUTINE NEMKCRSBK
*     ------------------------
*     
*     ( purpose )
*       Create NEUTCRS bank
*
*     ( input ) 
*       COMMON /NEUTCRS/
*
*     ( output ) 
*       NEUT BANK
*
*     ( creation date and author )
*       1995.03.16 ; First version by K.Kaneyuki
*       2011.11.13 ; Modified for NEUTCRS bank by P.de Perio
*
***********************************************************************
#include "neutcrs.h"

      INTEGER IDATA(500)
      REAL    RDATA(500)
      CHARACTER*32 CNAME
      EQUIVALENCE (RDATA(1),IDATA(1))
C
C     -- MAKE EMB FOR VERTEX BANK
C
C     LEM+1 : # OF VERTICES

      RDATA(1)=CRSX
      RDATA(2)=CRSY
      RDATA(3)=CRSZ
      RDATA(4)=CRSPHI

      CNAME='NEUTCRS'

      NWMAST=4
      NSEGM=1
      LSEGM=4

      CALL KZBCR1(CNAME,NWMAST,'4F',NSEGM,LSEGM,'4F',IERR)      
      IF (IERR.NE.0) GOTO 1000
      CALL KZREP1(CNAME,0,NWMAST,RDATA(1))      

      RETURN
C     
C     ++ ERROR MESSAGE
C     
 1000 WRITE(6,900)
  900 FORMAT(' *** ERROR *** ( BANK CREATION IN NEMKCRSBK )')
      RETURN
      END
