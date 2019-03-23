***********************************************************************
*     -----------------
      SUBROUTINE VCMKDETGM
*     -----------------
*
*     ( purpose )
*       Create FXDETGEO bank
*
*     ( input ) 
*       none
*
*     ( output ) 
*       none
*
*     ( creation date and author )
*       2003.5.30 ; Y.Hayato
*
***********************************************************************

#include "fxdetgeo.h"

      INTEGER IDATA(20)
      REAL    RDATA(20)
      EQUIVALENCE(IDATA(1),RDATA(1))
      CHARACTER*32 CNAME
C
C     -- MAKE EMB
C
      CNAME='FXDETGEO'

C-- Version 
      IDATA(1)=1
      
      NWMAST=1
      NSEGM=1
      LSEGM=9

      CALL KZBCR1(CNAME,NWMAST,'1I',NSEGM,LSEGM,'1I 8F ',IERR)      
      IF (IERR.NE.0) GOTO 1000
      CALL KZREP1(CNAME,0,NWMAST,IDATA)      

C     
C     -- SET DETECTOR INFORMATION 
C     
C        IND+ 1 : Detector ID
C           + 2 : CENTER of the BEAM(X)
C           + 3 : CENTER of the BEAM(Y)
C           + 4 : CENTER of the DETECTOR(X)
C           + 5 : CENTER of the DETECTOR(Y)
C           + 6 : CENTER of the DETECTOR(Z)
C           + 7 : Size(half) of the DETECTOR(X)
C           + 8 : Size(half) of the DETECTOR(Y)

      IDATA( 1)=IDFXDET
      RDATA( 2)=BMCNT(1)
      RDATA( 3)=BMCNT(2)
      RDATA( 4)=FXDETCNT(1)
      RDATA( 5)=FXDETCNT(2)
      RDATA( 6)=FXDETCNT(3)
      RDATA( 7)=FXDETSZ(1)
      RDATA( 8)=FXDETSZ(2)
      
      ISEG=1
      CALL KZREP1(CNAME,ISEG,LSEGM,IDATA(1))         

      RETURN
C     
C     ++ ERROR MESSAGE
C     
 1000 WRITE(6,900)
  900 FORMAT(' *** ERROR *** ( BANK CREATION IN FXDETGEO )')
      RETURN
      END
