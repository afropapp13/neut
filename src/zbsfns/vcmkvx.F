***********************************************************************
*     -----------------
      SUBROUTINE VCMKVX
*     -----------------
*
*     ( purpose )
*       Create MCVERTEX bank
*
*     ( input ) 
*       none
*
*     ( output ) 
*       none
*
*     ( creation date and author )
*       1992.6.5 ; First version by K.Kaneyuki
*       1994.9.19; mod. for zbs V1.5
*
***********************************************************************

#include "vcwork.h"
#include "vcvrtx.h"

      INTEGER IDATA(20)
      REAL    RDATA(20)
      EQUIVALENCE(IDATA(1),RDATA(1))
      CHARACTER*32 CNAME
C
C     -- MAKE EMB
C
C     LEM+1 : # OF VERTEX
C
      IDATA(1)=NVTXVC

      CNAME='MCVERTEX'

      NWMAST=1
      NSEGM=NVTXVC
      LSEGM=6

      CALL KZBCR1(CNAME,NWMAST,'1I',NSEGM,LSEGM,'3F 2I 1F ',IERR)      
      IF (IERR.NE.0) GOTO 1000
      CALL KZREP1(CNAME,0,NWMAST,IDATA)      

      DO 10 I=1,NVTXVC
C     
C     -- SET VERTEX INFORMATION 
C     
C        IND+ 1 : X
C           + 2 : Y
C           + 3 : Z
C           + 4 : KIND OF VERTEX
C           + 5 : PARENT PARTICLE
C           + 6 : TIME

         RDATA( 1)=PVTXVC(1,I)
         RDATA( 2)=PVTXVC(2,I)
         RDATA( 3)=PVTXVC(3,I)
         IDATA( 4)=IFLVVC(I)
         IDATA( 5)=IPARVC(I)
         RDATA( 6)=TIMVVC(I)
C     
C     -- MAKE I-TH SEGMENT
C     
         ISEG=I
         CALL KZREP1(CNAME,ISEG,6,IDATA(1))         
 10   CONTINUE

      RETURN
C     
C     ++ ERROR MESSAGE
C     
 1000 WRITE(6,900)
  900 FORMAT(' *** ERROR *** ( BANK CREATION IN VCMKVX )')
      RETURN
      END
