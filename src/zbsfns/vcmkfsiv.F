***********************************************************************
*     -----------------
      SUBROUTINE VCMKFSIV
*     -----------------
*
*     ( purpose )
*       Create FSIVTX bank
*
*     ( input ) 
*       none
*
*     ( output ) 
*       none
*
*     ( creation date and author )
*       2010.9.26 ; First version by Y.Hayato
*
***********************************************************************

      IMPLICIT NONE

#include "fsihist.h"

      INTEGER IDATA(20)
      REAL    RDATA(20)
      EQUIVALENCE(IDATA(1),RDATA(1))
      CHARACTER*32 CNAME

      INTEGER*4 NWMAST,NSEGM,LSEGM,I,ISEG,IERR

C
C     -- MAKE EMB
C
C     LEM+1 : # OF VERTEX
C
      IDATA(1)=NVERT

      CNAME='FSIVTX'

      NWMAST=1
      NSEGM=NVERT
      LSEGM=6

      CALL KZBCR1(CNAME,NWMAST,'1I',NSEGM,LSEGM,'3F 2I 1F ',IERR)      
      IF (IERR.NE.0) GOTO 1000
      CALL KZREP1(CNAME,0,NWMAST,IDATA)      

      DO 10 I=1,NVERT
C     
C     -- SET VERTEX INFORMATION 
C     
C        IND+ 1 : X
C           + 2 : Y
C           + 3 : Z
C           + 4 : KIND OF VERTEX
C          (+ 5 : N/A Paranet particle ??)
C          (+ 6 : dummy )

         RDATA( 1)=POSVERT(1,I)
         RDATA( 2)=POSVERT(2,I)
         RDATA( 3)=POSVERT(3,I)
         IDATA( 4)=IFLGVERT(I)
         IDATA( 5)=0
         RDATA( 6)=0.
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
  900 FORMAT(' *** ERROR *** ( BANK CREATION IN VCMKFSIV )')
      RETURN
      END
