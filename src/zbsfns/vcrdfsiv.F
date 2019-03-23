***********************************************************************
*     -----------------
      SUBROUTINE VCRDFSIV
*     -----------------
*
*     ( purpose )
*       Read FSIVTX bank
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

      CHARACTER*32 CNAME
      INTEGER IDATA(20)
      REAL    RDATA(20)
      EQUIVALENCE(IDATA(1),RDATA(1))

      INTEGER*4 IPBANK,ISEG,I,NDATA
C
C     BANK LOCATION
C
      CNAME='FSIVTX'
      CALL KZBLOC(CNAME,IPBANK)
      IF (IPBANK.EQ.0) GOTO 1000
C
C     -- READ EMB
C
C     LEM+1 : # OF VERTEX
C
      ISEG=0
      CALL KZGET1(CNAME,ISEG,NDATA,IDATA)
      NVERT=IDATA(1)
C      WRITE(6,*) NVERT

      DO 10 I=1,NVERT
C     
C     -- READ VERTEX INFORMATION 
C     
C        IND+ 1 : X
C           + 2 : Y
C           + 3 : Z
C           + 4 : KIND OF VERTEX
C          (+ 5 : N/A Paranet particle ??)
C          (+ 6 : dummy )

         ISEG=I
         CALL KZGET1(CNAME,ISEG,NDATA,IDATA)
         IF (NDATA.EQ.0) GOTO 2000

         POSVERT(1,I) =RDATA( 1)
         POSVERT(2,I) =RDATA( 2)
         POSVERT(3,I) =RDATA( 3)
         IFLGVERT(I)  =IDATA( 4)
C         IDATA( 5)
C         RDATA( 6)

 10   CONTINUE

      RETURN
C     
C     ++ ERROR MESSAGE
C     
 1000 WRITE(6,900)
 900  FORMAT(' *** ERROR *** ( FSIVTX DOES NOT EXIST IN VCRDFSIV )')
      RETURN
 2000 WRITE(6,910) ISEG
 910  FORMAT(' *** ERROR *** ( FSIVTX BANK : ',I2,
     $     'DOES NOT EXIST IN VCRDFSIV )')
      RETURN
      END
