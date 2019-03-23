***********************************************************************
*     -----------------
      SUBROUTINE VCRDVX
*     -----------------
*
*     ( purpose )
*       Read MCVERTEX bank
*
*     ( input ) 
*       none
*
*     ( output ) 
*       none
*
*     ( creation date and author )
*       1992.6.30 ; First version by K.Kaneyuki
*
***********************************************************************
#include "vcwork.h"
#include "vcvrtx.h"

      CHARACTER*32 CNAME
      INTEGER IDATA(20)
      REAL    RDATA(20)
      EQUIVALENCE(IDATA(1),RDATA(1))
C
C     BANK LOCATION
C
      CNAME='MCVERTEX'
      CALL KZBLOC(CNAME,IPBANK)
      IF (IPBANK.EQ.0) GOTO 1000
C
C     -- READ EMB
C
C     LEM+1 : # OF VERTEX
C
      ISEG=0
      CALL KZGET1(CNAME,ISEG,NDATA,IDATA)
      NVTXVC=IDATA(1)
C      WRITE(6,*) NVTXVC

      DO 10 I=1,NVTXVC
C     
C     -- READ VERTEX INFORMATION 
C     
C        IND+ 1 : X
C           + 2 : Y
C           + 3 : Z
C           + 4 : KIND OF VERTEX
C           + 5 : PARENT PARTICLE
C           + 6 : TIME

         ISEG=I
         CALL KZGET1(CNAME,ISEG,NDATA,IDATA)
         IF (NDATA.EQ.0) GOTO 2000

         PVTXVC(1,I)=RDATA(1)
         PVTXVC(2,I)=RDATA(2)
         PVTXVC(3,I)=RDATA(3)
         IFLVVC(I)  =IDATA(4)
         IPARVC(I)  =IDATA(5)
         TIMVVC(I)  =RDATA(6)
         DO 20 J=1,NVC
            IF (IVTIVC(J).EQ.I) 
     $           CALL UCOPY(PVTXVC(1,I),POSIVC(1,J),3)
            IF (IVTFVC(J).EQ.I) 
     $           CALL UCOPY(PVTXVC(1,I),POSFVC(1,J),3)
 20      CONTINUE
 10   CONTINUE

      RETURN
C     
C     ++ ERROR MESSAGE
C     
 1000 WRITE(6,900)
 900  FORMAT(' *** ERROR *** ( MCVERTEX DOES NOT EXIST IN VCRDVX )')
      RETURN
 2000 WRITE(6,910) ISEG
 910  FORMAT(' *** ERROR *** ( MCVERTEX BANK : ',I2,
     $     'DOES NOT EXIST IN VCRDVX )')
      RETURN
      END
