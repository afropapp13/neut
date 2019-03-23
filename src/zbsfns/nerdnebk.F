***********************************************************************
*     ------------------------
      SUBROUTINE NERDNEBK(POS)
*     ------------------------
*     
*     ( purpose )
*       Read NEUT bank
*
*     ( input ) 
*       none
*
*     ( output ) 
*       POS : INTERACTION POSITION
*       COMMON /NEWORK/      
*
*     ( creation date and author )
*       1995.03.16 ; First version by K.Kaneyuki
*
***********************************************************************
#include "nework.h"

      REAL POS(3)
      CHARACTER*32 CNAME
      INTEGER IDATA(20)
      REAL    RDATA(20)
      EQUIVALENCE (IDATA(1),RDATA(1))

      CNAME='NEUT'
      CALL KZBLOC(CNAME,IPBANK)
      IF (IPBANK.EQ.0) GOTO 1000
C
C     -- READ EMB
C
C     LEM+1 : # OF PARTICLE
C     LEM+2 : MODE OF INTERACTION
C     LEM+3 : ISOL   ( FLAG OF SOLACT )
C     LEM+4 : SOLACT ( SOLAR ACTIVITY )
C
      ISEG=0
      CALL KZGET1(CNAME,ISEG,NDATA,IDATA)
      IF (NDATA.EQ.0) GOTO 2000
      NUMNE =IDATA(1)
      MODENE=IDATA(2)
      POS(1)=RDATA(5)
      POS(2)=RDATA(6)
      POS(3)=RDATA(7)

      DO 10 I=1,NUMNE
C     
C     -- READ VECTOR INFORMATION OF M.C. ON I-TH PARTICLE
C
C        IND+ 1 : PARTICLE CODE
C           + 2 : PX ( MEV )
C           + 3 : PY ( MEV )
C           + 4 : PZ ( MEV )
C           + 5 : PARENT PARTICLE
C           + 6 : FLAG OF FINAL STATE
C           + 7 : CHERENKOV FLAG
C
         ISEG=I
         CALL KZGET1(CNAME,ISEG,NDATA,IDATA)

         IPNE(I)  =IDATA( 1)
         PNE(1,I) =RDATA(2)
         PNE(2,I) =RDATA(3)
         PNE(3,I) =RDATA(4)
         IORGNE(I)=IDATA(5)
         IFLGNE(I)=IDATA(6)
         ICRNNE(I)=IDATA(7)
 10   CONTINUE

      RETURN
C     
C     ++ ERROR MESSAGE
C     
 1000 WRITE(6,900)
 900  FORMAT(' *** ERROR *** ( NEUT DOES NOT EXIST IN NERDNEBK )')
      RETURN
 2000 WRITE(6,910) ISEG
 910  FORMAT(' *** ERROR *** ( NEUT BANK : ',I2,
     $     'DOES NOT EXIST IN NERDNEBK )')
      RETURN
      END
