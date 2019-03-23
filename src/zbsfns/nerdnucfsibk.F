***********************************************************************
*     ------------------------
      SUBROUTINE NERDNUCFSIBK
*     ------------------------
*     
*     ( purpose )
*       Read NEUT NUCFSIHIST bank
*
*     ( input ) 
*       COMMON /NUCLEONFSIHIST/
*
*     ( output ) 
*       NEUT BANK
*
*     ( creation date and author )
*       1995.03.16 ; First version by K.Kaneyuki
*       2010.10.24 ; Modified for FSIHIST bank by P.de Perio
*       2014.11.24 ; Modified for NUCFSIHIST bank by Y. Hayato
*
***********************************************************************
#include "nucleonfsihist.h"

      INTEGER IDATA(500)
      REAL    RDATA(500)
      CHARACTER*32 CNAME
      EQUIVALENCE (RDATA(1),IDATA(1))
C
C     -- MAKE EMB FOR VERTEX BANK
C
C     LEM+1 : # OF VERTICES


      CNAME='NUCFSIHIST'
      CALL KZBLOC(CNAME,IPBANK)
      IF (IPBANK.EQ.0) GOTO 1000

      NWMAST=1

      CALL KZGET1(CNAME,NWMAST,NDATA,IDATA)
      IF (NDATA.EQ.0) GOTO 1000
      NFNVERT=IDATA(1)

      DO 10 I=1,NFNVERT
C     
C     -- SET VECTOR INFORMATION OF M.C. ON I-TH VERTEX
C     
C        IND+ 1 : INTERACTION FLAG OF VERTEX (see nucleonfsihist.h)
C           + 2 : first step index of this track 
C           + 3 : POS X (fm)
C           + 4 : POS Y (fm)
C           + 5 : POS Z (fm)
C           + 6 : Mom. X ()
C           + 7 : Mom. Y ()
C           + 8 : Mom. Z ()
C           + 9 : Energy ()
C
         ISEG=I
         CALL KZGET1(CNAME,ISEG,NDATA,RDATA)         
         if (NDATA.eq.0) goto 1000

         NFIFLAG(I)    =IDATA( 1)
         NFfirststep(i)=IDATA( 2)
         NFx(I)        =RDATA( 3)
         NFy(I)        =RDATA( 4)
         NFz(I)        =RDATA( 5)
         NFpx(I)       =RDATA( 6)
         NFpy(I)       =RDATA( 7)
         NFpz(I)       =RDATA( 8)
         NFe(I)        =RDATA( 9)

C     
C     -- MAKE I-TH SEGMENT
C     

 10   CONTINUE

C
C     -- MAKE EMB FOR INTERMEDIATE PARTICLE BANK
C
C     LEM+1 : # OF STEPs

      CNAME='NUCFSIHISTSTEP'

      CALL KZBLOC(CNAME,IPBANK)
      IF (IPBANK.EQ.0) GOTO 1000
      ISEG = 0

      CALL KZGET1(CNAME,ISEG,NDATA,IDATA(1))      
      IF (NDATA.EQ.0) GOTO 1000
      NFNSTEP = IDATA(1)

      DO 20 I=1,NFNSTEP
C     
C     -- SET VECTOR INFORMATION OF M.C. ON I-TH INTERMEDIATE PARTICLE
C     
C        IND+ 1 : CMS energy at k-th step
C           + 2 : total probability at kth step

         ISEG=I
         CALL KZGET1(CNAME,ISEG,NDATA,RDATA)         

         NFecms2(i)=RDATA( 1)
         NFptot(i) =RDATA( 2)

 20   CONTINUE

      RETURN
C     
C     ++ ERROR MESSAGE
C     
 1000 WRITE(6,900)
  900 FORMAT(' *** ERROR *** ( BANK CREATION IN NEMKNUCFSIBK )')
      RETURN
      END
