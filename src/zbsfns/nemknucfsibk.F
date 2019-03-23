***********************************************************************
*     ------------------------
      SUBROUTINE NEMKNUCFSIBK
*     ------------------------
*     
*     ( purpose )
*       Create NEUT NUCFSIHIST bank
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

      IDATA(1)=NFNVERT

      CNAME='NUCFSIHIST'

      NWMAST=1
      NSEGM=NFNVERT
      LSEGM=4

      CALL KZBCR1(CNAME,NWMAST,'1I',NSEGM,LSEGM,'2I 7F',IERR)      
      IF (IERR.NE.0) GOTO 1000
      CALL KZREP1(CNAME,0,NWMAST,IDATA(1))      

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
         IDATA( 1)=NFIFLAG(I)
         IDATA( 2)=NFfirststep(i)
         RDATA( 3)=NFx(I)
         RDATA( 4)=NFy(I)
         RDATA( 5)=NFz(I)
         RDATA( 6)=NFpx(I)
         RDATA( 7)=NFpy(I)
         RDATA( 8)=NFpz(I)
         RDATA( 9)=NFe(I)

C     
C     -- MAKE I-TH SEGMENT
C     
         ISEG=I
         CALL KZREP1(CNAME,ISEG,LSEGM,RDATA)         

 10   CONTINUE


C
C     -- MAKE EMB FOR INTERMEDIATE PARTICLE BANK
C
C     LEM+1 : # OF STEPs

      IDATA(1)=NFNSTEP

      CNAME='NUCFSIHISTSTEP'

      NWMAST=1
      NSEGM=NFNSTEP
      LSEGM=2

      CALL KZBCR1(CNAME,NWMAST,'1I',NSEGM,LSEGM,'2F',IERR)      
      IF (IERR.NE.0) GOTO 1000
      CALL KZREP1(CNAME,0,NWMAST,IDATA(1))      

      DO 20 I=1,NFNSTEP
C     
C     -- SET VECTOR INFORMATION OF M.C. ON I-TH INTERMEDIATE PARTICLE
C     
C        IND+ 1 : CMS energy at k-th step
C           + 2 : total probability at kth step

         RDATA( 1)=NFecms2(i)
         RDATA( 2)=NFptot(i)
C     
C     -- MAKE I-TH SEGMENT
C     
         ISEG=I
         CALL KZREP1(CNAME,ISEG,LSEGM,RDATA)         

 20   CONTINUE


      RETURN
C     
C     ++ ERROR MESSAGE
C     
 1000 WRITE(6,900)
  900 FORMAT(' *** ERROR *** ( BANK CREATION IN NEMKNUCFSIBK )')
      RETURN
      END
