***********************************************************************
*     -----------------
      SUBROUTINE VCMKFSIP
*     -----------------
*     
*     ( purpose )
*       Create FSIPART bank
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

      INTEGER IDATA(100)
      REAL    RDATA(100)
      CHARACTER*32 CNAME
      EQUIVALENCE (IDATA(1),RDATA(1))

      INTEGER*4 I,NWMAST,NSEGM,LSEGM,ISEG,IERR
      
      REAL*4  AMASS
C
C     -- MAKE EMB
C
C     LEM+1 : # OF PARTICLE
C
      IDATA(1)=NVCVERT
      RDATA(2)=FSIPROB

      CNAME='FSIPART'

      NWMAST=2
      NSEGM=NVCVERT
      LSEGM=15

      CALL KZBCR1(CNAME,NWMAST,'1I 1F',NSEGM,LSEGM,'1I 9F 3I 3R',IERR)      
      IF (IERR.NE.0) GOTO 1000
      CALL KZREP1(CNAME,0,NWMAST,IDATA(1))      


CC    *** NOT SUPPORTED ***

      DO 10 I=1,NVCVERT

         CALL MCMASS(IPVERT(I),AMASS)

C     
C     -- SET VECTOR INFORMATION OF M.C. ON I-TH PARTICLE
C     
C        IND+ 1 : PARTICLE CODE
C           + 2 : Mass ( MeV / c^2 )
C           + 3 : Momentum in Lab frame ( MEV/C )
C           + 4 : X direction in Lab frame
C           + 5 : Y direction in Lab frame
C           + 6 : Z direction in Lab frame
C           + 7 : Momentum in nucleon rest frame ( MEV/C )
C           + 8 : X direction in nucleon rest
C           + 9 : Y direction in nucleon rest
C           +10 : Z direction in nucleon rest
C           +11 : VERTEX # ( Initial )
C           +12 : VERTEX # ( FINAL )

         IDATA( 1)=IPVERT(I)

         RDATA( 2)=AMASS

         RDATA( 3)=ABSPVERT(I)

         RDATA( 4)=DIRVERT(1,I)
         RDATA( 5)=DIRVERT(2,I)
         RDATA( 6)=DIRVERT(3,I)

         RDATA( 7)=ABSTPVERT(I)

         RDATA( 8)=0.
         RDATA( 9)=0.
         RDATA(10)=0.

         IDATA(11)=IVERTI(I)
         IDATA(12)=IVERTF(I)

         IDATA(13)=0

         RDATA(14)=0.
         RDATA(15)=0.
         RDATA(16)=0.

C     
C     -- MAKE I-TH SEGMENT
C     
         ISEG=I
         CALL KZREP1(CNAME,ISEG,8,IDATA)         

 10   CONTINUE

      RETURN
C     
C     ++ ERROR MESSAGE
C     
 1000 WRITE(6,900)
  900 FORMAT(' *** ERROR *** ( BANK CREATION IN VCMKVC )')
      RETURN
      END
