***********************************************************************
*     -----------------
      SUBROUTINE VCRDFSIP
*     -----------------
*     
*     ( purpose )
*       Read FSIPART bank
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
      EQUIVALENCE (IDATA(1),RDATA(1))

      REAL*4    AMASS
      INTEGER*4 I,IPBANK,ISEG,NDATA

      CNAME='FSIPART'
      CALL KZBLOC(CNAME,IPBANK)
      IF (IPBANK.EQ.0) GOTO 1000
C
C     -- READ EMB
C
C     LEM+1 : # OF PARTICLE
C
      ISEG=0
      CALL KZGET1(CNAME,ISEG,NDATA,IDATA)
      IF (NDATA.EQ.0) GOTO 2000
      NVCVERT=IDATA(1)
      FSIPROB=RDATA(2)
C      WRITE(6,*) NVCVERT

      DO 10 I=1,NVCVERT
C     
C     -- READ VECTOR INFORMATION OF M.C. ON I-TH PARTICLE
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
C
         ISEG=I
         CALL KZGET1(CNAME,ISEG,NDATA,IDATA)

         IPVERT(I)   =IDATA( 1)
                               
         AMASS       =RDATA( 2)
                               
         ABSPVERT(I) =RDATA( 3)
                               
         DIRVERT(1,I)=RDATA( 4)
         DIRVERT(2,I)=RDATA( 5)
         DIRVERT(3,I)=RDATA( 6)
                               
         ABSTPVERT(I)=RDATA( 7)
                               
C         0.          =RDATA( 8)
C         0.          =RDATA( 9)
C         0.          =RDATA(10)
                               
         IVERTI(I)   =IDATA(11)
         IVERTF(I)   =IDATA(12)
                               
C         0           =IDATA(13)
                               
C         0.          =RDATA(14)
C         0.          =RDATA(15)
C         0.          =RDATA(16)

 10   CONTINUE

      RETURN
C     
C     ++ ERROR MESSAGE
C     
 1000 WRITE(6,900)
 900  FORMAT(' *** ERROR *** ( FSIPART DOES NOT EXIST IN VCRDFSIP )')
      RETURN
 2000 WRITE(6,910) ISEG
 910  FORMAT(' *** ERROR *** ( FSIPAR BANK : ',I2,
     $     'DOES NOT EXIST IN VCRDFSIP )')
      RETURN
      END
