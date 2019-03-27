***********************************************************************
*     -----------------
      SUBROUTINE VCMKVCBF
*     -----------------
*
*     ( purpose )
*       Create MCVECTBF bank
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
*       199x.?.??; new bank MCVECTBF to store information before
*                  nucleon rescattering
*
***********************************************************************

      implicit none
#include "vcwork.h"
#include "vcvrtx.h"

      INTEGER IDATA(20)
      REAL    RDATA(20)
      CHARACTER*32 CNAME
      EQUIVALENCE (IDATA(1),RDATA(1))

      integer i,iseg
      integer nwmast,nsegm,lsegm,ierr

C
C     -- MAKE EMB
C
C     LEM+1 : # OF PARTICLE
C
      IDATA(1)=NVC

      CNAME='MCVECTBF'

      NWMAST=1
      NSEGM=NVC
      LSEGM=10

      CALL KZBCR1(CNAME,NWMAST,'1I',NSEGM,LSEGM,'1I 4F 5I',IERR)
      IF (IERR.NE.0) GOTO 1000
      CALL KZREP1(CNAME,0,NWMAST,IDATA(1))


CC    *** NOT SUPPORTED ***

      DO 10 I=1,NVC
C
C     -- SET VECTOR INFORMATION OF M.C. ON I-TH PARTICLE
C
C        IND+ 1 : PARTICLE CODE
C           + 2 : ENERGY ( MEV )
C           + 3 : PX ( MEV/C )
C           + 4 : PY ( MEV/C )
C           + 5 : PZ ( MEV/C )
C           + 6 : PARENT PARTICLE
C           + 7 : VERTEX # ( INITIAL )
C           + 8 : CHERENKOV FLAG
C           + 9 : FINAL STATE FLAG
C           +10 : VERTEX # ( FINAL )

         IDATA( 1)=IPVC(I)
         RDATA( 2)=SQRT(AMASVC(I)**2+PVC(1,I)**2+PVC(2,I)**2
     $                  +PVC(3,I)**2)
         RDATA( 3)=PVC(1,I)
         RDATA( 4)=PVC(2,I)
         RDATA( 5)=PVC(3,I)
         IDATA( 6)=IORGVC(I)
         IDATA( 7)=IVTIVC(I)
         IDATA( 8)=ICRNVC(I)
         IDATA( 9)=IFLGVC(I)
         IDATA(10)=IVTFVC(I)
C
C     -- MAKE I-TH SEGMENT
C
         ISEG=I
         CALL KZREP1(CNAME,ISEG,10,IDATA)

 10   CONTINUE

      RETURN
C
C     ++ ERROR MESSAGE
C
 1000 WRITE(6,900)
  900 FORMAT(' *** ERROR *** ( BANK CREATION IN VCMKVCBK )')
      RETURN
      END
