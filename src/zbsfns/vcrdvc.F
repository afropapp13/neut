***********************************************************************
*     -----------------
      SUBROUTINE VCRDVC
*     -----------------
*     
*     ( purpose )
*       Read MCVECT bank
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
      EQUIVALENCE (IDATA(1),RDATA(1))

      CNAME='MCVECT'
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
      NVC=IDATA(1)
C      WRITE(6,*) NVC

      DO 10 I=1,NVC
C     
C     -- READ VECTOR INFORMATION OF M.C. ON I-TH PARTICLE
C     
C        IND+ 1 : PARTICLE CODE
C           + 2 : ENERGY ( MEV )
C           + 3 : PX ( MEV )
C           + 4 : PY ( MEV )
C           + 5 : PZ ( MEV )
C           + 6 : PARENT PARTICLE
C           + 7 : VERTEX # ( INITIAL )
C           + 8 : CHERENKOV FLAG
C           + 9 : FINAL STATE FLAG
C           +10 : VERTEX # ( FINAL )
C
         ISEG=I
         CALL KZGET1(CNAME,ISEG,NDATA,IDATA)

         IPVC(I)  =IDATA( 1)
         CALL MCMASS(IPVC(I),AMASS)
         AMASVC(I)=AMASS
         PVC(1,I) =RDATA( 3)
         PVC(2,I) =RDATA( 4)
         PVC(3,I) =RDATA( 5)
         IORGVC(I)=IDATA( 6)
         IVTIVC(I)=IDATA( 7)
         ICRNVC(I)=IDATA( 8)
         IFLGVC(I)=IDATA( 9)
         IVTFVC(I)=IDATA(10)
 10   CONTINUE

      RETURN
C     
C     ++ ERROR MESSAGE
C     
 1000 WRITE(6,900)
 900  FORMAT(' *** ERROR *** ( MCVECT DOES NOT EXIST IN VCRDVC )')
      RETURN
 2000 WRITE(6,910) ISEG
 910  FORMAT(' *** ERROR *** ( MCVECT BANK : ',I2,
     $     'DOES NOT EXIST IN VCRDVC )')
      RETURN
      END
