************************************************************************
*
*     -----------------------------
      SUBROUTINE KZEXB1(LBIR,ISEGM)
*     -----------------------------
*
*     (Purpose)
*       Expand bank data records for fixed segment bank
*
*     (Input)
*       LBIR  ; Address of the bank index record
*       ISEGM ; Target segment# to be added
*
*     (Output)
*       none
*
*     (Creation Date and Author)
*       1994.07.07 ; First version by J.Kanzaki
*
************************************************************************


#include "zbs.h"
#include "zbssys.h"


      NSEGM =IZBS(LBIR+INSEGM)
      NEXTSG=IZBS(LBIR+INEXSG)
      NEXTSG=NEXTSG+1
      ISDIFF=ISEGM-NSEGM
      IF(ISDIFF.LE.0) RETURN
      NSEGM1=IZBS(LBIR+INSGM1)
      IF(ISDIFF.GT.NSEGM1) NSEGM1=(ISDIFF/NSEGM1+1)*NSEGM1
      IF(MOD(NEXTSG,2).EQ.0) NSEGM1=2*NSEGM1
      LSEGM = IZBS(LBIR+ILSEGF)
      LDTR=LZBS(LBIR-ILDATR)
      CALL MZPUSH(IXSTOR,LDTR,0,LSEGM*NSEGM1,' ')
C DWC Must restore local bank index LBIR after calling
C DWC MZPUSH
      LBIR=LZBS(LDTR+1)
      NSEGM=NSEGM+NSEGM1
      IZBS(LBIR+INSEGM)=NSEGM
      IZBS(LBIR+INEXSG)=NEXTSG
      IZBS(LBIR+INSGM1)=NSEGM1

      RETURN
      END
