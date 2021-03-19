************************************************************************
*
*     -----------------------------
      SUBROUTINE KZEXB0(LBIR,ISEGM)
*     -----------------------------
*
*     (Purpose)
*       Expand bank index record
*
*     (Input)
*       LBIR  ; Address of the bank index record
*       ISEGM ; Target segment# to be added
*
*     (Output)
*       none
*
*     (Creation Date and Author)
*       1991.12.25 ; First version by J.Kanzaki
*
************************************************************************


#include "zbs.h"
#include "zbssys.h"


      NSEGM =IZBS(LBIR+INSEGM)
      NEXTSG=IZBS(LBIR+INEXSG)
      NEXTSG=NEXTSG+1
      ISDIFF=ISEGM-NSEGM
      NSEGM1=IZBS(LBIR+INSGM1)
      IF(ISDIFF.GT.NSEGM1) NSEGM1=(ISDIFF/NSEGM1+1)*NSEGM1
      IF(MOD(NEXTSG,2).EQ.0) NSEGM1=2*NSEGM1
      CALL MZPUSH(IXSTOR,LBIR,0,LINDSG*NSEGM1,' ')
      NSEGM=NSEGM+NSEGM1
      IZBS(LBIR+INSEGM)=NSEGM
      IZBS(LBIR+INEXSG)=NEXTSG
      IZBS(LBIR+INSGM1)=NSEGM1

      RETURN
      END
