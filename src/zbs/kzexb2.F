************************************************************************
*
*     -----------------------------
      SUBROUTINE KZEXB2(LBIR,ISEGM)
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
      ISDIFF=ISEGM-NSEGM
      IF(ISDIFF.LE.0) RETURN
      CALL MZPUSH(IXSTOR,LBIR,0,LINDSG*ISDIFF,' ')
      IZBS(LBIR+INSEGM)=ISEGM
      DO 100 I=NSEGM+1,ISEGM
      LIND=LBIR+LINDHD+I*LINDSG
      IF(I.EQ.0) THEN
         IZBS(LIND+1)=0
         IZBS(LIND+2)=0
      ELSE
         IZBS(LIND+1)=-1
         IZBS(LIND+2)=0
      END IF
  100 CONTINUE

      RETURN
      END
