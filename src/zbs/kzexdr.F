************************************************************************
*
*     -----------------------------
      SUBROUTINE KZEXDR(LBIR,LDATA)
*     -----------------------------
*
*     (Purpose)
*       Expand data record.
*
*     (Input)
*       LBIR  ; Address of the bank index record
*       LDATA ; Target data length to be added.
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


      LDTR=LZBS(LBIR-ILDATR)

      LDTREC=IZBS(LBIR+ILDTRC)
      LFILLD=IZBS(LBIR+ILFILL)

      NEXTDT=IZBS(LBIR+INEXDT)
      NEXTDT=NEXTDT+1

      LDDIFF=LFILLD+LDATA-LDTREC

      LDTRC1=IZBS(LBIR+ILDTR1)
      IF(LDDIFF.GT.LDTRC1) LDTRC1=(LDDIFF/LDTRC1+1)*LDTRC1
      LDTRC1=2*LDTRC1
CCC   IF(MOD(NEXTDT,2).EQ.0) LDTRC1=2*LDTRC1
      CALL MZPUSH(IXSTOR,LDTR,0,LDTRC1,' ')
C DWC Must restore local bank index LBIR after MZPUSH,
C DWC in case of garbage collection
      LBIR=LZBS(LDTR+1)
      IZBS(LBIR+ILDTRC)=LDTREC+LDTRC1
      IZBS(LBIR+INEXDT)=NEXTDT
      IZBS(LBIR+ILDTR1)=LDTRC1

      RETURN
      END
