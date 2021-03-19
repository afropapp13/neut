************************************************************************
*
*     ----------------------------------------
      SUBROUTINE KZGET(CNAME,ISEGM,NDATA,IDATA)
*     ----------------------------------------
*
*     (Purpose)
*       Get the contents of the bank with the specified segment# 
*
*     (Input)
*       CNAME ; The name of the bank
*       ISEGM ; Segment number
*
*     (Output)
*       NDATA ; Length ot the data
*       IDATA(i) ; Contents of the segment
*
*     (Creation Date and Author)
*       1994.07.27 ; First version by J.Kanzaki
*
************************************************************************


#include "zbs.h"
#include "zbssys.h"

      CHARACTER*(*) CNAME

      INTEGER*4 IDATA(*)


      NDATA = 0

      CALL KZBLOC(CNAME, LBIR)

      IF(LBIR.LE.0) RETURN

      IF(ISEGM.LT.0) RETURN

      IPDATA=LZBS(LBIR-ILDATR)
      IF(LBIR.LE.0 .OR. IPDATA.LE.0) RETURN

      NSEGM=IZBS(LBIR+INSEGM)
      IF(NSEGM.LT.ISEGM) RETURN

      ITYPE=IZBS(LBIR+IBKTYP)

      IF(ITYPE.EQ.0) THEN

         CALL KZGET0(CNAME,ISEGM,NDATA,IDATA)

      ELSE IF(ITYPE.EQ.1) THEN

         CALL KZGET1(CNAME,ISEGM,NDATA,IDATA)

      END IF

      RETURN
      END
