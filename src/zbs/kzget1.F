************************************************************************
*
*     -----------------------------------------
      SUBROUTINE KZGET1(CNAME,ISEGM,NDATA,IDATA)
*     -----------------------------------------
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
*       1991.12.25 ; First version by J.Kanzaki
*
************************************************************************


#include "zbs.h"
#include "zbssys.h"

      CHARACTER*(*) CNAME

      INTEGER*4 LBIR,IDATA(*)


      NDATA = 0

      CALL KZBLOC(CNAME, LBIR)

      IF(LBIR.LE.0) RETURN

      IF(ISEGM.LT.0) RETURN

      IPDATA=LZBS(LBIR-ILDATR)
      IF(LBIR.LE.0 .OR. IPDATA.LE.0) RETURN

      ITYPE=IZBS(LBIR+IBKTYP)
      IF(ITYPE.NE.1) RETURN

      NWMAST = IZBS(LBIR+INMAST)
      NSEGM  = IZBS(LBIR+INSEGM)
      LSEGM  = IZBS(LBIR+ILSEGF)

      IF(ISEGM.LT.0 .OR. ISEGM.GT.NSEGM) THEN
         CALL MSGOUT('KZGET1','Illegal Segment# specified',
     &                ISEGM)
         RETURN
      END IF

      IF(ISEGM.EQ.0) THEN
         NDATA = NWMAST
         LDATA = IPDATA
      ELSE
         NDATA = LSEGM
         LDATA = IPDATA + NWMAST + LSEGM * (ISEGM-1)
      END IF
      DO 100 IDAT=1,NDATA
      IDATA(IDAT) = IZBS(LDATA+IDAT)
  100 CONTINUE

      RETURN
      END
