************************************************************************
*
*     -----------------------------------------------
      SUBROUTINE KZGET2(CNAME,ISEGM,ISSEG,NDATA,IDATA)
*     -----------------------------------------------
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
*       1992.12.22 ; First version by J.Kanzaki
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

      ITYPE=IZBS(LBIR+IBKTYP)
      IF(ITYPE.NE.2) RETURN

      NSEGM=IZBS(LBIR+INSEGM)
      IF(NSEGM.LT.ISEGM) RETURN

      NWMAST = IZBS(LBIR+INMAST)
      LSEGM  = IZBS(LBIR+ILSEGF)

      IF(ISEGM.LT.0 .OR. ISEGM.GT.NSEGM) THEN
         CALL MSGOUT('KZGETP','Illegal Segment# specified',
     &                ISEGM)
         RETURN
      END IF

      IF(ISEGM.EQ.0) THEN
         NDATA = NWMAST
         LDATA = IPDATA
      ELSE
         IPSEG=IZBS(LBIR+LINDHD+ISEGM*LINDSG+1)
         NSSEG=IZBS(LBIR+LINDHD+ISEGM*LINDSG+2)
         IF(ISSEG.LE.0 .OR. ISSEG.GT.NSSEG) RETURN
         NDATA = LSEGM
         LDATA = IPDATA + IPSEG + LSEGM * (ISSEG-1)
      END IF

      DO 100 IDAT=1,NDATA
      IDATA(IDAT) = IZBS(LDATA+IDAT)
  100 CONTINUE

      RETURN
      END
