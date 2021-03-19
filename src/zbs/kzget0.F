************************************************************************
*
*     -----------------------------------------
      SUBROUTINE KZGET0(CNAME,ISEGM,NDATA,IDATA)
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

      IF(ITYPE.NE.0) RETURN

      IPSEGM=IZBS(LBIR+LINDHD+ISEGM*LINDSG+1)
      LFILSG=IZBS(LBIR+LINDHD+ISEGM*LINDSG+2)
      IF(LFILSG.EQ.0) IPSEGM=-1

      LDAT=IPDATA+IPSEGM
CCCCC PRINT *,'KZGETP;IPSEGM = ',IPSEGM
      IF(IPSEGM.LT.0) RETURN

      LIOWD = 1
      NDATA0 = 0

CCCCC PRINT *,'KZGETP;LFILSG = ',LFILSG
CCCCC PRINT *,'KZGETP;LDAT = ',LDAT

  100 CONTINUE

      IOWORD = IZBS(LDAT+LIOWD)

      NDAT = IOWORD/16
      IT   = MOD(IOWORD,16)

CCCCC PRINT *,'KZGETP;IOWORD = ',IOWORD
CCCCC PRINT *,'KZGETP;NDAT,IT = ',NDAT,IT
      DO 110 IDAT=1,NDAT
      IDATA(NDATA0+IDAT) = IZBS(LDAT+LIOWD+IDAT)
  110 CONTINUE

      LIOWD = LIOWD+NDAT+1
      NDATA0 = NDATA0+NDAT

      IF(LIOWD.LE.LFILSG) GO TO 100

      NDATA = NDATA0

      RETURN
      END
