************************************************************************
*
*     -----------------------------------------
      SUBROUTINE KZWPOS(LBIR,ISEGM,NWPOS, LWPOS)
*     -----------------------------------------
*
*     (Purpose)
*       Get the contents of the bank with the specified segment# 
*
*     (Input)
*       LBIR  ; Address of the bank index record
*       ISEGM ; Segment number
*       NWPOS ; The position ( no. of words ) where the data
*             ; should be inserted.
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


      LWPOS = -1

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
      IF(IPSEGM.LT.0) RETURN

      LIOWD = 1
      NDATA0 = 0

  100 CONTINUE

      IOWORD = IZBS(LDAT+LIOWD)

      NDAT = IOWORD/16
      IT   = MOD(IOWORD,16)

      IF(NWPOS.LE.NDATA0+NDAT) THEN

         LWPOS = LIOWD + NWPOS - NDATA0

      ELSE

         LIOWD = LIOWD+NDAT+1
         NDATA0 = NDATA0+NDAT

         IF(LIOWD.LE.LFILSG) GO TO 100

      END IF
  
      RETURN
      END
