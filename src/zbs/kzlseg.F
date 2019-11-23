************************************************************************
*
*     ------------------------------------
      SUBROUTINE KZLSEG(CNAME,ISEGM,LSEGM)
*     ------------------------------------
*
*     (Purpose)
*       Get the length of the specified segment
*
*     (Input)
*       CNAME  ; The name of the bank
*       ISEGM  ; the segment number
*
*     (Output)
*       LSEGM  ; segment length
*
*     (Creation Date and Author)
*       1993.07.19 ; First version by J.Kanzaki
*
************************************************************************


#include "zbs.h"
#include "zbssys.h"

      CHARACTER*(*) CNAME

      LSEGM = 0

      CALL KZBLOC(CNAME, LBIR)

      IF(LBIR.LE.0) RETURN

      NSEGM=IZBS(LBIR+INSEGM)
      IF(NSEGM.LT.ISEGM) RETURN

      ITYPE=IZBS(LBIR+IBKTYP)

      IF(ITYPE.EQ.0) THEN

         LSEGM = IZBS(LBIR+LINDHD+ISEGM*LINDSG+2)

      ELSE

         IF(ISEGM.EQ.0) THEN

            LSEGM = IZBS(LBIR+INMAST)

         ELSE

            LSEGM = IZBS(LBIR+ILSEGF)

         END IF

      END IF

      RETURN
      END
