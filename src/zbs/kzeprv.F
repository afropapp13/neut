************************************************************************
*
*     ------------------------
      SUBROUTINE KZEPRV(ICOND)
*     ------------------------
*
*     (Purpose)
*        Select previous event in the event buffer for the
*        current event
*
*     (Creation Date)
*        1994.06.14 : First version by J.Kanzaki.
*
************************************************************************

#include "zbs.h"
#include "zbssys.h"


      ICOND = 0

      LEVC = LZBS(IPEVC)

      LELS = LZBS(LEVC-ILEMLS)
      LCUR = LZBS(LEVC-ILEMCR)

      IF(LCUR.NE.LELS) THEN

         LPRV = LZBS(LCUR+2)
         LZBS(LEVC-ILEMCR) = LPRV

      ELSE

         ICOND = 1

      END IF

      CALL KZMKDT

      RETURN
      END
