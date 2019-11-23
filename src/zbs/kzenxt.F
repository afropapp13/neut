************************************************************************
*
*     ------------------------
      SUBROUTINE KZENXT(ICOND)
*     ------------------------
*
*     (Purpose)
*        Select next event in the event buffer for the
*        current event
*
*     (Creation Date)
*        1994.06.15 : First version by J.Kanzaki.
*
************************************************************************

#include "zbs.h"
#include "zbssys.h"


      ICOND = 0

      LEVC = LZBS(IPEVC)

      LELS = LZBS(LEVC-ILEMLS)
      LLST = LZLAST(IXSTOR,LELS)

      LCUR = LZBS(LEVC-ILEMCR)

      IF(LCUR.NE.LLST) THEN

         LNXT = LZBS(LCUR)
         LZBS(LEVC-ILEMCR) = LNXT

      ELSE

         ICOND = 1

      END IF

      CALL KZMKDT

      RETURN
      END
