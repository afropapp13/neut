************************************************************************
*
*     ------------------------
      SUBROUTINE KZMNXT(ICOND)
*     ------------------------
*
*     (Purpose)
*        Select next sub-event in the event buffer for the
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
      LECR = LZBS(LEVC-ILEMCR)
      LEHLS = LZBS(LECR-ILEHLS)
      LLST = LZLAST(IXSTOR,LEHLS)

      LCUR = LZBS(LECR-ILEHCR)

      IF(LCUR.NE.LLST) THEN

         LNXT = LZBS(LCUR)
         LZBS(LECR-ILEHCR) = LNXT

      ELSE

         ICOND = 1

      END IF

      CALL KZMKDT

      RETURN
      END
