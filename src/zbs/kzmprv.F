************************************************************************
*
*     ------------------------
      SUBROUTINE KZMPRV(ICOND)
*     ------------------------
*
*     (Purpose)
*        Select previous sub-event in the event buffer for the
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
      LELS = LZBS(LECR-ILEHLS)
      LCUR = LZBS(LECR-ILEHCR)

      IF(LCUR.NE.LELS) THEN

         LPRV = LZBS(LCUR+2)
         LZBS(LECR-ILEHCR) = LPRV

      ELSE

         ICOND = 1

      END IF

      CALL KZMKDT

      RETURN
      END
