************************************************************************
*
*     -----------------
      SUBROUTINE KZMLST
*     -----------------
*
*     (Purpose)
*        Set the last sub-event as the current event
*
*     (Inputs)
*
*     (Output)
*
*     (Creation Date)
*        1994.06.15 : First version by J.Kanzaki.
*
************************************************************************


#include "zbssys.h"
#include "zbs.h"

      LEVC = LZBS(IPEVC)
      LECR = LZBS(LEVC-ILEMCR)
      LLST = LZLAST(IXSTOR,LZBS(LECR-ILEHLS))

      LZBS(LECR-ILEHCR) = LLST

      CALL KZMKDT

      RETURN
      END
