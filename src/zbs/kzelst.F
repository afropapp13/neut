************************************************************************
*
*     -----------------
      SUBROUTINE KZELST
*     -----------------
*
*     (Purpose)
*        Set the last event as the current event
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

      LLST = LZLAST(IXSTOR,LZBS(LEVC-ILEMLS))

      LZBS(LEVC-ILEMCR) = LLST

      CALL KZMKDT

      RETURN
      END
