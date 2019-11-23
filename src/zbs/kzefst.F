************************************************************************
*
*     -----------------
      SUBROUTINE KZEFST
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

      LFST = LZBS(LEVC-ILEMLS)

      LZBS(LEVC-ILEMCR) = LFST

      CALL KZMKDT

      RETURN
      END
