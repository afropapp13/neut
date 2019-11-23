************************************************************************
*
*     -----------------
      SUBROUTINE KZMFST
*     -----------------
*
*     (Purpose)
*        Set the first sub-event as the current event
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
      LFST = LZBS(LECR-ILEHLS)

      LZBS(LECR-ILEHCR) = LFST

      CALL KZMKDT

      RETURN
      END
