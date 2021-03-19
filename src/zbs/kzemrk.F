************************************************************************
*
*     -----------------
      SUBROUTINE KZEMRK
*     -----------------
*
*     (Purpose)
*        Set event mark on the current event
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

      LCUR = LZBS(LEVC-ILEMCR)

      IZBS(LCUR+IEVMRK) = 1

      RETURN
      END
