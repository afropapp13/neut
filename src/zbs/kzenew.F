************************************************************************
*
*     -----------------
      SUBROUTINE KZENEW
*     -----------------
*
*     (Purpose)
*        Create new event without clearing the event buffer
*
*     (Creation Date)
*        1994.07.06 : First version by J.Kanzaki.
*
************************************************************************

#include "zbs.h"
#include "zbssys.h"


      LEVC = LZBS(IPEVC)

      LZBS(LEVC-ILEMCR) = 0

      CALL KZCLDT

      RETURN
      END
