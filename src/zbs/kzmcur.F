************************************************************************
*
*     -------------------------------
      SUBROUTINE KZMCUR(ISBEVT,ICOND)
*     -------------------------------
*
*     (Purpose)
*        Set a sub-events specified as the current event
*
*     (Inputs)
*        ISBEVT ; Sub-event number of the events which should
*               ; be made to the current event
*
*     (Output)
*        ICOND     ; return condition flag
*              = 0 ; normal return
*              = 1 ; Specified event does not exist
*
*     (Creation Date)
*        1994.06.15 : First version by J.Kanzaki.
*
************************************************************************


#include "zbssys.h"
#include "zbs.h"


      ICOND = 0

      LEVC = LZBS(IPEVC)
      LECR = LZBS(LEVC-ILEMCR)
      LCUR = LZBS(LECR-ILEHCR)

      CALL KZMNUM(NSBTOT,ISBCUR)

      LLAST = LZLAST(IXSTOR,LZBS(LECR-ILEHLS))

      IF(NSBTOT.LE.0 .OR. ISBEVT.GT.NSBTOT) THEN

         ICOND = 1
         RETURN

      END IF

      LEVT = LZBS(LECR-ILEHLS)
      IEVE = 0

 1000 CONTINUE

      IEVE = IEVE+1

      IF(IEVE .EQ. ISBEVT) THEN

         LZBS(LECR-ILEHCR) = LEVT
         CALL KZMKDT
         RETURN

      END IF

      LEVT = LZBS(LEVT)
      IF(LEVT.GT.0) GO TO 1000

      ICOND = 2

      RETURN
      END
