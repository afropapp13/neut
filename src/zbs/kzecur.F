************************************************************************
*
*     -------------------------------
      SUBROUTINE KZECUR(IEVENT,ICOND)
*     -------------------------------
*
*     (Purpose)
*        Set a events specified as the current event
*
*     (Inputs)
*        IEVENT ; event number of the events which should
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

      LCUR = LZBS(LEVC-ILEMCR)

      CALL KZENUM(NEVTOT,IEVCUR)

      LLAST = LZLAST(IXSTOR,LZBS(LEVC-ILEMLS))

      IF(NEVTOT.LE.0 .OR. IEVENT.GT.NEVTOT) THEN

         ICOND = 1
         RETURN

      END IF

      LEVT = LZBS(LEVC-ILEMLS)
      IEVE = 0

 1000 CONTINUE

      IEVE = IEVE+1

      IF(IEVE .EQ. IEVENT) THEN

         LZBS(LEVC-ILEMCR) = LEVT
         CALL KZMKDT
         RETURN

      END IF

      LEVT = LZBS(LEVT)
      IF(LEVT.GT.0) GO TO 1000

      ICOND = 2

      RETURN
      END
