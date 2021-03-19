************************************************************************
*
*     -------------------------------
      SUBROUTINE KZMJMP(NEVJMP,ICOND)
*     -------------------------------
*
*     (Purpose)
*        Skip some sub-events and change the current event
*        in the event buffer
*
*     (Inputs)
*        NEVJMP ; no. of events to be skipped
*
*     (Output)
*        ICOND  ; error condition flag
*               ; = 0 normal return
*               ; = 1 if skipped events exceeds the buffer
*
*     (Creation Date)
*        1994.06.15 : First version by J.Kanzaki.
*
************************************************************************


#include "zbssys.h"
#include "zbs.h"


      ICOND = 0

      CALL KZMNUM(NSBEVT,ISBEVT)

      LEVC = LZBS(IPEVC)
      LEMCR = LZBS(LEVC-ILEMCR)
      LEVT = LZBS(LEMCR-ILEHCR)

      IF(NEVJMP.GT.0) THEN

         IF(ISBEVT+NEVJMP .GT. NSBEVT) THEN

            ICOND = 1

         ELSE

            DO 10 I=1,NEVJMP
            LEVT = LZBS(LEVT)
   10       CONTINUE
            LZBS(LEMCR-ILEHCR) = LEVT
            CALL KZMKDT

         END IF

      ELSE

         IF(ISBEVT+NEVJMP .LE. 0) THEN

            ICOND = 1

         ELSE

            DO 20 I=1,-NEVJMP
            LEVT = LZBS(LEVT+2)
   20       CONTINUE
            LZBS(LEMCR-ILEHCR) = LEVT
            CALL KZMKDT

         END IF

      END IF

      RETURN
      END
