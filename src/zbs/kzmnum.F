************************************************************************
*
*     --------------------------------
      SUBROUTINE KZMNUM(NSBEVT,ISBEVT)
*     --------------------------------
*
*     (Purpose)
*        Returns the total number of sub-events and
*        the sub-event number of the current event
*        in the event buffer
*
*     (Inputs)
*
*     (Output)
*        NSBEVT ; Total number of sub-events
*        ISBEVT ; Sub-event number 
*
*     (Creation Date)
*        1994.06.15 : First version by J.Kanzaki.
*
************************************************************************


#include "zbssys.h"
#include "zbs.h"


***** Get pointer to the event control bank
      LEVC = LZBS(IPEVC)
      LEMCR = LZBS(LEVC-ILEMCR)

      NSBEVT = NZBANK(IXSTOR,LZBS(LEMCR-ILEHLS))

      LEVT = LZBS(LEMCR-ILEHLS)
      LCUR = LZBS(LEMCR-ILEHCR)
      
      ISBEVT = 1

  100 CONTINUE

      IF(LEVT.EQ.LCUR) THEN

         IF(ISBEVT.LE.0 .OR. ISBEVT.GT.NSBEVT) 
     &      ISBEVT = -1
         RETURN

      END IF
            
      LEVT = LZBS(LEVT)
      ISBEVT = ISBEVT+1

      IF(ISBEVT.LE.NSBEVT) GO TO 100

      ISBEVT = -1

      RETURN
      END
