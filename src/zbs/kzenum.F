************************************************************************
*
*     --------------------------------
      SUBROUTINE KZENUM(NEVENT,IEVENT)
*     --------------------------------
*
*     (Purpose)
*        Returns the total number of events and
*        the event number of the current event
*        in the event buffer
*
*     (Inputs)
*
*     (Output)
*        NEVENT ; Total number of events
*        IEVENT ; Event number 
*
*     (Creation Date)
*        1994.06.15 : First version by J.Kanzaki.
*
************************************************************************


#include "zbssys.h"
#include "zbs.h"


***** Get pointer to the event control bank
      LEVC = LZBS(IPEVC)

      NEVENT = NZBANK(IXSTOR,LZBS(LEVC-ILEMLS))

      LEVT = LZBS(LEVC-ILEMLS)
      LCUR = LZBS(LEVC-ILEMCR)
      
      IEVENT = 1

  100 CONTINUE

      IF(LEVT.EQ.LCUR) THEN

         IF(IEVENT.LE.0 .OR. IEVENT.GT.NEVENT) 
     &      IEVENT = -1
         RETURN

      END IF
            
      LEVT = LZBS(LEVT)
      IEVENT = IEVENT+1

      IF(IEVENT.LE.NEVENT) GO TO 100

      IEVENT = -1

      RETURN
      END
