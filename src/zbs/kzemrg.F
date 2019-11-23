************************************************************************
*
*     -----------------
      SUBROUTINE KZEMRG
*     -----------------
*
*     (Purpose)
*        Merge marked events into the current event
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

      LEMCR = LZBS(LEVC-ILEMCR)
      LSTCR = LZLAST(IXSTOR,LZBS(LEMCR-ILEHLS))
      LEHCR = LZBS(LEMCR-ILEHLS)

      LEMB = LZBS(LEVC-ILEMLS)
      IEVENT = 0

 1000 CONTINUE

      IEVENT = IEVENT+1

      IF(LEMB.NE.LEMCR .AND.
     &   IZBS(LEMB+IEVMRK).EQ.1) THEN

         LEHD = LZBS(LEMB-ILEHLS)
         CALL ZSHUNT(IXSTOR,LEHD,LSTCR,0,1)

      END IF

      LEMB = LZBS(LEMB)
      IF(LEMB.GT.0) GO TO 1000


      RETURN
      END
