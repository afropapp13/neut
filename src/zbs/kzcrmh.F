************************************************************************
*
*     -----------------
      SUBROUTINE KZCRMH
*     -----------------
*
*     (Purpose)
*        Create Event Mother Bank and Event Header Bank
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


C==== Search for Event Mother Bank

      LEVC = LZBS(IPEVC)
      LEMB = LZBS(LEVC-ILEMCR)

      IF(LEMB.LE.0) THEN

C==== If EMB does not exist, create it.

         CALL KZCREM(LEVC,LEMB)

      END IF

      LEHD = LZBS(LEMB-ILEHCR)

      IF(LEHD.LE.0) THEN

C==== If EHD does not exist, create it.

         CALL KZCREH(LEMB,LEHD)

      END IF

      RETURN
      END
