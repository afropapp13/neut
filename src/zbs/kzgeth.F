************************************************************************
*
*     ------------------------
      SUBROUTINE KZGETH(IDATA)
*     -------------------------
*
*     (Purpose)
*       Get the contents of the data record of Event Mother Bank
*
*     (Input)
*
*     (Output)
*       IDATA ; Contents of the data record of Event Mother Bank
*
*     (Creation Date and Author)
*       1994.07.19 ; First version by J.Kanzaki
*
************************************************************************


#include "zbs.h"
#include "zbssys.h"

      INTEGER*4 IDATA


      IDATA = 0
C==== Search for Event Mother Bank

      LEVC = LZBS(IPEVC)
      LEMB = LZBS(LEVC-ILEMCR)

      IF(LEMB.LE.0) RETURN

      IDATA = IZBS(LEMB+1)

      RETURN
      END
