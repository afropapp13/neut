************************************************************************
*
*     ------------------------------
      SUBROUTINE KZNSEG(CNAME,NSEGM)
*     ------------------------------
*
*     (Purpose)
*       Get the number of segments of the specified bank.
*
*     (Input)
*       CNAME  ; The name of the bank
*
*     (Output)
*       NSEGM  ; No. of segments
*
*     (Creation Date and Author)
*       1993.05.19 ; First version by J.Kanzaki
*
************************************************************************


#include "zbs.h"
#include "zbssys.h"

      CHARACTER*(*) CNAME


      CALL KZBLOC(CNAME, LBIR)

      NSEGM = -1

      IF(LBIR.LE.0) RETURN

      NSEGM = IZBS(LBIR+INSGMX)

      RETURN
      END
