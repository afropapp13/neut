************************************************************************
*
*     -----------------------------------
      SUBROUTINE KZNSSG(CNAME,ISEGM,NSSEG)
*     -----------------------------------
*
*     (Purpose)
*       Get the number of sub-segments of the specified segment
*       for 2-dim. fixed banks.
*
*     (Input)
*       CNAME ; The name of the bank
*       ISEGM ; Segment #
*
*     (Output)
*       NSSEG  ; No. of sub-segments
*
*     (Creation Date and Author)
*       1993.07.16 ; First version by J.Kanzaki
*
************************************************************************


#include "zbs.h"
#include "zbssys.h"

      CHARACTER*(*) CNAME


      CALL KZBLOC(CNAME, LBIR)

      NSSEG = -1

      IF(LBIR.LE.0) RETURN

      ITYPE = IZBS(LBIR+IBKTYP)

      IF(ITYPE.NE.2) RETURN

      NSEGM=IZBS(LBIR+INSEGM)

      IF(ISEGM.LE.0 .OR. ISEGM.GT.NSEGM) RETURN

      LIND=LBIR+LINDHD+ISEGM*LINDSG
      NSSEG=IZBS(LIND+2)

      RETURN
      END
