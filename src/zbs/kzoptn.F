***********************************************************************
*
*     -----------------------
      SUBROUTINE KZOPTN(COPT)
*     -----------------------
*
*     (Purpose)
*       Set up various environments for ZBS
*
*     (History)
*       Jul. 14, 1994.    First version created.        J. Kanzaki
*
***********************************************************************


#include "zbs.h"
#include "zbssys.h"

      CHARACTER*(*) COPT

      IF(INDEX(COPT,'NO-ZEBRA-INIT').NE.0) THEN

         IFZEBI = 0

      END IF

      RETURN
      END
