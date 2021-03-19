      SUBROUTINE KZLDAT(CNAME,LDAT)

#include "zbs.h"
#include "zbssys.h"

      CHARACTER*(*) CNAME
      INTEGER LDAT

      CALL KZBLOC(CNAME,LBIR)

      LDTR = LZBS(LBIR-ILDATR)

      LDAT = LDTR

      RETURN
      END

