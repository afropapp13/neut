************************************************************************
*
*     -----------------------------
      SUBROUTINE KZBCRE(CNAME,IERR)
*     -----------------------------
*
*     (Purpose)
*        Create new ZBS bank.
*
*     (Inputs)
*        CNAME  ; Name of the new bank ( <= 32 characters ).
*
*     (Output)
*        IPBANK ; Pointer to the created bank.
*        IERR      ; Error flag.
*             =  0 ; NORMAL RETURN
*
*     (Creation Date)
*        1991.12.25 : First version by J.Kanzaki.
*
************************************************************************


      CHARACTER*(*) CNAME

#include "zbssys.h"
#include "zbs.h"


      CALL KZBCR0(CNAME,IERR)

      RETURN
      END
