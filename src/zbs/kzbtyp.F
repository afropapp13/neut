************************************************************************
*
*     ------------------------------
      SUBROUTINE KZBTYP(CNAME,ITYPE)
*     ------------------------------
*
*     (Purpose)
*       Get the type of the specified bank.
*
*     (Input)
*       CNAME  ; The name of the bank
*
*     (Output)
*       ITYPE  ; The type of the bank
*
*     (Creation Date and Author)
*       1994.07.27 ; First version by J.Kanzaki
*
************************************************************************


#include "zbs.h"
#include "zbssys.h"

      CHARACTER*(*) CNAME


      CALL KZBLOC(CNAME, LBIR)

      IF(LBIR.LE.0) RETURN

      ITYPE=IZBS(LBIR+IBKTYP)

      RETURN
      END
