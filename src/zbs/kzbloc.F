************************************************************************
*     -------------------------------------
      SUBROUTINE KZBLOC(CNAME,LBANK)
*     ------------------------------
*
*     (Purpose)
*        Locate bank with "CNAME"
*
*     (Inputs)
*        CNAME  ; Bank name to be located.
*
*     (Output)
*        LBANK  ; Pointer to the bank.
*
*     (Creation Date)
*        1991.12.25 : First version by J.Kanzaki.
*        1996.04.26 : Add declarations. J.Kanzaki
*
************************************************************************

#include "zbsdeclr.h"
#include "zbssys.h"
#include "zbs.h"

      CHARACTER*(*) CNAME
      INTEGER LBANK

      INTEGER IPDIR,LDIR,LPNT,LEXF
      INTEGER IEXIF,IEBANK,IPBANK


      LBANK=0
      CALL KZPDIR(CNAME,IPDIR)

      LDIR = LZBS(IPDIRT)

      LEXF = LZBS(LDIR-ILBEXF)
      IEXIF = IZBS(LEXF+IPDIR)

      IF(IEXIF.EQ.0) RETURN

      LPNT=LZBS(LDIR-ILBPNT)
      IEBANK=IZBS(LPNT+IPDIR)

      IF(IEBANK.LE.0) RETURN

      IPBANK=LZBS(IEBANK)

      IF(IPBANK.EQ.0) RETURN

      LBANK = IPBANK

      RETURN
      END
