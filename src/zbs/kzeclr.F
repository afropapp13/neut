************************************************************************
*
*     -----------------
      SUBROUTINE KZECLR
*     -----------------
*
*     (Purpose)
*        Delete (ZDROP'S and one ZGARB) all existing banks
*        and clear bank directory
*
*     (Creation Date)
*        1992.06.04 : First version by J.Kanzaki.
*
************************************************************************

#include "zbs.h"
#include "zbssys.h"

      IXEVTD = IXDIV(2)

      LDIR = LZBS(IPDIRT)

      NBANK = IZBS(LDIR+INBANK)
      LPNT=LZBS(LDIR-ILBPNT)
      LEXF=LZBS(LDIR-ILBEXF)
      LIOF=LZBS(LDIR-ILBIOF)

      LEVC = LZBS(IPEVC)

      CALL MZDROP(IXSTOR,LZBS(LEVC-ILEMLS),'L')

      LZBS(LEVC-ILEMCR) = 0

      DO 100 I=1,NBANK
      IF(IZBS(LEXF+I).EQ.0) GO TO 100
      IEBANK=IZBS(LPNT+I)
      IF(IEBANK.LE.0) GO TO 100
      IPBANK=LZBS(IEBANK)
      IF(IPBANK.LE.0) GO TO 100
      CALL MZDROP(IXSTOR,IPBANK,' ')
      IZBS(LEXF+I)=0
      IZBS(LIOF+I)=0
  100 CONTINUE

      CALL MZGARB(IXEVTD,0)

      RETURN
      END
