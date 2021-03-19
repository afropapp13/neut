************************************************************************
*
*     ------------------------
      SUBROUTINE KZBDEL(CNAME)
*     ------------------------
*
*     (Purpose)
*        Delete bank specified by CNAME.
*
*     (Creation Date)
*        1992.12.22 : First version by J.Kanzaki.
*
************************************************************************

      CHARACTER*(*) CNAME

#include "zbs.h"
#include "zbssys.h"


      LDIR = LZBS(IPDIRT)

      LPNT=LZBS(LDIR-ILBPNT)
      LEXF=LZBS(LDIR-ILBEXF)
      LIOF=LZBS(LDIR-ILBIOF)
      LSRT=LZBS(LDIR-ILBSRT)

      IFIND=KZNAME(CNAME)

      IPDIR=IZBS(LSRT+IFIND)
      IF(IPDIR.LE.0) RETURN

      IEBANK=IZBS(LPNT+IPDIR)
      IF(IEBANK.LE.0) RETURN

      IPBANK=LZBS(IEBANK)
      IF(IPBANK.LE.0) RETURN

      CALL MZDROP(IXSTOR,IPBANK,' ')

C     bug fixed by M.Shiozawa, 98/7/30
C      IZBS(IEBANK)=0
      LZBS(IEBANK)=0
      IZBS(LEXF+IPDIR)=0
      IZBS(LIOF+IPDIR)=0

      RETURN
      END
