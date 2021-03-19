************************************************************************
*     ------------------------------
      SUBROUTINE KZPDIR(CNAME,IPDIR)
*     ------------------------------
*
*     (Purpose)
*       Returns entry point in the pointer table for a given
*       bank with a name "CNAME".
*
*     (Input)
*       CNAME ; Name of the bank.
*
*     (Output)
*       IPDIR ; Entry point of the pointer table.
*
*     (Creation Date and Author)
*       1994.06.15 ; First version by J.Kanzaki
*
************************************************************************


#include "zbssys.h"
#include "zbs.h"

      CHARACTER*(*) CNAME


      IPDIR=0

      IPFIND=KZNAME(CNAME)

      IF(IPFIND.EQ.0) RETURN

      LDIR = LZBS(IPDIRT)

      LPNT=LZBS(LDIR-ILBPNT)
      LSRT=LZBS(LDIR-ILBSRT)

      IF(IPFIND.GT.0) THEN

         IPDIR=IZBS(LSRT+IPFIND)

      ELSE

         CALL KZADDN(IPFIND,CNAME,IPDIR)

      END IF

      RETURN
      END
