************************************************************************
*
*     ---------------------------------------
      SUBROUTINE KZCRD0(IBANK,CNAME,LDTR,IERR)
*     ---------------------------------------
*
*     (Purpose)
*       Create data record.
*
*     (Input)
*       LBIR  ; Address of the bank index record of the current bank
*       CNAME ; Name of the bank
*
*     (Output)
*       LDTR  ; Address of the created data record
*
*     (Creation Date and Author)
*       1991.12.25 ; First version by J.Kanzaki
*
************************************************************************


#include "zbs.h"
#include "zbssys.h"

      CHARACTER*(*) CNAME

      PARAMETER (NLDTR= 0,NSDTR= 0,IODTR=7)

     
      IERR = 0

      NDDTR = LDTRC0

      IXEVTD = IXDIV(2)

      CALL MZBOOK(IXEVTD,LDTR,LZBS(IBANK),-1,
     &            CNAME(1:4),NLDTR,NSDTR,NDDTR,IODTR,0)

      RETURN
      END
