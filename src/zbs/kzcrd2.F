************************************************************************
*
*     -----------------------------------------------------
      SUBROUTINE KZCRD2(IBANK,CNAME,NWMAST,NSEGM,LSEGM,IOX)
*     -----------------------------------------------------
*
*     (Purpose)
*       Create data record. (fixed 2-dim. case)
*
*     (Input)
*       CNAME  ; Name of the bank
*       NWMAST ; Length of the master record
*       NSEGM  ; No. of segments
*       LSEGM  ; Legnth of th segment
*
*     (Output)
*       LDTR  ; Address of the created data record
*
*     (Creation Date and Author)
*       1992.12.22 ; First version by J.Kanzaki
*
************************************************************************


#include "zbs.h"
#include "zbssys.h"

      CHARACTER*(*) CNAME

      PARAMETER (NLDTR= 0,NSDTR= 0)


      NDDTR = NWMAST + LDTRC0

      IXEVTD = IXDIV(2)

      CALL MZBOOK(IXEVTD,LDTR,LZBS(IBANK),-1,
     &            CNAME(1:4),NLDTR,NSDTR,NDDTR,IOX,0)

      RETURN
      END
