************************************************************************
*
*     -----------------------------------------------------
      SUBROUTINE KZCRD1(IBANK,CNAME,NWMAST,NSEGM,LSEGM,IOX)
*     -----------------------------------------------------
*
*     (Purpose)
*       Create data record.
*
*     (Input)
*       IBANK  ; Entry of the Address of the bank in the index table
*       CNAME  ; Name of the bank
*       NWMAST ; Length of the master record
*       NSEGM  ; No. of segments
*       LSEGM  ; Legnth of th segment
*       IOX    ; I/O characteristic word
*
*     (Output)
*
*     (Creation Date and Author)
*       1991.12.25 ; First version by J.Kanzaki
*
************************************************************************


#include "zbs.h"
#include "zbssys.h"

      CHARACTER*(*) CNAME

      PARAMETER (NLDTR= 0,NSDTR= 0)


      LBIR = LZBS(IBANK)

      NSEGM1 = IZBS(LBIR+INSEGM)

      NDDTR = NWMAST + NSEGM1 * LSEGM 

      IXEVTD = IXDIV(2)

      CALL MZBOOK(IXEVTD,LDTR,LZBS(IBANK),-1,
     &            CNAME(1:4),NLDTR,NSDTR,NDDTR,IOX,0)

      RETURN
      END
