************************************************************************
*     $Id: kzcrd3.F,v 1.1 1997/10/27 18:58:33 mcgrew Exp $
*     -----------------------------------------------------
      SUBROUTINE KZCRD3(IBANK,CNAME,NLINKS,NWMAST,NSEGM,LSEGM,IOX)
*     -----------------------------------------------------
*
*     (Purpose)
*       Create data record.
*
*     (Input)
*       IBANK  ; Entry of the Address of the bank in the index table
*       CNAME  ; Name of the bank
*       NLINKS ; The number of supporting links
*       NWMAST ; Length of the master record
*       NSEGM  ; No. of segments
*       LSEGM  ; Length of the segment
*       IOX    ; I/O characteristic word
*
*     (Output)
*
*     (Creation Date and Author)
*       1997.10.21; First version by C. McGrew
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
     &            CNAME(1:4),NLINKS,NLINKS,NDDTR,IOX,0)

      RETURN
      END

