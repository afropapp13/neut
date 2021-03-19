************************************************************************
*
*     ------------------------
      SUBROUTINE KZDSTO(IFLAG)
*     ------------------------
*
*     (Purpose)
*       Print out contents of the specified bank.
*
*     (Input)
*       CNAME ; Name of the bank
*       IFLAG ; 
*
*     (Output)
*       none
*
*     (Creation Date and Author)
*       1992.12.18 ; First version by J.Kanzaki
*
************************************************************************


#include "zbssys.h"
#include "zbs.h"

      INTEGER*4 LTEXT
      CHARACTER*50 CTEXT/'Dump the store structure'/


      LTEXT = KZLENG(CTEXT)

      CALL DZSTOR(CTEXT(1:LTEXT),IXSTOR)

      RETURN
      END
