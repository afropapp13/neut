************************************************************************
*
*     ------------------------
      SUBROUTINE KZDDIV(IFLAG)
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
      CHARACTER*50 CTEXT/'Dump the division structure'/

      LTEXT = KZLENG(CTEXT)

      CALL DZSNAP(CTEXT(1:LTEXT),IXDIV(2),'F')
C     CALL DZSNAP(CTEXT(1:LTEXT),IXDIV(3),'F')

      RETURN
      END
