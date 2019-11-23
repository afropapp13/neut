************************************************************************
*     ---------------------
      SUBROUTINE KZEND(LUN)
*     ---------------------
*
*     (Purpose)
*       Close the file specified by the LUN.
*
*     (Input)
*       LUN ; Logical unit number of the file to be closed.
*
*     (Output)
*       none
*
*     (Creation Date and Author)
*       1991.12.25 ; First version by J.Kanzaki
*
************************************************************************

#include "zbssys.h"
#include "zbs.h"

      CALL FZENDO(LUN,'T')

      RETURN
      END
