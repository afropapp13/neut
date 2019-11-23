************************************************************************
*
*     ----------------------
      SUBROUTINE KZFILE(LUN)
*     ----------------------
*
*     (Purpose)
*       Define files for data I/O
*
*     (Input)
*       LUN  ; Logical unit number
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

C     CALL FZFILE(LUN,0,'XO')
      CALL FZFILE(LUN,0,'XIO')
C     CALL FZFILE(LUN,1620,'XIO')

CCCCC PRINT *,'KZFILE;IQUEST(1) = ',IQUEST(1)

      RETURN
      END
