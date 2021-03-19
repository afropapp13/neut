************************************************************************
*                                                                      *
*     -----------------
      SUBROUTINE KZCLDT
*     -----------------
*                                                                      *
*     (Purpose)                                                        *
*       Clear bank directory table
*                                                                      *
*     (Input)                                                          *
*                                                                      *
*     (Output)                                                         *
*                                                                      *
*     (Comment)                                                        *
*                                                                      *
*     (Creation Date and Author)
*       1994.07.06 ; First version by J.Kanzaki
*                                                                      *
************************************************************************

#include "zbssys.h"
#include "zbs.h"


***** Get pointers to the bank directory table
      LDIR = LZBS(IPDIRT)

      NBANK=IZBS(LDIR+INBANK)
      IF(NBANK.LE.0) RETURN

      LEXF=LZBS(LDIR-ILBEXF)
      IF(LEXF.LE.0) RETURN

      DO 100 IBANK=1,NBANK
      IZBS(LEXF+IBANK) = 0
  100 CONTINUE

      RETURN
      END
