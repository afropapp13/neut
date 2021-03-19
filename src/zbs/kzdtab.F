************************************************************************
*
*     ----------------------
      SUBROUTINE KZDTAB(LUN)
*     ----------------------
*
*     (Purpose)
*        Dump contents of the bank directory table
*
*     (Inputs)
*
*     (Output)
*
*     (Creation Date)
*        1993.05.19 : First version by J.Kanzaki.
*        1997.8.26 : Fix byte order problem by C.McGrew
*
************************************************************************


#include "zbssys.h"
#include "zbs.h"

      CHARACTER*32 CNAME

      INTEGER*4 INAME(8),LNAME(8)
#ifdef INTEGER_COPY
      EQUIVALENCE (INAME(1),CNAME(1:4))
#endif

C     Get pointers to the bank directory

      LDIR=LZBS(IPDIRT)

      DO 10 I=1,8
      LNAME(I)=LZBS(LDIR-I)
   10 CONTINUE

      LPNT=LZBS(LDIR-ILBPNT)
      LEXF=LZBS(LDIR-ILBEXF)
      LIOF=LZBS(LDIR-ILBIOF)
      LSRT=LZBS(LDIR-ILBSRT)

      NBANK=IZBS(LDIR+INBANK)
      IF(NBANK.LE.0) RETURN

      WRITE(LUN,600)
  600 FORMAT(2(/),1H ,5X,'<< BANK DIRECTORY >>'//
     &            1H ,5X,'BANK NAME ',25X,'  PNT  EXF  IOF  SRT',/)

      DO 100 IBANK=1,NBANK

      DO 110 I=1,8
      INAME(I)=IZBS(LNAME(I)+IBANK)
  110 CONTINUE

#ifndef INTEGER_COPY
      CALL KZITOC(INAME,CNAME,8)
#endif

      IPNT = IZBS(LPNT+IBANK)
      IEXF = IZBS(LEXF+IBANK)
      IIOF = IZBS(LIOF+IBANK)

      DO 120 JBANK=1,NBANK
      JSRT = IZBS(LSRT+JBANK)
      IF(JSRT.EQ.IBANK) THEN
         ISRT = JBANK
         GO TO 130
      END IF
  120 CONTINUE

      CALL MSGOUT('KZDTAB','Bank directory matching error',1)
      ISRT = 0

  130 CONTINUE

      WRITE(LUN,610) IBANK,CNAME,IPNT,IEXF,IIOF,ISRT
  610 FORMAT(1H ,I3,2X,A32,3X,4I5)

  100 CONTINUE

      WRITE(LUN,620)
  620 FORMAT(2(/))

      RETURN
      END
