************************************************************************
*
*     ------------------------
      SUBROUTINE KZEXDI(ICODE)
*     ------------------------
*
*     (Purpose)
*       Expand the bank directory table.
*
*     (Input)
*       none
*
*     (Output)
*       ICODE ; Error code.
*
*     (Creation Date and Author)
*       1991.12.25 ; First version by J.Kanzaki
*
************************************************************************

#include "zbs.h"
#include "zbssys.h"


      ICODE=0

      LDIR=LZBS(IPDIRT)

      DO 10 I=1,LLNKHD

C DWC    An index in a bank may not be passed as an argument to MZPUSH
C DWC    see Zebra manual...
C DWC    Also, the bank index LDIR may change if a garbage collection
C DWC    occurs on this call, and therefore must be restored after each
C DWC    call in the loop
C      CALL MZPUSH(IXSTOR,LZBS(LDIR-I),0,LBDIR1,' ')
      LPUSH = LZBS(LDIR-I)
      CALL MZPUSH(IXSTOR,LPUSH,0,LBDIR1,' ')
      LDIR = LZBS(IPDIRT)
      IF (LZBS(LDIR-I).NE.LPUSH) THEN
         print *,'kzexdi: Whoa - index LDIR-I incorrect after MZPUSH!'
         LZBS(LDIR-I) = LPUSH
      ENDIF
   10 CONTINUE

      IZBS(LDIR+INBTAB)=IZBS(LDIR+INBTAB)+LBDIR1

      RETURN
      END

