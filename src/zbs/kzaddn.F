************************************************************************
*
*     -------------------------------------
      SUBROUTINE KZADDN(IPFIND,CNAME,IPDIR)
*     -------------------------------------
*
*     (Purpose)
*        Add new entry into the directory table of ZBS.
*
*     (Inputs)
*        IPFIND ; The position where the entry should be added.
*        CNAME  ; Name of the entry ( <= 32 characters).
*
*     (Output)
*        IPDIR  ; Pointer to the new entry.
*
*     (Creation Date)
*        1991.12.25 : First version by J.Kanzaki.
*        1997.8.26 : Fix byte order problem by C.McGrew
*
************************************************************************

      CHARACTER*(*) CNAME
      CHARACTER*32 CNAME1

      INTEGER*4 INAME(8),LNAME(8)
#ifdef INTEGER_COPY
      EQUIVALENCE (INAME(1),CNAME1(1:4))
#endif

#include "zbssys.h"
#include "zbs.h"


      LDIR = LZBS(IPDIRT)

      LSRT=LZBS(LDIR-ILBSRT)

      NBANK=IZBS(LDIR+INBANK)
      NBANK=NBANK+1
      LBDIR=IZBS(LDIR+INBTAB)
      IF(LBDIR.LT.NBANK) THEN
         CALL KZEXDI(ICODE)
         LDIR = LZBS(IPDIRT)
      END IF
      IZBS(LDIR+INBANK)=NBANK

      CNAME1=CNAME
#ifndef INTEGER_COPY
      CALL KZCTOI(CNAME1,INAME,32)
#endif

      DO 10 I=1,8
      LNAME(I)=LZBS(LDIR-I)
   10 CONTINUE

      DO 20 I=1,8
      IZBS(LNAME(I)+NBANK)=INAME(I)
   20 CONTINUE

      IPSORT=ABS(IPFIND)
      IF(IPSORT.LT.NBANK) THEN
         DO 30 I=NBANK,IPSORT+1,-1
         IZBS(LSRT+I)=IZBS(LSRT+I-1)
   30    CONTINUE
      END IF
      IZBS(LSRT+IPSORT)=NBANK
      IPDIR=NBANK

      RETURN
      END
