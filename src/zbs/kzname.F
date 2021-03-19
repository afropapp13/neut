************************************************************************
*
*     ----------------------
      FUNCTION KZNAME(CNAME)
*     ----------------------
*
*     (Purpose)
*       Searches for the entry with a specified name, CNAME,
*       from the bank directory table.
*
*     (Input)
*       CNAME ; Name of the bank
*
*     (Output)
*       KZNAME ; Pointer to the bank sort table.
*              ; If the name, CNAME, has been already registered into the
*              ; name table,
*              ;  = position of that bank in the sort table
*              ; If not,
*              ;  = - ( position in the sort table just after the bank,
*              ;        CNAME, to be stored )
*
*     (Creation Date and Author)
*       1991.12.25 ; First version by J.Kanzaki
*        1997.8.26 : Fix byte order problem by C.McGrew
*
************************************************************************

#include "zbs.h"
#include "zbssys.h"

      CHARACTER*(*) CNAME
      CHARACTER*32 CNAME1
      INTEGER*4 INAME(8),LNAME(8)
#ifdef INTEGER_COPY
      EQUIVALENCE (INAME(1),CNAME1(1:4))
#endif

      KZNAME=0

      LDIR=LZBS(IPDIRT)

      NBANK=IZBS(LDIR+INBANK)
C     PRINT *,'KZNAME;NBANK = ',NBANK
      IF(NBANK.LE.0) THEN
         KZNAME=-1
         RETURN
      END IF

      CNAME1=CNAME
#ifndef INTEGER_COPY
      CALL KZCTOI(CNAME1,INAME,32)
#endif

C     DO 5 I=1,8
CCCCC PRINT *,'KZNAME;I,INAME = ',I,INAME(I)
C   5 CONTINUE

      DO 10 I=1,8
      LNAME(I)=LZBS(LDIR-I)
   10 CONTINUE
      LSRT=LZBS(LDIR-ILBSRT)

      IPLOW=1
      IPHIG=NBANK

   20 IF(IPLOW.GT.IPHIG) GO TO 100
      IP=(IPLOW+IPHIG)/2
C     PRINT *,'KZNAME;IP = ',IP
      K=IZBS(LSRT+IP)
      ICOMP=KZCMPN(LNAME,K,INAME)
C     PRINT *,'KZNAME;ICOMP = ',ICOMP
      IF(ICOMP.LE.0) IPHIG=IP-1
      IF(ICOMP.GE.0) IPLOW=IP+1

      GO TO 20

  100 K=IZBS(LSRT+IP)
      ICOMP=KZCMPN(LNAME,K,INAME)
      IF(ICOMP.EQ.0) THEN
         KZNAME=IP
      ELSE IF(ICOMP.LT.0) THEN
         KZNAME=-IP
      ELSE
         KZNAME=-IP-1
      END IF

      RETURN
      END
