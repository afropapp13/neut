************************************************************************
*
*     -----------------------------
      SUBROUTINE KZBCR0(CNAME,IERR)
*     -----------------------------
*
*     (Purpose)
*        Create new ZBS bank.
*
*     (Inputs)
*        CNAME  ; Name of the new bank ( <= 32 characters ).
*
*     (Output)
*        IERR      ; Error flag.
*             =  0 ; NORMAL RETURN
*
*     (Creation Date)
*        1991.12.25 : First version by J.Kanzaki.
*
************************************************************************


      CHARACTER*(*) CNAME
      CHARACTER*32 CNAME1

#include "zbssys.h"
#include "zbs.h"


      CALL KZCRMH

C==== SEARCH FOR BANK "CNAME" FROM THE DIRECTORY

      CNAME1 = CNAME

      CALL KZPDIR(CNAME1,IPDIR)

      LDIR = LZBS(IPDIRT)
      LEXF=LZBS(LDIR-ILBEXF)
      IEXIF=IZBS(LEXF+IPDIR)

      IF(IEXIF.NE.0) THEN
         IERR=1
         RETURN
      END IF

C     STORE INFORMATIONS ON THE NEW BANK INTO THE DIRECTORY

      LSRT=LZBS(LDIR-ILBSRT)

      LPNT=LZBS(LDIR-ILBPNT)
      LIOF=LZBS(LDIR-ILBIOF)

CCCCC PRINT *,'LPNT = ',LPNT
CCCCC PRINT *,'LIOF = ',LIOF

      IERR=0

      IZBS(LEXF+IPDIR)=1
      IZBS(LIOF+IPDIR)=1

C==== CREATE BANK RECORDS

      IBANK=IZBS(LPNT+IPDIR)

      CALL KZCRB0(IBANK,CNAME1,LBIR)

      CALL KZCRD0(IBANK,CNAME1,LDTR,IERDTR)

      RETURN
      END
