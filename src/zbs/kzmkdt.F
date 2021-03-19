************************************************************************
*     -----------------
      SUBROUTINE KZMKDT
*     -----------------
*
*     (Purpose)
*       Create bank directory table for the current event
*
*     (Input)
*
*     (Output)
*
*     (Comment)
*
*     (Creation Date and Author)
*       1994.06.14 ; First version by J.Kanzaki
*       1996.04.26 : Add declarations. J.Kanzaki
*
************************************************************************

#include "zbsdeclr.h"
#include "zbssys.h"
#include "zbs.h"

      CHARACTER*32 CNAME

#ifdef INTEGER_COPY
      INTEGER INAME(8)
      EQUIVALENCE (INAME(1),CNAME(1:4))
#endif

      INTEGER LEVC,LEMB,LEHD,LBIR

      INTEGER LDIR,LEXF,LSRT,LPNT,LIOF

      INTEGER IPFIND,IPDIR,NBANK,IBANK

      INTEGER KZNAME

      INTEGER I

***** Get pointers to the bank directory table
      LDIR = LZBS(IPDIRT)

      LEXF=LZBS(LDIR-ILBEXF)
      IF(LEXF.LE.0) RETURN

      NBANK=IZBS(LDIR+INBANK)
c      PRINT *,'KZMKDT;NBANK = ',NBANK
      IF(NBANK.GT.0) THEN

         DO 10 IBANK=1,NBANK
         IZBS(LEXF+IBANK) = 0
   10    CONTINUE

      END IF

***** Get pointer to the event control bank

      LEVC = LZBS(IPEVC)

      LEMB = LZBS(LEVC-ILEMCR)
      LEHD = LZBS(LEMB-ILEHCR)

      LBIR = LZBS(LEHD-ILBILS)

      IF(LBIR.LE.0) RETURN

***** Read next record. This should be bank header or event trailer.
  100 CONTINUE

#ifdef INTEGER_COPY
      DO 110 I=1,8
      INAME(I) = IZBS(LBIR+I)
  110 CONTINUE
#else
      CALL KZITOC(IZBS(LBIR+1),CNAME,8)
#endif


#ifdef DEBUG_ZBS
      PRINT *,'KZMKDT;CNAME = ',CNAME
#endif
***** SEARCH FOR BANK "CNAME" FROM THE DIRECTORY
      IPFIND=KZNAME(CNAME)
      IF(IPFIND.EQ.0) RETURN

***** GET POINTERS TO THE SYSTEM MOTHER BANK AND THE BANK DIRECTOR
      LDIR = LZBS(IPDIRT)

      LEXF=LZBS(LDIR-ILBEXF)
      IF(LEXF.LE.0) RETURN
      LSRT=LZBS(LDIR-ILBSRT)
      IF(LSRT.LE.0) RETURN

***** "CNAME" is found in the directory table.
      IF(IPFIND.GT.0) THEN

         IPDIR=IZBS(LSRT+IPFIND)

***** IF BANK OF SPEICIFIED NAME "CNAME" IS NOT FOUND IN THE DIRECTORY,
*     CREATE NEW ENTRY IN THE BANK DIRECTORY.

      ELSE

        CALL KZADDN(IPFIND,CNAME,IPDIR)

      END IF

***** STORE INFORMATIONS ON THE NEW BANK INTO THE DIRECTORY

      LEXF=LZBS(LDIR-ILBEXF)
      LSRT=LZBS(LDIR-ILBSRT)

      LPNT=LZBS(LDIR-ILBPNT)
      IF(LPNT.LE.0) RETURN
      LIOF=LZBS(LDIR-ILBIOF)
      IF(LIOF.LE.0) RETURN

      IZBS(LEXF+IPDIR)=1
      IZBS(LIOF+IPDIR)=1

***** Pointer# to store the address of bank
      IBANK = IZBS(LPNT+IPDIR)

      LZBS(IBANK) = LBIR

      LBIR = LZBS(LBIR)


      IF(LBIR.GT.0) GO TO 100

      RETURN
      END
