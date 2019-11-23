************************************************************************
*     $Id: kzcrb3.F,v 1.1 1997/10/27 18:58:32 mcgrew Exp $
*     -------------------------------------------------
      SUBROUTINE KZCRB3(IBANK,CNAME,NLINKS,NWMAST,NSEGM,LSEGM)
*     -------------------------------------------------
*
*     (Purpose)
*        Create bank index record.
*
*     (Inputs)
*        IBANK  ; Entry position of bank directory
*        CNAME  ; Name of the entry ( <= 32 characters).
*        NLINKS ; The number of supporting links.
*        NWMAST ; The length of the fixed section.
*        NSEGM  ; The length of the variable section
*        LSEGM  ; Number of repeats in variable section.
*
*     (Output)
*
*     (Creation Date)
*        1997.10.21 : First version by C.McGrew
*
************************************************************************

#include "zbs.h"
#include "zbssys.h"

      INTEGER*4 LEMB,IBANK
      CHARACTER*(*) CNAME

      INTEGER*4 LBIR

      PARAMETER (NLBIR=1,NSBIR=1,IOBIR=2)

      CHARACTER*32 CNAME1

#ifdef INTEGER_COPY
      INTEGER*4 INAME(8)
      EQUIVALENCE (INAME(1),CNAME1(1:4))
#endif

      PARAMETER (ITYPE=1)


      CNAME1 = CNAME

CCCCC PRINT *,'KZCRBI;CNAME1 = ',CNAME1

      LEVC = LZBS(IPEVC)
      LEMB = LZBS(LEVC-ILEMCR)
      LEHD = LZBS(LEMB-ILEHCR)

      LLAST = LZLAST(IXSTOR,LZBS(LEHD-ILBILS))

      IXEVTD = IXDIV(2)
      NDBIR = LINDHD

CCCCC PRINT *,'KZCRBI;LLAST = ',LLAST

      IF(LLAST.NE.0) THEN

         CALL MZBOOK(IXEVTD,LZBS(IBANK),LLAST,0,
     &               CNAME(1:4),NLBIR,NSBIR,NDBIR,IOBIR,0)

      ELSE

         CALL MZBOOK(IXEVTD,LZBS(IBANK),LEHD,-ILBILS,
     &               CNAME(1:4),NLBIR,NSBIR,NDBIR,IOBIR,0)

      END IF

      LBIR = LZBS(IBANK)

#ifdef INTEGER_COPY
      DO 10 I=1,8
      IZBS(LBIR+I) = INAME(I)
   10 CONTINUE
#else
      CALL KZCTOI(CNAME1,IZBS(LBIR+1),32)
#endif

C**** space for a bank type

      IZBS(LBIR+IBKTYP)=ITYPE
      IZBS(LBIR+ILDTRC)=0
      IZBS(LBIR+ILFILL)=0
      IZBS(LBIR+INEXDT)=0
      IZBS(LBIR+ILDTR1)=0
      IF(NSEGM.GT.0) THEN
         IZBS(LBIR+INSEGM) = NSEGM
         IZBS(LBIR+INSGM1) = MAX(NINT(FSEGM1*REAL(NSEGM)),1)
      ELSE IF(LSEGM.GT.0) THEN
         IZBS(LBIR+INSEGM) = NSEGM0
         IZBS(LBIR+INSGM1) = NINT(FSEGM1*REAL(NSEGM0))
      ELSE
         IZBS(LBIR+INSEGM) = 0
         IZBS(LBIR+INSGM1) = 0
      END IF
      IZBS(LBIR+INSGMX)=-1
      IZBS(LBIR+IISLST)=0
      IZBS(LBIR+INEXSG)=0

      IZBS(LBIR+INMAST) = NWMAST
      IZBS(LBIR+ILSEGF) = LSEGM

      RETURN
      END
