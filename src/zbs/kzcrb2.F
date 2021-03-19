************************************************************************
*
*     -------------------------------------------------
      SUBROUTINE KZCRB2(IBANK,CNAME,NWMAST,NSEGM,LSEGM)
*     -------------------------------------------------
*
*     (Purpose)
*        Create bank index record.
*
*     (Inputs)
*        LEMB   ; Address of the event mother bank.
*        IBANK  ; Entry position of bank directory
*        CNAME  ; Name of the entry ( <= 32 characters).
*
*     (Output)
*        LBIR   ; Address of the created bank index record.
*
*     (Creation Date)
*        1992.12.18 : First version by J.Kanzaki.
*        1997.8.26 : Fix byte order problem by C.McGrew
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

      PARAMETER (ITYPE=2)


      CNAME1 = CNAME

      LEVC = LZBS(IPEVC)
      LEMB = LZBS(LEVC-ILEMCR)
      LEHD = LZBS(LEMB-ILEHCR)

      LLAST = LZLAST(IXSTOR,LZBS(LEHD-ILBILS))

      IXEVTD = IXDIV(2)
      NDBIR = LINDHD + LINDSG*(NSEGM+1)

      IF(LLAST.NE.0) THEN

         CALL MZBOOK(IXEVTD,LZBS(IBANK),LLAST,0,
     &               CNAME(1:4),NLBIR,NSBIR,NDBIR,IOBIR,0)

      ELSE

         CALL MZBOOK(IXEVTD,LZBS(IBANK),LEHD,-ILBILS,
     &               CNAME(1:4),NLBIR,NSBIR,NDBIR,IOBIR,0)

      END IF

      LBIR = LZBS(IBANK)

      DO 100 I=0,NSEGM
      LIND=LBIR+LINDHD+I*LINDSG
      IF(I.EQ.0) THEN
         IZBS(LIND+1)=0
         IZBS(LIND+2)=0
      ELSE
         IZBS(LIND+1)=-1
         IZBS(LIND+2)=0
      END IF
  100 CONTINUE

#ifdef INTEGER_COPY
      DO 10 I=1,8
      IZBS(LBIR+I) = INAME(I)
   10 CONTINUE
#else
      CALL KZCTOI(CNAME1,IZBS(LBIR+1),32)
#endif

C**** space for a bank type

      IZBS(LBIR+IBKTYP)=ITYPE
      IZBS(LBIR+ILDTRC)=NWMAST+LDTRC0
      IZBS(LBIR+ILFILL)=NWMAST
      IZBS(LBIR+INEXDT)=0
      IZBS(LBIR+ILDTR1)=NINT(FDTRC1*REAL(LDTRC0))
      IZBS(LBIR+INSEGM) = NSEGM
      IZBS(LBIR+INSGMX)=-1
      IZBS(LBIR+IISLST)=0
      IZBS(LBIR+INEXSG)=0
      IZBS(LBIR+INSGM1)=0

      IZBS(LBIR+INMAST) = NWMAST
      IZBS(LBIR+ILSEGF) = LSEGM

      RETURN
      END
