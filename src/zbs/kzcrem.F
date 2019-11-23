************************************************************************
*
*     ----------------------------
      SUBROUTINE KZCREM(LEVC,LEMB)
*     ----------------------------
*
*     (Purpose)
*       Create event mother bank.
*
*     (Input)
*       LEVC  ; Address of Event Control Bank
*
*     (Output)
*       LEMB  ; Address of created event mother bank
*
*     (Creation Date and Author)
*       1991.12.25 ; First version by J.Kanzaki
*
************************************************************************


#include "zbssys.h"
#include "zbs.h"

      CHARACTER*4 CNEMB/'*EMB'/
      PARAMETER (NLEMB= 2,NSEMB= 1,IOEMB=2)


      IXEVTD = IXDIV(2)

      LLAST = LZLAST(IXSTOR,LZBS(LEVC-ILEMLS))

      NDEMB = LEMBDT

      IF(LLAST.NE.0) THEN

         CALL MZBOOK(IXEVTD,LEMB,LLAST,0,
     &               CNEMB,NLEMB,NSEMB,NDEMB,IOEMB,0)

      ELSE

         CALL MZBOOK(IXEVTD,LEMB,LEVC,-ILEMLS,
     &               CNEMB,NLEMB,NSEMB,NDEMB,IOEMB,0)

      END IF
CCCCC PRINT *,'KZCREM;LEMB = ',LEMB

C**** Temporary for single event version

      LZBS(LEVC-ILEMCR) = LEMB

      IZBS(LEMB+IEVMRK) = 0

      RETURN
      END
