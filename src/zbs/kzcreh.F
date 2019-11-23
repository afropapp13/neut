************************************************************************
*
*     ----------------------------
      SUBROUTINE KZCREH(LEMB,LEHD)
*     ----------------------------
*
*     (Purpose)
*       Create event header bank.
*
*     (Input)
*       LEMB  ; Address of Event Mother Bank
*
*     (Output)
*       LEHD  ; Address of created event header bank
*
*     (Creation Date and Author)
*       1991.12.25 ; First version by J.Kanzaki
*
************************************************************************


#include "zbssys.h"
#include "zbs.h"

      CHARACTER*4 CNEHD/'*EHD'/
      PARAMETER (NLEHD= 1,NSEHD= 1,NDEMB=1,IOEHD=2)


      IXEVTD = IXDIV(2)

      LLAST = LZLAST(IXSTOR,LZBS(LEMB-ILEHLS))

      IF(LLAST.NE.0) THEN

         CALL MZBOOK(IXEVTD,LEHD,LLAST,0,
     &               CNEHD,NLEHD,NSEHD,NDEHD,IOEHD,0)

      ELSE

         CALL MZBOOK(IXEVTD,LEHD,LEMB,-1,
     &               CNEHD,NLEHD,NSEHD,NDEHD,IOEHD,0)

      END IF
CCCCC PRINT *,'KZCREM;LEMB = ',LEMB

C**** Temporary for single event version

      LZBS(LEMB-ILEHCR) = LEHD

      RETURN
      END
