************************************************************************
*
*     ----------------------
      SUBROUTINE KZDEVT(LUN)
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
*        1994.06.14 : First version by J.Kanzaki.
*
************************************************************************


#include "zbssys.h"
#include "zbs.h"

      CHARACTER*1 CCUR,CMRK
      CHARACTER*1 CSTR/'*'/,CPLS/'+'/,CBLK/' '/

***** Get pointer to the event control bank
      LEVC = LZBS(IPEVC)

      LCUR = LZBS(LEVC-ILEMCR)

      CALL KZENUM(NEVENT,IEVENT)

      WRITE(LUN,600) NEVENT
  600 FORMAT(2(/),5X,'<< Event Table >>'//
     &            5X,'No. of events in memory = ',I5/)

      IF(NEVENT.LE.0) RETURN

      WRITE(LUN,'(5X,5X,''      Address''/)')

      LEMB = LZBS(LEVC-ILEMLS)
      IEVENT = 0

 1000 CONTINUE

      IEVENT = IEVENT+1

      IF(LEMB.EQ.LCUR) THEN
         CCUR = CPLS
      ELSE
         CCUR = CBLK
      END IF

      IF(IZBS(LEMB+IEVMRK).EQ.1) THEN
         CMRK = CSTR
      ELSE
         CMRK = CBLK
      END IF

      NSBEVT = NZBANK(IXSTOR,LZBS(LEMB-ILEHLS))
C     PRINT *,'KZDEVT;LZBS(LEMB-ILEHCR) = ',LZBS(LEMB-ILEHCR)
C     PRINT *,'KZDEVT;LZBS(LEMB-ILEHLS) = ',LZBS(LEMB-ILEHLS)

      WRITE(LUN,'(5X,I4,2X,A1,2X,A1,I10,I5)') IEVENT,CCUR,CMRK,LEMB,
     &     NSBEVT

      LEMB = LZBS(LEMB)
      IF(LEMB.GT.0) GO TO 1000

      WRITE(LUN,'(2(/))')

      RETURN
      END
