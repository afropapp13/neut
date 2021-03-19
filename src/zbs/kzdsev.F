************************************************************************
*
*     ----------------------
      SUBROUTINE KZDSEV(LUN)
*     ----------------------
*
*     (Purpose)
*        Dump information on the sub-events
*
*     (Inputs)
*
*     (Output)
*
*     (Creation Date)
*        1994.06.15 : First version by J.Kanzaki.
*
************************************************************************


#include "zbssys.h"
#include "zbs.h"

      CHARACTER*1 CCUR,CMRK
      CHARACTER*1 CSTR/'*'/,CPLS/'+'/,CBLK/' '/

***** Get pointer to the event control bank
      LEVC = LZBS(IPEVC)
      LECR = LZBS(LEVC-ILEMCR)
      LCUR = LZBS(LECR-ILEHCR)

      CALL KZMNUM(NSBEVT,ISBEVT)

      WRITE(LUN,600) NSBEVT
  600 FORMAT(2(/),5X,'<< Sub-event Table >>'//
     &            5X,'No. of sub-events in memory = ',I5/)

      IF(NSBEVT.LE.0) RETURN

      WRITE(LUN,'(5X,5X,''      Address''/)')

      LEHD = LZBS(LECR-ILEHLS)
      ISBEVT = 0

 1000 CONTINUE

      ISBEVT = ISBEVT+1

      IF(LEHD.EQ.LCUR) THEN
         CCUR = CPLS
      ELSE
         CCUR = CBLK
      END IF

C     IF(IZBS(LEMB+IEVMRK).EQ.1) THEN
C        CMRK = CSTR
C     ELSE
C        CMRK = CBLK
C     END IF

C     PRINT *,'KZDEVT;LZBS(LEMB-ILEHCR) = ',LZBS(LEMB-ILEHCR)
C     PRINT *,'KZDEVT;LZBS(LEMB-ILEHLS) = ',LZBS(LEMB-ILEHLS)

      WRITE(LUN,'(5X,I4,2X,A1,I10)') ISBEVT,CCUR,LEHD

      LEHD = LZBS(LEHD)
      IF(LEHD.GT.0) GO TO 1000

      WRITE(LUN,'(2(/))')

      RETURN
      END
