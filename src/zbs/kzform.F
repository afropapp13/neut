************************************************************************
*
*     -----------------------------------------------
      SUBROUTINE KZFORM(CNAME,CFORMM,CFORMS,IOX,IERR)
*     -----------------------------------------------
*
*     (Purpose)
*       Create format character word.
*
*     (Input)
*       CNAME  ; Name of the bank
*       CFORMM ; Format of the master segment
*       CFORMS ; Format of the data segment
*
*     (Output)
*       IOX    ; Format character word.
*       IERR   ; Error code.
*
*     (Creation Date and Author)
*       1993.10.23 ; First version by J.Kanzaki
*
************************************************************************


#include "zbs.h"
#include "zbssys.h"

      CHARACTER*(*) CNAME

      CHARACTER*(*) CFORMM,CFORMS
      PARAMETER (LFORM=100)
      CHARACTER*100 CFORM
      CHARACTER*1 CBLK/' '/


      IERR = 0

      LFORMM = KZLENG(CFORMM)
      LFORMS = KZLENG(CFORMS)
      IF(LFORMM+LFORMS.GT.LFORM+3) THEN
         IERR = 1
         CALL MSGOUT('KZFORM','Too long format words',1)
         RETURN
      END IF

      IFMAST = 1
      IF(LFORMM.LE.0 .OR. CFORMM(1:LFORMM).EQ.CBLK) IFMAST = 0
      IFSEGM = 1
      IF(LFORMS.LE.0 .OR. CFORMS(1:LFORMS).EQ.CBLK) IFSEGM = 0

      IF(IFMAST.EQ.1) THEN
         IF(IFSEGM.EQ.1) THEN
            CFORM = CFORMM(1:LFORMM)//' / '//CFORMS(1:LFORMS)
         ELSE
            CFORM = CFORMM(1:LFORMM)
         END IF
      ELSE
         IF(IFSEGM.EQ.1) THEN
            CFORM = '/ '//CFORMS(1:LFORMS)
         ELSE
            IERR = 2
            CALL MSGOUT('KZFORM','Invalid Format Specifications',2)
            RETURN
         END IF
      END IF

      CALL MZFORM(CNAME(1:4),CFORM,IOX)

      RETURN
      END
