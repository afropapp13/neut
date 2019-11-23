************************************************************************
*
*     -----------------------------------------------
      SUBROUTINE KZPUT1(CNAME,ISEGM,NWPOS,NDATA,IDATA)
*     -----------------------------------------------
*
*     (Purpose)
*       Put data words into the given segment of the bank.
*
*     (Input)
*       CNAME ; The name of the bank
*       ISEGM ; Segment#
*       NWPOS ; The position ( no. of words ) where the data
*             ; should be replaced.
*       CFORM ; I/O format word in a character
*       NDATA ; Length of the data
*       IDATA ; Data words to be added
*
*     (Output)
*       none
*
*     (Creation Date and Author)
*       1991.12.25 ; First version by J.Kanzaki
*
************************************************************************

#include "zbs.h"
#include "zbssys.h"

      CHARACTER*(*) CNAME

      INTEGER*4 LBIR,ISEGM,NDATA,IDATA(*)


      CALL KZBLOC(CNAME, LBIR)

      ITYPE = IZBS(LBIR+IBKTYP)
      IF(ITYPE.NE.1) RETURN

C     fixed and 1-dim. case

      NWMAST = IZBS(LBIR+INMAST)
      NSEGM  = IZBS(LBIR+INSEGM)
      LSEGM  = IZBS(LBIR+ILSEGF)
      IF(ISEGM.LT.0) THEN
         CALL MSGOUT('KZPUT1',
     &        'Illegal segment# is specified',1)
         RETURN
      END IF

      IF(ISEGM.GT.0 .AND. NWPOS+NDATA.GT.LSEGM .OR.
     &   ISEGM.EQ.0 .AND. NWPOS+NDATA.GT.NWMAST) THEN
         CALL MSGOUT('KZPUT1',
     &      'DATA length is larger than pre-defined value',2)
         RETURN
      END IF

      IF(ISEGM.GT.NSEGM) CALL KZEXB1(LBIR,ISEGM)

      IF(ISEGM.GT.IZBS(LBIR+INSGMX)) IZBS(LBIR+INSGMX)=ISEGM

      LDTR=LZBS(LBIR-ILDATR)

      IF(ISEGM.EQ.0) THEN
         LDAT = LDTR + NWPOS
      ELSE
         LDAT = LDTR + NWMAST + LSEGM * (ISEGM-1) + NWPOS
      END IF

      DO 100 I=1,NDATA
      IZBS(LDAT+I)=IDATA(I)
  100 CONTINUE

      RETURN
      END
