************************************************************************
*
*     -----------------------------------------------
      SUBROUTINE KZREP2(CNAME,ISEGM,ISSEG,NDATA,IDATA)
*     -----------------------------------------------
*
*     (Purpose)
*       Put data words into the given segment of the bank.
*
*     (Input)
*       CNAME ; The name of the bank
*       ISEGM ; Segment#
*       ISSEG ; Sub-segment#
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

      INTEGER*4 LDATA


      CALL KZBLOC(CNAME, LBIR)

      IF(LBIR.LE.0) RETURN

      LDATA = NDATA

      ITYPE = IZBS(LBIR+IBKTYP)

      IF(ITYPE.NE.2) RETURN

      NWMAST = IZBS(LBIR+INMAST)
      NSEGM  = IZBS(LBIR+INSEGM)
      LSEGM  = IZBS(LBIR+ILSEGF)
      IF(ISEGM.LT.0) THEN
         CALL MSGOUT('KZREP2',
     &        'Illegal segment# is specified',1)
         RETURN
      END IF

      IF(ISEGM.GT.0 .AND. NDATA.NE.LSEGM .OR.
     &   ISEGM.EQ.0 .AND. NDATA.NE.NWMAST) THEN
         CALL MSGOUT('KZREP2',
     &      'DATA length is larger than pre-defined value',2)
         RETURN
      END IF

      IF(ISEGM.GT.NSEGM) CALL KZEXB2(LBIR,ISEGM)

      IF(ISEGM.GT.IZBS(LBIR+INSGMX)) IZBS(LBIR+INSGMX)=ISEGM

      CALL KZWPS2(CNAME,ISEGM,ISSEG,0,LSEGM, LDAT)

      DO 200 I=1,NDATA
      IZBS(LDAT+I)=IDATA(I)
  200 CONTINUE

      RETURN
      END
