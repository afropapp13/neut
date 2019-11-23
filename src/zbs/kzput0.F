************************************************************************
*
*     -----------------------------------------------------
      SUBROUTINE KZPUT0(CNAME,ISEGM,NWPOS,CFORM,NDATA,IDATA)
*     -----------------------------------------------------
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
      CHARACTER*1 CFORM

      INTEGER*4 LDATA,IOWORD


      CALL KZBLOC(CNAME, LBIR)

      LDATA = NDATA+1
      IOWORD = 16*NDATA+KZIOFM(CFORM)

      ITYPE = IZBS(LBIR+IBKTYP)

      IF(ITYPE.NE.0) RETURN

      NSEGM=IZBS(LBIR+INSEGM)

      IF(ISEGM.LT.0) THEN
         CALL MSGOUT('KZPUT0',
     &        'Illegal segment# is specified',1)
         RETURN
      END IF

      IF(ISEGM.GT.NSEGM) CALL KZEXB0(LBIR,ISEGM)

      IF(ISEGM.GT.IZBS(LBIR+INSGMX)) IZBS(LBIR+INSGMX)=ISEGM

      LDTR=LZBS(LBIR-ILDATR)

      LDTREC=IZBS(LBIR+ILDTRC)
      LFILLD=IZBS(LBIR+ILFILL)

      IF(LFILLD+LDATA.GT.LDTREC) CALL KZEXDR(LBIR,LDATA)

      LDTR=LZBS(LBIR-ILDATR)

      CALL KZWPOS(LBIR,ISEGM,NWPOS, LPOS)

      IF(LPOS.LT.0) THEN
         CALL MSGOUT('KZPUT0',
     &        'Insert position location error',2)
         RETURN
      END IF

      DO 100 I=1,NDATA
      IZBS(LDTR+LPOS+I)=IDATA(I)
  100 CONTINUE

      RETURN
      END
