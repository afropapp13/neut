************************************************************************
*     ------------------------------------------
      SUBROUTINE KZCHG0(CNAME,ISEGM,NWPOS,IDATA)
*     -------------------------------------------
*
*     (Purpose)
*       Exchange one word of data in the given segment of the bank.
*
*     (Input)
*       CNAME ; The name of the bank
*       ISEGM ; Segment#
*       NWPOS ; The position where the data should be replaced.
*       IDATA ; Data words to be added
*
*     (Output)
*       none
*
*     (Creation Date and Author)
*       1996.05.22 ; First version by J.Kanzaki
*
************************************************************************

#include "zbs.h"
#include "zbssys.h"

      CHARACTER*(*) CNAME

      INTEGER LBIR,ISEGM,IDATA


      CALL KZBLOC(CNAME, LBIR)

      LDATA = NDATA+1
      IOWORD = 16*NDATA+KZIOFM(CFORM)

      ITYPE = IZBS(LBIR+IBKTYP)

      IF(ITYPE.NE.0) RETURN

      NSEGM=IZBS(LBIR+INSEGM)

      IF(ISEGM.LT.0 .OR. 
     &   ISEGM.GT.NSEGM .OR.
     &   ISEGM.GT.IZBS(LBIR+INSGMX)) THEN
         CALL MSGOUT('KZCHG0',
     &        'Illegal segment# is specified',1)
         RETURN
      END IF

      LDTR=LZBS(LBIR-ILDATR)

      LDTREC=IZBS(LBIR+ILDTRC)
      LFILLD=IZBS(LBIR+ILFILL)

      IF(LFILLD+LDATA.GT.LDTREC) CALL KZEXDR(LBIR,LDATA)

      LDTR=LZBS(LBIR-ILDATR)

      CALL KZWPOS(LBIR,ISEGM,NWPOS, LPOS)

      IF(LPOS.LT.0) THEN
         CALL MSGOUT('KZCHG0',
     &        'Insert position location error',2)
         RETURN
      END IF

      IZBS(LDTR+LPOS)=IDATA

      RETURN
      END
