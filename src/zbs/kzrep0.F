************************************************************************
*
*     -----------------------------------------------
      SUBROUTINE KZREP0(CNAME,ISEGM,CFORM,NDATA,IDATA)
*     -----------------------------------------------
*
*     (Purpose)
*       Put data words into the given segment of the bank.
*
*     (Input)
*       CNAME ; The name of the bank
*       ISEGM ; Segment#
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

      IF(LBIR.LE.0) RETURN

      LDATA = NDATA+1
      IOWORD = 16*NDATA+KZIOFM(CFORM)

      ITYPE = IZBS(LBIR+IBKTYP)

      IF(ITYPE.NE.0) RETURN

      NSEGM=IZBS(LBIR+INSEGM)

      IF(ISEGM.LT.0) THEN
         CALL MSGOUT('KZREP0',
     &        'Illegal segment# is specified',3)
         RETURN
      END IF

      IF(ISEGM.GT.NSEGM) CALL KZEXB0(LBIR,ISEGM)

      IF(ISEGM.GT.IZBS(LBIR+INSGMX)) IZBS(LBIR+INSGMX)=ISEGM

      LDTREC=IZBS(LBIR+ILDTRC)
      LFILLD=IZBS(LBIR+ILFILL)

      IF(LFILLD+LDATA.GT.LDTREC) CALL KZEXDR(LBIR,LDATA)

      LDTR=LZBS(LBIR-ILDATR)

      IZBS(LBIR+ILFILL)=LFILLD+LDATA

      ISLAST=IZBS(LBIR+IISLST)

      LIND=LBIR+LINDHD+ISEGM*LINDSG

      LSEGM=IZBS(LIND+1)
      LNSEGM=IZBS(LIND+2)

CCCCC PRINT *,'KZPUTP;ISLAST = ',ISLAST
CCCCC PRINT *,'KZPUTP;LNSEGM = ',LNSEGM

      IF(ISEGM.EQ.ISLAST) THEN

         LDAT=LDTR+LSEGM+LNSEGM
         IZBS(LIND+2)=LNSEGM+LDATA
      ELSE IF(LNSEGM.GT.0) THEN

         ILAST=LDTR+IZBS(LBIR+LINDHD+ISLAST*LINDSG+1)
     &             +IZBS(LBIR+LINDHD+ISLAST*LINDSG+2)
         ISTART=LDTR+LSEGM+LNSEGM
         LSHIFT=ILAST-ISTART
         DO 100 I=LSHIFT,1,-1
         IZBS(ISTART+LDATA+I) = IZBS(ISTART+I)
  100    CONTINUE
         LDAT=ISTART
         IZBS(LIND+2)=LNSEGM+LDATA
         DO 110 I=0,NSEGM
         IPS=IZBS(LBIR+LINDHD+I*LINDSG+1)
         ILS=IZBS(LBIR+LINDHD+I*LINDSG+2)
         IF(IPS.GT.LSEGM .AND. ILS.NE.0)
     &      IZBS(LBIR+LINDHD+I*LINDSG+1)=IPS+LDATA
  110    CONTINUE
      ELSE

         LSEGM=0
         IF(ISLAST.NE.-1)
     &   LSEGM=IZBS(LBIR+LINDHD+ISLAST*LINDSG+1)
     &        +IZBS(LBIR+LINDHD+ISLAST*LINDSG+2)
         LDAT=LDTR+LSEGM
         IZBS(LBIR+IISLST)=ISEGM
         IZBS(LIND+1)=LSEGM
         IZBS(LIND+2)=LDATA
      END IF

      IZBS(LDAT+1) = IOWORD
      DO 200 I=1,NDATA
      IZBS(LDAT+I+1)=IDATA(I)
  200 CONTINUE

      RETURN
      END
