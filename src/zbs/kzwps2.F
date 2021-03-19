************************************************************************
*
*     -----------------------------------------------------
      SUBROUTINE KZWPS2(CNAME,ISEGM,ISSEG,NWPOS,NDATA, LDAT)
*     -----------------------------------------------------
*
*     (Purpose)
*       Put data words into the given segment of the bank.
*
*     (Input)
*       CNAME ; The name of the bank
*       ISEGM ; Segment#
*       ISSEG ; Sub-segment#
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

      INTEGER*4 LBIR,ISEGM,ISSEG,NDATA


      CALL KZBLOC(CNAME, LBIR)

      ITYPE = IZBS(LBIR+IBKTYP)
      IF(ITYPE.NE.2) RETURN

C     fixed and 2-dim. case

      NWMAST = IZBS(LBIR+INMAST)
      NSEGM  = IZBS(LBIR+INSEGM)
      LSEGM  = IZBS(LBIR+ILSEGF)

      LDTR=LZBS(LBIR-ILDATR)

      LDTREC=IZBS(LBIR+ILDTRC)
      LFILLD=IZBS(LBIR+ILFILL)

      ISLAST=IZBS(LBIR+IISLST)

      LIND=LBIR+LINDHD+ISEGM*LINDSG

      IPSEG=IZBS(LIND+1)
      NSSEG=IZBS(LIND+2)

      IF(ISEGM.EQ.0) THEN

         LDAT = LDTR + NWPOS

      ELSE IF(IPSEG.EQ.-1) THEN

         LADD = LSEGM*ISSEG

         IF(LFILLD+LADD.GT.LDTREC)
     &      CALL KZEXDR(LBIR,LADD)

         CALL KZBLOC(CNAME, LBIR)

         LDTR=LZBS(LBIR-ILDATR)

         IZBS(LIND+1) = LFILLD
         IZBS(LIND+2) = ISSEG

         LDAT = LDTR+LFILLD+LSEGM*(ISSEG-1)+NWPOS

         IZBS(LBIR+ILFILL) = LFILLD+LADD

         IZBS(LBIR+IISLST) = ISEGM

      ELSE

         IF(ISSEG.EQ.0 .OR. ISSEG.GT.NSSEG) THEN

            IF(ISSEG.EQ.0) THEN

               LADD = LSEGM

            ELSE

               LADD = LSEGM*(ISSEG-NSSEG)

            END IF

            IF(LFILLD+LADD.GT.LDTREC)
     &         CALL KZEXDR(LBIR,LADD)

            CALL KZBLOC(CNAME, LBIR)

            LDTR=LZBS(LBIR-ILDATR)

            IZBS(LBIR+ILFILL)=LFILLD+LADD

            ILAST=LDTR+IZBS(LBIR+LINDHD+ISLAST*LINDSG+1)
     &          + LSEGM*IZBS(LBIR+LINDHD+ISLAST*LINDSG+2)
            ISTART=LDTR+IPSEG+LSEGM*(ISSEG-1)
            LSHIFT=ILAST-ISTART
            DO 100 I=LSHIFT,1,-1
            IZBS(ISTART+LSEGM+I)=IZBS(ISTART+I)
  100       CONTINUE
            LDAT=ISTART+NWPOS
            IZBS(LIND+2)=ISSEG
            DO 110 I=1,NSEGM
            IPS=IZBS(LBIR+LINDHD+I*LINDSG+1)
            IF(IPS.GT.IPSEG) THEN
               IZBS(LBIR+LINDHD+I*LINDSG+1)=IPS+LADD
            END IF
  110       CONTINUE

         ELSE

            LDAT = LDTR+IPSEG+LSEGM*(ISSEG-1)+NWPOS

         END IF

      END IF

      RETURN
      END
