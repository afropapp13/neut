************************************************************************
*
*     ----------------------
      SUBROUTINE KZWRIT(LUN)
*     ----------------------
*
*     (Purpose)
*       Output an event to the file specified by LUN.
*
*     (Input)
*       LUN   ; Logical unit number assigned to the file.
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

      PARAMETER (IOD=2,NUH=10)
      INTEGER*4 IUHEAD(NUH)



      IXEVTD = IXDIV(2)

      CALL KZMNUM(NSBEVT,ISBEVT)
      IF(NSBEVT.LE.0) RETURN

C==== Search for Current Event

      LEVC = LZBS(IPEVC)
      LEMB = LZBS(LEVC-ILEMCR)
      LEHD = LZBS(LEMB-ILEHCR)

      LBIR = LZBS(LEHD-ILBILS)

  100 CONTINUE

      IF(LBIR.LE.0) GO TO 1000

      ITYPE = IZBS(LBIR+IBKTYP)

      IF(ITYPE.EQ.1) THEN

         NSEGM = IZBS(LBIR+INSEGM)
         NSEGML = IZBS(LBIR+INSGMX)
         LSEGM = IZBS(LBIR+ILSEGF)

         IF(NSEGML.GT.0 .AND. NSEGM.GT.NSEGML) THEN
            LPUSH = LSEGM*(NSEGM-NSEGML)
            IZBS(LBIR+INSEGM) = NSEGML
            LDTR=LZBS(LBIR-ILDATR)
            CALL MZPUSH(IXSTOR,LDTR,0,-LPUSH,' ')
C DWC  Must restore local link LBIR in case MZPUSH triggered garbage collection
            LBIR = LZBS(LDTR+1)
         END IF

      ELSE

         NSEGM  = IZBS(LBIR+INSEGM)
         NSEGML = IZBS(LBIR+INSGMX)

         IF(NSEGML.GE.0) THEN
            LPUSH = LINDSG*(NSEGM-NSEGML)
            IF(LPUSH.GT.0) THEN
               IZBS(LBIR+INSEGM) = NSEGML
               CALL MZPUSH(IXSTOR,LBIR,0,-LPUSH,' ')
            END IF
         END IF


         LDTREC = IZBS(LBIR+ILDTRC)
         LFILLD = IZBS(LBIR+ILFILL)
         LPUSH  = LDTREC - LFILLD

         LDTR = LZBS(LBIR-ILDATR)
         IF(LPUSH.GT.0 .AND. LDTR.GT.0) THEN
            IZBS(LBIR+ILDTRC) = LFILLD
            CALL MZPUSH(IXSTOR,LDTR,0,-LPUSH,' ')
         END IF
C  DWC  Must restore local link copy LBIR in case MZPUSH triggered garbage
C  DWC  collection
         LBIR = LZBS(LDTR+1)
      END IF

      LBIR = LZBS(LBIR)

      GO TO 100

 1000 CONTINUE

C DWC   Must restore local link copy LEMB in case anything above triggered
C DWC   garbage collection

      LEVC = LZBS(IPEVC)
      LEMB = LZBS(LEVC-ILEMCR)

      CALL FZOUT(LUN,IXEVTD,LEMB,1,' ',IOD,NUH,IUHEAD)

      RETURN
      END
