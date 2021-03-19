************************************************************************
*
*     ------------------------------------------
      SUBROUTINE KZBCR2(CNAME,NWMAST,CFORMM,
     &                  NSEGM,LSEGM,CFORMS,IERR) 
*     ------------------------------------------
*
*     (Purpose)
*       Create new ZBS bank with the fixed length segment
*
*     (Input)
*       CNAME ; Name of the new bank ( <= 32 characters ).
*
*     (Output)
*       IERR     ; Error flag
*            = 0 ; Normal return
*            > 0 ; Error return
*
*     (Creation Date and Author)
*       1991.12.25 ; First version by J.Kanzaki
*
************************************************************************


      CHARACTER*(*) CNAME
      CHARACTER*32 CNAME1

      CHARACTER*(*) CFORMM,CFORMS

#include "zbssys.h"
#include "zbs.h"


      CALL KZCRMH

C==== SEARCH FOR BANK "CNAME" FROM THE DIRECTORY

      CNAME1 = CNAME

      CALL KZPDIR(CNAME1,IPDIR)

      LDIR = LZBS(IPDIRT)
      LEXF=LZBS(LDIR-ILBEXF)

      IEXIF=IZBS(LEXF+IPDIR)

      IF(IEXIF.NE.0) THEN
         IERR=1
         RETURN
      END IF

C     STORE INFORMATIONS ON THE NEW BANK INTO THE DIRECTORY

      LEXF=LZBS(LDIR-ILBEXF)
      LSRT=LZBS(LDIR-ILBSRT)

      LPNT=LZBS(LDIR-ILBPNT)
      LIOF=LZBS(LDIR-ILBIOF)
      LIOX=LZBS(LDIR-ILBIOX)
      LEVH=LZBS(LDIR-ILBEVH)

      IERR=0

      IZBS(LEXF+IPDIR)=1
      IZBS(LIOF+IPDIR)=1

C==== CREATE BANK RECORDS

      IBANK=IZBS(LPNT+IPDIR)

      CALL KZCRB2(IBANK,CNAME1,NWMAST,NSEGM,LSEGM)

      IOX = IZBS(LIOX+IPDIR)
      IEVH = IZBS(LEVH+IPDIR)

      CALL KZGETH(IEVH0)

      IF(IOX.EQ.0 .OR. IEVH0.NE.IEVH) THEN

         CALL KZFORM(CNAME1,CFORMM,CFORMS,IOX,IERIOX)

         IF(IERIOX.NE.0) THEN
            IERR = IERIOX
            RETURN
         END IF

         IZBS(LIOX+IPDIR) = IOX
         IZBS(LEVH+IPDIR) = IEVH0

      END IF

      CALL KZCRD2(IBANK,CNAME1,NWMAST,NSEGM,LSEGM,IOX)

CCCCC IPBANK = LZBS(IBANK)

      RETURN
      END
