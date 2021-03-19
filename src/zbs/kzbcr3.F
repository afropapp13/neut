************************************************************************
*     $Id: kzbcr3.F,v 1.1 1997/10/27 18:58:31 mcgrew Exp $
*     -----------------------------------------
      SUBROUTINE KZBCR3(CNAME,NLINKS,NWMAST,CFORMM,
     $                  NSEGM,LSEGM,CFORMS,IERR) 
*     -----------------------------------------
*
*     (Purpose)
*       Create new ZBS bank with down links, a fixed section, followed by
*       a repeated section.  This creates a bank just like KZBCR1, except
*       This bank can support a down structure.
*
*     (Input)
*        CNAME ; Name of the new bank ( <= 32 characters ).
*        NLINKS ; The number of links to create.
*        NWMAST ; The length of the fixed section.
*        CFORMM ; IO format string of fixed section.
*        NSEGM  ; The length of the variable section
*        LSEGM  ; Number of repeats in variable section.
*        CFORMS ; IO format string of variable section.
*
*     (Output)
*       IERR     ; Error flag
*            = 0 ; Normal return
*            > 0 ; Error return
*
*     (Creation Date and Author)
*       1997.10.21 ; First version by C.McGrew
*
************************************************************************


      CHARACTER*(*) CNAME
      CHARACTER*32 CNAME1

      CHARACTER*(*) CFORMM,CFORMS

      INTEGER NLINKS

#include "zbssys.h"
#include "zbs.h"


      CALL KZCRMH

      CNAME1 = CNAME

      CALL KZPDIR(CNAME,IPDIR)

      LDIR = LZBS(IPDIRT)
      LEXF=LZBS(LDIR-ILBEXF)

      IEXIF=IZBS(LEXF+IPDIR)

      IF(IEXIF.NE.0) THEN
         IERR=1
         RETURN
      END IF

C     STORE INFORMATIONS ON THE NEW BANK INTO THE DIRECTORY

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

      CALL KZCRB3(IBANK,CNAME1,NLINKS,NWMAST,NSEGM,LSEGM)

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

      CALL KZCRD3(IBANK,CNAME1,NLINKS,NWMAST,NSEGM,LSEGM,IOX)

CCCCC IPBANK = LZBS(IBANK)

      RETURN
      END
