c
c    $Id: kzreda.F,v 1.1 1997/10/08 17:23:47 superk Exp $
c    $Log: kzreda.F,v $
c    Revision 1.1  1997/10/08 17:23:47  superk
c    KZREAD with A option (KS)
c
c
************************************************************************
*     ----------------------------
      SUBROUTINE KZREDA (LUN,IERR)
*     ----------------------------
*
*     (Purpose)
*       Read one event from sequential file and restore it in the
*       Venus Data Structure.  MODIFIED VERSION with "A" option
*       to FZIN, for use with zootsuit (K. Scholberg)
*
*     (Input)
*       LUN    : Logical unit# of the input device.
*
*     (Output)
*       IERR =  1 : Find EOF while reading data.
*       IERR =  2 : Ilegal sequence of record occured while reading.
*
*     (Comment)
*       For the format of data in a seqential file, see the comment
*       in SUBROUTINE ZVWRIT.
*
*     (History)
*       1991.12.25 : First version by J.Kanzaki
*       1996.04.25 :  J.Kanzaki, Add declarations.
*       1997.10.08 : K. Scholberg, make into KZREADA
*
************************************************************************

#include "zbsdeclr.h"
#include "zbssys.h"
#include "zbs.h"

      INTEGER LUN,IERR

      INTEGER NUH0
      PARAMETER (NUH0=100)
      INTEGER IUHEAD(NUH0)

      INTEGER IXEVTD
      INTEGER LEVC,LEMB,LLAST

      INTEGER NUH

      INTEGER LZLAST


***** Initial set up.
      IERR   = 0

      IXEVTD = IXDIV(2)

***** Get pointer to the event control bank
      LEVC = LZBS(IPEVC)
      LLAST = LZLAST(IXSTOR,LZBS(LEVC-ILEMLS))

      NUH = NUH0

      IF(LLAST.NE.0) THEN

         CALL FZIN(LUN,IXEVTD,LLAST,0,'A',NUH,IUHEAD)
         LEMB = LZBS(LLAST)

      ELSE

         CALL FZIN(LUN,IXEVTD,LEVC,-1,'A',NUH,IUHEAD)
         LEMB = LZBS(LEVC-ILEMLS)

      END IF

      IF(IQUEST(1).GE.2) THEN
         IERR = 1
         RETURN
      END IF

      LZBS(LEVC-ILEMCR) = LEMB

C     CALL KZCLDT

      CALL KZMKDT

      RETURN
      END
