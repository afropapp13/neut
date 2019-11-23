************************************************************************
*     -----------------
      SUBROUTINE KZITOC(INAME,CNAME,NWORD)
*     -----------------
*
*     (Purpose)
*       Move characters from an integer into a character variable.
*
*     (Input)
*       INAME     An integer array with NWORD elements.
*       NWORD     The number of elements in INAME.
*
*     (Output)
*       CNAME     An character variable with 4*NWORD characters.
*                 This will contain the unpacked characters.
*
*     (Comment)
*       This allows ZBS to work on both big and little endian machines.
*
*     (Creation Date and Author)
*       1997.08.21 ; First version by C.McGrew
*
************************************************************************

      IMPLICIT NONE

      CHARACTER*(*) CNAME

      INTEGER INAME(*)

      INTEGER NWORD

      CHARACTER*4 CFROMI

      INTEGER I

      DO 110 I=1,NWORD
         CNAME(4*(I-1)+1:4*I) = CFROMI(INAME(I))
 110  CONTINUE

      RETURN
      END
