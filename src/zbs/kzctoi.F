************************************************************************
*     -----------------
      SUBROUTINE KZCTOI(CNAME,INAME,NCHAR)
*     -----------------
*
*     (Purpose)
*       Move characters into an integer in a well defined order.
*
*     (Input)
*       CNAME     A character variable with NCHAR characters.
*       NCHAR     The number of characters in CNAME
*
*     (Output)
*       INAME     An array of integers*4 with at least NCHAR/4 elements.
*                 This will contain the packed characters.
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

      INTEGER NCHAR

      INTEGER IFROMC

      INTEGER I, J

      DO 110 I=1,NCHAR,4
         J = I + 3
         IF (J.GT.NCHAR) J = NCHAR
         INAME(I/4 + 1) = IFROMC(CNAME(I:J))
 110  CONTINUE
      
      RETURN
      END
