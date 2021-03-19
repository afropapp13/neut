************************************************************************
*
*     -------------------
      FUNCTION KZLENG(CSTR)
*     -------------------
*
*     (Purpose)
*       Returns a length of the given character string.
*       It does not count blanks.
*
*     (Input)
*       CSTR ; Character string
*
*     (Output)
*       KZLENG ; Length of the string
*
*     (Creation Date and Author)
*       1991.12.25 ; First version by J.Kanzaki
*
************************************************************************


      CHARACTER*(*) CSTR
      CHARACTER*1 CBLK/' '/

      LSTR=LEN(CSTR)

      IF(LSTR.GE.0) THEN
         DO 10 I=LSTR,1,-1
            IF(CSTR(I:I).NE.CBLK) THEN
               KZLENG=I
               RETURN
            END IF
 10      CONTINUE
      END IF

      KZLENG=0

      RETURN
      END
