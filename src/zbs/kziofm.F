************************************************************************
*
*     ------------------------------
      INTEGER FUNCTION KZIOFM(CFORM)
*     ------------------------------
*
*     (Purpose)
*       Create I/O format word from a given format in a character form.
*
*     (Input)
*       CFORM ; Format in a character string.
*
*     (Output)
*       KZIOFM ; I/O format word
*
*     (Creation Date and Author)
*       1991.12.25 ; First version by J.Kanzaki
*
************************************************************************


      CHARACTER*1 CFORM


      IF(CFORM.EQ.'B') THEN

         KZIOFM = 1

      ELSE IF(CFORM.EQ.'I') THEN

         KZIOFM = 2

      ELSE IF(CFORM.EQ.'F') THEN

         KZIOFM = 3

      ELSE IF(CFORM.EQ.'D') THEN

         KZIOFM = 4

      ELSE IF(CFORM.EQ.'H') THEN

         KZIOFM = 5

      ELSE

         KZIOFM = 0

      END IF

      RETURN
      END
