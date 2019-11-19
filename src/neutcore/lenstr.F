************************************************************************
*     ---------------------------------
      INTEGER*4 FUNCTION LENSTR(S)
*     ---------------------------------
*
*     (Purpose)
*       Return length of the string S
*
*     (Input)
*       S     : input STRING
*
*     (Output)
*       LENSTR: Actual length of the string
*
************************************************************************
      implicit none

      character*(*) S
      integer*4 i

      DO 100 I=1,LEN(S)
         if (S(I:I).eq.' ') then
            goto 110
         endif
 100  continue

 110  continue
      lenstr  = i - 1
      return
      end
