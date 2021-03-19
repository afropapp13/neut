************************************************************************
*     ------------------------------
      FUNCTION KZCMPN(LNAME,K,INAME)
*     ------------------------------
*
*     (Purpose)
*        Compare two bank names; one at IZBS(LNAME(i)+K) and
*        the other INAME(i).
*
*     (Inputs)
*        LNAME  ; Pointers to the bank name tables.
*        K      ; Position of the name to be compared.
*        INAME  ; Array for the reference name.
*
*     (Output)
*        KCCMPN     ; Returns the order of the names.
*               < 0   LNAME(i)+K is before INAME(i).
*               = 0       "      is the same as   "  .
*               > 0       "      is after     "   .
*
*     (Creation Date)
*        1991.12.25 : First version by J.Kanzaki.
*
************************************************************************

      INTEGER*4 INAME(8),LNAME(8)

#include "zbssys.h"
#include "zbs.h"


      DO 10 I=1,8
      IF(IZBS(LNAME(I)+K).GT.INAME(I)) THEN
         KZCMPN=-1
         RETURN
      ELSE IF(IZBS(LNAME(I)+K).LT.INAME(I)) THEN
         KZCMPN=+1
         RETURN
      END IF
   10 CONTINUE

      KZCMPN = 0
      RETURN
      END
