***********************************************************************
*     -------------------
      FUNCTION RLU(IDUMMY)
*     -------------------
*     
*     ( purpose )
*       random number
*
*     ( input ) 
*       IDUMMY : DUMMY INPUT
*
*     ( output ) 
*       RLU : (0,1) rundum number using RANLUX
*
*     ( creation date and author )
*       2001. 2. 9 ; First version by Y.Hayato
*
*     ( comment )
*
***********************************************************************

      IMPLICIT NONE

      REAL*4    RLU
      INTEGER*4 IDUMMY

      REAL*4    RNUM

      CALL RANLUX(RNUM,1)

      RLU = RNUM
      RETURN
      END
