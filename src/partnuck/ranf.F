***********************************************************************
*     -------------------
      FUNCTION RANF()
*     -------------------
*     
*     ( purpose )
*       random number for kaon tracking
*
*     ( input ) 
*       IDUMMY : DUMMY INPUT
*
*     ( output ) 
*       RANF : (0,1) rundum number using RLU
*
*     ( creation date and author )
*       2001. 2. 9 ; First version by Y.Hayato
*
*     ( comment )
*
***********************************************************************

      IMPLICIT NONE

      REAL*4    RANF

      REAL*4    RLU
      external  RLU

      INTEGER*4 IDUMMY

      RANF = RLU(IDUMMY)

      RETURN
      END
