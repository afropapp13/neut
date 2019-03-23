***********************************************************************
*     -------------------
      SUBROUTINE VCMKFSICM
*     -------------------
*     
*     ( purpose )
*       make FSI Common
*
*     ( input ) 
*       none
*
*     ( output ) 
*       COMMON /VERTWORK/
*
*     ( creation date and author )
*       2010.09.19 ; First version by Y.Hayato
*
***********************************************************************
      IMPLICIT NONE
      CALL VCMKFSIV
      CALL VCMKFSIP
      RETURN
      END
