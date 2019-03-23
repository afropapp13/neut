***********************************************************************
*     -------------------
      SUBROUTINE VCRDFSICM
*     -------------------
*     
*     ( purpose )
*       Clear and read FSI Common
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
      CALL VCCLFSICM
      CALL VCRDFSIV
      CALL VCRDFSIP
      RETURN
      END
