***********************************************************************
*     -------------------
      SUBROUTINE VCRDVCCM
*     -------------------
*     
*     ( purpose )
*       Read MCVECT & MCVERTEX bank
*
*     ( input ) 
*       none
*
*     ( output ) 
*       COMMON /VCWORK/
*       COMMON /VCVRTX/
*
*     ( creation date and author )
*       1995.03.16 ; First version by K.Kaneyuki
*
***********************************************************************
      CALL VCCLCM
      CALL VCRDVC
      CALL VCRDVX
      RETURN
      END
