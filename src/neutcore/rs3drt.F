***********************************************************************
*     ---------------------------------------------------
      SUBROUTINE RS3DRT(VECTOR,SINTH,COSTH,SINPHI,COSPHI)
*     ---------------------------------------------------
*
*     ( purpose )
*     rotate 3-Dimension VECTOR with (sinth,costh) and (sinphi,cosphi)
*
*     ( input )
*       vector(3): initial vector
*       sinth    : sin(theta)
*       costh    : cos(theta)
*       sinphi   : sin(phi)
*       cosphi   : cos(phi)
*     ( output )
*       vector(3): final vector
*
*     ( creation date and author )
*       1995.Apr. ; First version by Y.Hayato
*
*     ( comment )
*
***********************************************************************
      
      
      IMPLICIT NONE
      REAL*4   VECTOR(3),SINTH,COSTH,SINPHI,COSPHI

      REAL*4   TMPVEC(3)

      TMPVEC(3) = COSTH*VECTOR(3)  -  SINTH*VECTOR(1)
      TMPVEC(1) = SINTH*VECTOR(3)  +  COSTH*VECTOR(1)
      TMPVEC(2) = VECTOR(2)

      VECTOR(1) = COSPHI*TMPVEC(1) - SINPHI*VECTOR(2)
      VECTOR(2) = SINPHI*TMPVEC(1) + COSPHI*VECTOR(2)
      VECTOR(3) = TMPVEC(3)

      RETURN
      END
      
