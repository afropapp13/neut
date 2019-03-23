************************************************************************
*     ---------------------------------
      REAL FUNCTION RSY00(THETA,PHI)
*     ---------------------------------
C
C     (Purpose)
C        Calculate spherical harmonics function Y00
C
C     (Input)
C        Theta : theta      
C        phi   : phi
C     (Output)    
C        Y00
C
C     (Creation Date and Author)
C        199?.??.?? ; Y.Hayato      
C
************************************************************************      
c      REAL*4 RSY00
      REAL*4 THETA,PHI,PI
      PI = 3.14159
      RSY00 = SQRT(1. /( 4. * PI ))
      RETURN
      END
