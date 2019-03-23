************************************************************************
*     ---------------------------------
      REAL FUNCTION RSY21(THETA,PHI)
*     ---------------------------------
C
C     (Purpose)
C        Calculate spherical harmonics function Y21
C
C     (Input)
C        Theta : theta      
C        phi   : phi
C     (Output)    
C        Y21
C
C     (Creation Date and Author)
C        199?.??.?? ; Y.Hayato      
C
************************************************************************      

c      REAL*4 RSY21
      REAL*4 THETA,PHI,PI
      PI = 3.14159
      RSY21 =
     $  -SQRT(15. /( 8. * PI )) * COS(THETA)*SIN(THETA) * COS(PHI)
      RETURN
      END
