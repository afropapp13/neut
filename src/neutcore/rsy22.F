************************************************************************
*     ---------------------------------
      REAL FUNCTION RSY22(THETA,PHI)
*     ---------------------------------
C
C     (Purpose)
C        Calculate spherical harmonics function Y22
C
C     (Input)
C        Theta : theta      
C        phi   : phi
C
C     (Output)    
C        Y22
C
C     (Creation Date and Author)
C        199?.??.?? ; Y.Hayato      
C
************************************************************************      

c      REAL*4 RSY22
      REAL*4 THETA,PHI,PI
      PI = 3.14159
      RSY22 = SQRT(15. /( 32. * PI )) * SIN(THETA)*SIN(THETA)
     $     * COS(2 * PHI)
      RETURN
      END
      
