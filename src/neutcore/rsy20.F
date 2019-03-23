************************************************************************
*     ---------------------------------
      REAL FUNCTION RSY20(THETA,PHI)
*     ---------------------------------
C
C     (Purpose)
C        Calculate spherical harmonics function Y20
C
C     (Input)
C        Theta : theta      
C        phi   : phi
C     (Output)    
C        Y20
C
C     (Creation Date and Author)
C        199?.??.?? ; Y.Hayato      
C
************************************************************************      

c      REAL*4 RSY20
      REAL*4 THETA,PHI,PI
      PI = 3.14159
      RSY20 = SQRT(5. /( 16. * PI )) * (2 * (COS(THETA)*COS(THETA))
     $     - (SIN(THETA)*SIN(THETA)))
      RETURN
      END
