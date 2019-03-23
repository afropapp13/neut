      REAL FUNCTION ABSNEPDK(A,B,C)
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C--  CALLED FROM -  GENEV
C     PDK = SQRT(A*A+(B*B-C*C)**2/(A*A) - 2.0*(B*B+C*C))/2.0
      
      IMPLICIT NONE

      REAL*4 A,B,C
      REAL*4 A2,B2,C2

      A2 = A*A
      B2 = B*B
      C2 = C*C
      ABSNEPDK = 0.5*SQRT(A2 + (B2-C2)**2/A2 - 2.0*(B2+C2))
      RETURN
      END
