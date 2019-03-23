      REAL FUNCTION ABSNEPDK(A,B,C)
C=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
C--  CALLED FROM -  GENEV
C     PDK = SQRT(A*A+(B*B-C*C)**2/(A*A) - 2.0*(B*B+C*C))/2.0
      
      IMPLICIT NONE

      REAL*4 A,B,C
      REAL*4 A2,B2,C2
      REAL*4 ARG

      A2 = A*A
      B2 = B*B
      C2 = C*C
      ARG = A2 + (B2-C2)**2/A2 - 2.0*(B2+C2)
      if (ARG.le.0.) then
        ABSNEPDK = 0.0
      else
        ABSNEPDK = 0.5*sqrt(ARG)
      end if
C      ABSNEPDK = 0.5*SQRT(A2 + (B2-C2)**2/A2 - 2.0*(B2+C2))
      RETURN
      END
