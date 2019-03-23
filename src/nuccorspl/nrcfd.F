      SUBROUTINE NRCFD(YI,N,Y)

      IMPLICIT NONE
      REAL*4 Y(*),YI(*)
      INTEGER*4 N

      INTEGER*4 I
      REAL*4 CF

      CF=0.0
      DO 1 I=1,N
         CF=CF+Y(I)
         YI(I+1)=CF
 1    CONTINUE
      YI(1)=0.0
      cf=1./cf
      DO 2 I=2,n+1
         YI(I)=YI(I)*cf
 2    CONTINUE
      RETURN
      END
