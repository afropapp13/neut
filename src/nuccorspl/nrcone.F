      SUBROUTINE nrcone(SA,SB,CTH,PHI)
      
      IMPLICIT NONE
      REAL*4 CTH,PHI

C     ASSIGNS TO SB(3) A VECTOR FORMING AN	ANGLE ACOS(CTH)	WITH SA(3)
C     IF PHI<0 THEN A RANDOM PHI IS TAKEN (FROM CONICAL SURFACE)
C     
      REAL*4 SA(3),SB(3),D(3,3),RK(3)
      REAL*4 R
      INTEGER*4 I,J

      CALL nrRANVE(1.,RK,CTH,PHI)
      R=SQRT(SA(1)**2+SA(2)**2)

      IF(R.Eq.0.)then
         SB(1)=RK(1)
         sb(2)=rk(2)
         sb(3)=rk(3)*sa(3)
      else
         D(1,1)=SA(1)*SA(3)/R
         D(1,2)=SA(2)*SA(3)/R
         D(1,3)=-R
         D(2,1)=-SA(2)/R
         D(2,2)=SA(1)/R
         D(2,3)=0.
         D(3,1)=SA(1)
         D(3,2)=SA(2)
         D(3,3)=SA(3)
         DO 1 I=1,3
            SB(I)=0.
            DO 1 J=1,3
 1          SB(I)=SB(I)+D(J,I)*RK(J)
       end if

 10   RETURN
      END
