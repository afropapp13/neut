      SUBROUTINE nrLOREN(C,P,Q,N,K)
      implicit none

      integer*4 n,k
      real*4 C(4),P(4,N),Q(4,N)
      real*4 S,W,Z,D,E,A
      integer*4 I,J

      S=0.
      DO 1 I=1,3
 1    S=S+C(I)**2
      IF(S.EQ.0.)GO TO 100
      S=C(4)**2-S
      W=SQRT(S)
      Z=W+C(4)
      DO 20 J=1,N
         D=C(4)*P(4,J)
         DO 21 I=1,3
 21      D=D-K*C(I)*P(I,J)
         E=D/W
         A=(P(4,J)+E)/Z
         Q(4,J)=E
         DO 22 I=1,3
 22      Q(I,J)=P(I,J)-A*K*C(I)
 20   CONTINUE
      GO TO 1000
 100  DO 101 J=1,N
         DO 101 I=1,4
 101  Q(I,J)=P(I,J)
      GO TO 1000
 1000 RETURN
      END
