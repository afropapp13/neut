      real FUNCTION     nrPRBIN(YI,N,XMIN,XSTEP,YR)

      implicit none
      real*4 YI(*),XMIN,XSTEP,YR,fr

      integer*4 N,MIN,max,m
      MIN=1
      MAX=N+1
      GO TO 3
 1    IF((MAX-MIN).EQ.1) GO TO  5
      IF(YI(M)-YR) 2,6,4
 2    MIN=M
 3    M=(MAX+MIN)/2
      GO TO 1
 4    MAX=M
      GO TO 3
 5    M=MIN
      FR=(YR-YI(MIN))/(YI(MAX)-YI(MIN))
      nrprbin=XMIN+XSTEP*(FLOAT(MIN)+FR-1.0)
      RETURN
 6    FR=0.0
      nrprbin=XMIN+XSTEP*(FLOAT(M)-1.0)
      RETURN
      END
