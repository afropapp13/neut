************************************************************************
*     ---------------------------------
      SUBROUTINE SORTDL( X,LIST,M0,N,K0)
*     ---------------------------------
C
C     (Purpose)
C        SORT USING THE ALGORITHM 201 OF CACM
C-----   SEE P.66,DESY INFO,74/C
C
C     (Input)
C        
C     (Output)    
C
C     (Creation Date and Author)
C        1978.02.03 ; S.Yamada
C        19??.??.?? ; A. SATO
C
************************************************************************
      REAL*8 X
      DIMENSION  X(M0,N),LIST(N)
C
      IF(N.LE.1) RETURN
      M = N
   10 M = M/2
      IF(M.EQ.0) RETURN
      K = N-M
        DO 50 J=1,K
        I = J
   20   IF(I.LT.1) GO TO 50
        IF(X(K0,LIST(I+M)).GE.X(K0,LIST(I))) GO TO 50
         LW = LIST(I)
        LIST(I) = LIST(I+M)
        LIST(I+M) = LW
        I = I-M
        GO TO 20
   50   CONTINUE
      GO TO 10
      END
