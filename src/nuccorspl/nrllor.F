      SUBROUTINE NRLLOR(PPRIME,EPRIME,P,E,BETA)

      IMPLICIT NONE
      REAL*4    E,EPRIME
      REAL*4    PPRIME(3),P(3),BETA(3)

      REAL*4    BSQ,BDP,GAM,SHIFT
      INTEGER*4 I

      BSQ=0.
      BDP=0.

      DO 10 I=1,3
         BSQ=BSQ+BETA(I)**2
 10   BDP=BDP+BETA(I)*P(I)

      IF(BSQ) 30,30,15

 15   GAM=SQRT(1./(1.-BSQ))
      SHIFT=(GAM-1.)*BDP/BSQ -GAM*E
      DO 20 I=1,3
 20   PPRIME(I)=P(I)+SHIFT*BETA(I)

      EPRIME=GAM*(E-BDP)
      RETURN

 30   DO 40 I=1,3
 40      PPRIME(I)=P(I)
      EPRIME=E
      RETURN
      END
