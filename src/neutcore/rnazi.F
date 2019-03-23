************************************************************************
*     ---------------------
      SUBROUTINE RNAZI(C,S)
*     ---------------------
*
*     (Purpose)
*       Give cosine and sine of random direction
*     
*     (Input)
*       NONE
*
*     (Output)
*       C      : COSINE OF RANDOM DIRECTION
*       S      : SINE OF RANDOM DIRECTION      
*
*     (Creation Date and Author)
*       1978.09.08 ; S.Yamada, A.Sato
*       1995.02.03 ; K. KANEYUKI  FOR S.K.
*                    RANAZI -> RNAZI
*     
************************************************************************

    1 XR1=RLU(DUM)
      XR2=RLU(DUM)
      IF(XR1.EQ.0. .AND. XR2.EQ.0.) GO TO 90
      XR1=2.*(XR1-.5)
      X=XR1*XR1+XR2*XR2
      IF ( X.GT.1.0 .OR.  X.LE.0.0) GO TO 1
C      IF(X-1.) 2,2,1
C    2 IF(X) 1,1,3
C    3 C=(XR1*XR1-XR2*XR2)/X
      C=(XR1*XR1-XR2*XR2)/X      
      S=2.*XR1*XR2/X
      RETURN
C     
   90 WRITE(6,690)
  690 FORMAT('0**** 2 SUCCESSIVE RANDOM NUMBERS ARE 0.IN NEWDR1.***',
     1     /'  EMERGENCY STOP')
      STOP
      END
