C
      real*4 an,zz,c,a,cnn,ann,dr,af,cc2,dfact
      real*4 rmsrad,WPARM
      common /eftarget/an,zz,c,cnn,af,cc2,dfact
     &       ,rmsrad,WPARM
      common /efpars/a,ann,dr

C     Parameters listed below are defined in nesettarg.F

C      parameter(AN=16.)
C      parameter(ZZ=8.)
C      parameter(C=2.69)
C      parameter(A=1.80)
C      parameter(CNN=2.69)
C      parameter(ANN=1.80)
C      parameter(DR=0.05)
C      parameter(AF=0.40961)
C      parameter(cc2=2.5*C)
C      parameter(dfact=0.9985962)
C      parameter(dfact=1./(1.+EXP(-C/AF)))
