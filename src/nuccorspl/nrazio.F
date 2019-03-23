      subroutine nrAZIO(COSINE,SINE)
      implicit none

      REAL*4 COSINE,SINE

      REAL*4 R1,R1SQ
      REAL*4 R2,R2SQ
      REAL*4 SUM

      real*4   ranf
      external ranf

    5 r1=ranf()
      R1SQ=R1*R1
C     XSQ
      r2=ranf()
      R2SQ=R2*R2
C     YSQ
      SUM=R1SQ+R2SQ
      IF(SUM.GT.1.0)GOTO5
      SUM=SUM*0.5
C     (XSQ+YSQ)/2
      COSINE=(SUM-R1SQ)/SUM
C     (YSQ-XSQ)/(XSQ+YSQ)
      SINE=(R1*R2)/SUM
C     (2*X*Y)/(XSQ+YSQ)
      r1=ranf()
      IF(R1.LT.0.5)GOTO10
      SINE=-SINE
   10 RETURN
      END
