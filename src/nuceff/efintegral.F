      program efintegral

C  Integral Routine (CERN lib)
      real*8 gpindp
      external gpindp
      real*8 dinteg
      external dinteg
      real*8 epsin,epsout
      real*8 CC
      real*8 result

      CC=2.27*2.5

      iop=1
      epsin=1.D-7
      result = gpindp(0.,CC,epsin,epsout,dinteg,iop)

      write(*,*) "result= ",result*4*3.141592
      write(*,*) "4piR**3/3=",(4*3.141592*2.27*2.27*2.27)/3
      end

C
      FUNCTION dinteg(R)
C
C ++ FUNCTION TO GIVE NUCLEAR DENSITY (UNIT : AT CENTER OF NULEI -- 1)
C
      real*8 dinteg
      real*8 R
      C=2.27
      AF=A*0.2275598
      X=(-C)/AF
      DFACT=1./(1.+EXP(X))
      X=(R-C)/AF
C      dinteg=R*R/(1.+EXP(X))/DFACT
      dinteg=R*R/(1.+EXP(X))
C      dinteg=R*R
      RETURN
      END
