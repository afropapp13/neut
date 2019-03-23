      SUBROUTINE nr2iso(AMD,AM1,AM2,PD,P1)
      implicit none
      real*4 PD(4),P1(4,2),PTEM(4,2),dir(3)
      real*4 amd,am1,am2,amd2,a2m1,a2m2,e1,e2,p
      real*4 cs,ss,pdd
      integer i

      AMD2=AMD**2
      A2M1=AM1**2
      A2M2=AM2**2
      E1=(AMD2+A2M1-A2M2)/AMD/2.
      E2=(AMD2+A2M2-A2M1)/AMD/2.
      P=E1**2-A2M1
      P=SQRT(P)
      call nrazio(cs,ss)
      pdd=sqrt(pd(1)**2+pd(2)**2+pd(3)**2)
      do i=1,3
         dir(i)=pd(i)/pdd
      end do
      CALL nrcone(dir,ptem,cs,-1.)
      do i=1,3
         ptem(i,1)=ptem(i,1)*p
         ptem(i,2)=-ptem(i,1)
      end do
      PTEM(4,1)=E1
      PTEM(4,2)=E2
      CALL nrLOREN(PD,PTEM,P1,2,-1)
      RETURN
      END
