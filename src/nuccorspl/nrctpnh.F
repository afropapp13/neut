      subroutine nrctpnh(dchn,dchnf,dchnb,ped,cst,ierr)
      implicit none
      integer ierr
      real*4 dchn(143),dchnf(9,4),dchnb(12,5)
      real*4 ped,cst,ped1
      integer i,j,i1
      real*4 r,ranf,frac,flti
c      equivalence (w(1)  ,dchn(1)),
c     &            (w(144),dchnf(1,1)),
c     &            (w(180),dchnb(1,1))

      ierr=0
      cst=-2
      ped1=ped-660.
      i1=ped1/20.+1.e-6
      frac=(ped1-i1*20.)/20.
      i1=i1+1
      flti=dchn(i1)+(dchn(i1+1)-dchn(i1))*frac
      r=ranf()
      if(r.ge.flti) go to 190   !190 is forward
c     scatter bacward
      if(dchnb(1,1).gt.ped) go to 100
      do i=1,5
         if(dchnb(1,i).gt.ped) go to 40
      end do
      go to 100
 40   continue
      i1=i
      r=ranf()
C      frac=(ped-dchnb(1,i1))/(dchnb(1,i1+1)-dchnb(1,i1))
      frac=(ped-dchnb(1,i1-1))/(dchnb(1,i1)-dchnb(1,i1-1))
      do j=2,12
C         flti=dchnb(j,i1)+(dchnb(j,i1+1)-dchnb(j,i1))*frac
         flti=dchnb(j,i1-1)+(dchnb(j,i1)-dchnb(j,i1-1))*frac
         if(flti.ge.r) go to 80
      end do
      go to 100
 80   continue
      i1=j
      r=ranf()
      flti=i1-2 + r
      if(i1.le.3) then
         cst = 0.01*flti - 1.
      else if(i1.eq.4) then
         cst = 0.02*r - 0.98
      else if(i1.eq.5) then
         cst = 0.04*r - 0.96
      else if(i1.le.7) then
         cst = 0.06*flti - 1.16
      else if(i1.eq.8) then
         cst = 0.08*r - 0.8
      else if(i1.eq.9) then
         cst = 0.1*r - 0.72
      else if(i1.eq.10) then
         cst = 0.12*r - 0.62
      else if(i1.eq.11) then
         cst = 0.2*r - 0.5
      else
         cst = 0.3*(r-1.)
      end if
      return
 190  continue     !scattering forward
      if(dchnf(1,1).gt.ped) go to 100
      do i=1,4
         if(dchnf(1,i).gt.ped) go to 50
      end do
      go to 100
 50   continue
      i1=i
      r=ranf()
C      frac=(ped-dchnf(1,i1))/(dchnf(1,i1+1)-dchnf(1,i1))
      frac=(ped-dchnf(1,i1-1))/(dchnf(1,i1)-dchnf(1,i1-1))
      do j=2,9
C         flti=dchnf(j,i1)+(dchnf(j,i1+1)-dchnf(j,i1))*frac
         flti=dchnf(j,i1-1)+(dchnf(j,i1)-dchnf(j,i1-1))*frac
         if(flti.ge.r) go to 90
      end do
      go to 100
 90   continue
      i1=j
      r=ranf()
      flti=i1-2 + r
      if(i1.le.5) then
         cst = 1. - 0.025*flti
      else if(i1.eq.6) then
         cst = 0.85 + 0.05*r
      else if(i1.eq.7) then
         cst = 0.70 + 0.15*r
      else if(i1.eq.8) then
         cst = 0.50 + 0.2*r
      else
         cst = 0.5*r
      end if
      return
 100  ierr=1
      return
      end
