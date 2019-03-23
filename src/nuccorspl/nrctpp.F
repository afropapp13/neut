      subroutine nrctpp(w,lk,ped,cst,ierr)
      implicit none
      integer lk,ierr
      real*4 w(lk,5)
      real*4 ped,cst
      integer i,j,i1
      real*4 r,ranf,frac,sum,flti

      ierr=0
      cst=-2

      do i=1,5
         if(w(1,i).gt.ped) go to 40
      end do
 100  ierr=1
      return
 40   continue
      sum=0
      i1=i-1
      r = ranf()
      frac=(ped - w(1,i1))/(w(1,i1+1) - w(1,i1))
      do j=2,lk
         sum = sum + w(j,i1)*(1.-frac) + w(j,i1+1)*frac
         if(r.lt.sum) go to 70
      end do
      go to 100
 70   continue
      r=ranf()
      i1=j+1
      if(lk.eq.11 .and. i1.le.3) then
         cst = 0.4 * r
         go to 150
      end if
      if(lk.eq.12) i1=i1-1
      flti = i1-2 + r
      if(i1.le.4) then
         cst = 0.2 * flti
      else if(i1.le.6) then
         cst = 0.3 + 0.1 * flti
      else if(i1.le.10) then
         cst = 0.6 + 0.04 * flti
      else
         cst = 0.78 + 0.02 * flti
      end if
      if(cst.gt.1.)cst=1.
 150  continue
      if(ranf().lt.0.5) cst = -cst
      return
      end
