      subroutine nrctpnl(wl,wh,ped,cst,ierr)
      implicit none
      real*4 wl(5,16),wh(5,23),ped,cst
      integer ierr,i1,i
      real*4 ws(5),r,ranf,val,fact

      val=ped
      if(ped.gt.300) then
         val=ped-300.
         i1=val/20.
         if (i1.ge.22) then
            i1 = 21
         endif
      else
         i1=val/20.
         if (i1.ge.15) then
            i1 = 14
         endif
      endif
      fact=(val - i1*20.)/20.
      i1=i1+1
      do i=1,5
         if(ped.le.300.) then
            ws(i)=wl(i,i1)+(wl(i,i1+1)-wl(i,i1))*fact
         else
            ws(i)=wh(i,i1)+(wh(i,i1+1)-wh(i,i1))*fact
         end if
      end do
      r=ranf()
      if(ranf().ge.ws(1)) then
         if(ranf().lt.ws(4)) then
            cst = r
         else
            do i=1,3
               fact=ranf()
               if(fact.gt.r) r=fact
            end do
            cst = r
         end if
      else
         if(ranf().lt.ws(2)) then
            cst = -r
         else
            do i=1,6
               fact=ranf()
               if(fact.gt.r) r=fact
            end do
            cst = -r
         end if
      end if
      return
      end
