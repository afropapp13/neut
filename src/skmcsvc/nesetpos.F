      subroutine nesetpos(rad,pos3,pos)

      implicit none
      real*4 pos(3), rad, pos3

      real*4   rlu
      external rlu

      integer*4 idum

      pos(3)=pos3

 10   continue
      
         pos(1)=(1.-rlu(idum)*2)*rad
         pos(2)=(1.-rlu(idum)*2)*rad

         if (sqrt(pos(1)**2+pos(2)**2).gt.rad) then
            goto 10
         endif

 100  continue

      end
      
   
