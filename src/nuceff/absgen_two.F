c ---------------------------------------------------------------
   
      subroutine absgen_two(tpir,ct)

      implicit none

      real*4 tpir,f,ct
C      real*8 tpid,x,pidpp
      real*8 x
      real*8 epsin/1.0d-7/
      real*8 xmin,xmax,dx,tot,r,v

      real*8   abspidpp, dgauss
      external abspidpp, dgauss

#include "abspidpp.h"
C      common /pi_energy/ tpid

      integer*4 i

      tpid = tpir

      xmin = 0.d0
      xmax = 1.d0

      tot = dgauss(abspidpp,xmin,xmax,epsin)

      call ranlux(f,1)

      dx   = 1./2.
      xmin = 0.d0
      xmax = xmin + dx
      v    = dgauss(abspidpp,xmin,xmax,epsin)
      r = v/tot

      do i = 1, 15
         r = v/tot
         dx = dx/2.d0
         if (f.lt.r) then
            xmax = xmax - dx
            v    = dgauss(abspidpp,xmin,xmax,epsin)
         else if (f.gt.r) then
            xmax = xmax + dx
            v    = dgauss(abspidpp,xmin,xmax,epsin)
         else if (f.eq.r) then
            go to 200
         end if
      end do

 200  ct = xmax

      end
