          subroutine nrstp(start,icont,mode)

          implicit none
          real*4 start(3)
          integer*4 icont,mode

#include "nrdis.h"
C         common /nucdis/ probmi(31),probdi(31),rhotab(30),pnorm

          real*4    rbin,rmin,rmax
          integer*4 nbin
          parameter ( rbin=.2, rmin=0, rmax=5.99, nbin=30)

          integer*4 i
          real*4    xx,r,phi

          real*4    ranf,nrprbin
          external  ranf,nrprbin

c     icont 0 = generate interaction point, icont 1 = start from 0,0,0
          if(icont.eq.1)then
                 do i=1,3
                        start(i)=0
                 end do
c     xx set to random value, then r is root(xx) fraction of rmax
          else if(icont.eq.0)then
                 xx=ranf()
c     mode = 1 is surface
c     NB. looks like everything runs in surface mode
c     probmi/probdi don't seem to be set anywhere in NEUT
                 if(mode.eq.1)then              !surface mode
                        r=rmax*sqrt(xx)
c     phi is a random number times 2pi
                        phi=6.283185*ranf()
c     set start from phi and r
                        start(1)=r*cos(phi)
                        start(2)=r*sin(phi)
                        start(3)=-rmax
                 else
c     mode = 2 is volume
                        if(mode.eq.2)then       !volume mode
                           r=nrprbin(probmi,nbin,rmin,rbin,xx)
c     mode = other is Dover mode 
c     not sure what probmi/probdi are
                        else                            !Dover mode
                           r=nrprbin(probdi,nbin,rmin,rbin,xx)
                        end if
c     nrranve provides a random vector of length r
                        call nrranve(r,start,2.,-1.)
                 end if
          end if
          return
          end
