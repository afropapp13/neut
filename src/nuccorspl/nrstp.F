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

          if(icont.eq.1)then
                 do i=1,3
                        start(i)=0
                 end do
          else if(icont.eq.0)then
                 xx=ranf()
                 if(mode.eq.1)then              !surface mode
                        r=rmax*sqrt(xx)
                        phi=6.283185*ranf()
                        start(1)=r*cos(phi)
                        start(2)=r*sin(phi)
                        start(3)=-rmax
                 else
                        if(mode.eq.2)then       !volume mode
                           r=nrprbin(probmi,nbin,rmin,rbin,xx)
                        else                            !Dover mode
                           r=nrprbin(probdi,nbin,rmin,rbin,xx)
                        end if
                        call nrranve(r,start,2.,-1.)
                 end if
          end if
          return
          end
