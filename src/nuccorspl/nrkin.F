          subroutine nrkin(vin,pin,am,pm,pout,pn,eout)
          implicit none

          real*4 pin,am,pm,pn,eout

          real*4 vin(3),pout(3)

#include "nrcms.h"
C         common/nuccms/rhon,pfermi,unucl,ecms2,up4cms(3),beta(3),prot
C         logical prot


          real*4 p1in(3),p2in(3),pdum(3)

          real*4 bm2
          parameter (bm2=939.**2)

          integer*4 i
          real*4    e1,e2,edum

          do i=1,3
                 p1in(i)=pin*vin(i)
                 p2in(i)=-p1in(i)
          end do
          e1=sqrt(pin**2+am**2)
          e2=sqrt(pin**2+bm2)
          call nrllor(pout,eout,p1in,e1,beta)
          call nrllor(pdum,edum,p2in,e2,beta)
          pn=sqrt(eout**2-am**2)
          pm=sqrt(edum**2-bm2)
          return
          end
