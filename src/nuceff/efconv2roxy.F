********************************************************************
*
*     -----------------------------------------------------
      real function efconv2roxy(rr)
*     -----------------------------------------------------
*     
*     ( purpose )
*       Convert r(target) into r(16O) by assuming
*       rho(r(target)) = rho(r(16O))
*
*     ( input ) 
*       rr       : radius(fm)
*       NUMATOM  : atomic number(from efpion.h)
*
*     ( output ) 
*       efconv2roxy : converted radius
*
*     (Creation Date and Author)
*       2007.11.05 ; G.Mitsuka
*       2009.02.13 ; Y.Hayato, G.Mitsuka (Fixed bug)
*
********************************************************************
      implicit none
#include "efpion.h"
      real rr,cc,alpha,ww
      real dfactoxy,afoxy,ccoxy ! parameter for 16O
      parameter (dfactoxy = 0.9985962)
      parameter (afoxy = 0.40961)
      parameter (ccoxy = 2.69)
      real buf1,buf2,buf3,buf4

C     get parameter from efpion.h
      cc    = CNN
      alpha = AF
      ww    = WPARM

      if (rr.ge.(max(ccoxy,cc)*2.5-0.01)) then
         efconv2roxy = 7.0 - 1.E-4
         return
      endif

C     Charge density for 16O is defined using fitted 2-parameter
C     fermi model, but for other target use 3-parameter Fermi model
      buf1=1.+exp((rr-cc)/alpha)
      buf2=1.+ww*(rr/cc)**2

      if (buf2.le.0) then
         efconv2roxy = 7.0 - 1.E-4
         return
      endif

      buf3=buf1*dfactoxy/buf2-1.

      if(buf3.gt.0.) then
         buf4=log(buf3)
         efconv2roxy=afoxy*buf4+ccoxy
      else
c         efconv2roxy=7.0-1.E-4
         efconv2roxy=1.e-6
      endif

C      write(*,*) "rr=",rr,"/efconv=",efconv2roxy
      if(efconv2roxy.lt.0.) efconv2roxy=1.e-6
      if(efconv2roxy.ge.7.) efconv2roxy=7.0-1.E-4

      return
      end
