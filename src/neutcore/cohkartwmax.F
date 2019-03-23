******************************************
      real function cohkartwmax(Q2,s)
*
*     (Purpose)
*       Return the kinematically allowed maximum W value
*
*     (Input)
*       Q2    : Momentum transfer from leptonic sector
*       s     : Mn**2 + 2*Mn*Enu
*         where Mn nuclus mass, Enu incident neutrino energy
*
*     (Output)
*       cohkartwmax:The kinematically allowed maximum W value
*
*     (Creation Date and Author)
*       2007.01.08 ; G.Mitsuka -- first creation
*       2007.11.05 ; G.Mitsuka -- split this function to one program
*
      implicit none

#include "nework.h"
#include "necard.h"
      real buf1,buf2,buf3
      real s,Mn,Mlep,Q2,Mpr,Mnt
      
      CALL MCMASSGV(2212,Mpr)
      CALL MCMASSGV(2112,Mnt)
      Mn = NUMBNDP*Mpr + NUMBNDN*Mnt

      CALL MCMASSGV(IPNE(3),Mlep)

C     In NC, Mlep = 0
      if(abs(modene).eq.36) then
         cohkartwmax = s*(1.-Q2/(s-Mn**2))
C     In CC, lepton mass is considered
      else
         buf1 = (1./4.)*s**2*(1.-Mn**2/s)**2*(1.-Mlep**2/s)
         buf2 = (Q2-(s/2.)*(1.-Mn**2/s)+(Mlep**2/2.)*(1.+Mn**2/s))**2
         buf3 = (1.-Mn**2/s)*(Q2+Mlep**2)      
         cohkartwmax = (buf1-buf2)/buf3
      endif
      
      return
      end
