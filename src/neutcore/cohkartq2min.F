******************************************
      real function cohkartq2min(s,Wmin)
*
*     (Purpose)
*       Return the kinematically allowed minimum Q2 value
*
*     (Input)
*       s     : Mn**2 + 2*Mn*Enu
*       Wmin  : (Mn+Mpi)**2
*         where Mn nuclus mass, Enu incident neutrino energy, and
*         Mpi pion mass.
*
*     (Output)
*       cohkartq2min:The kinematically allowed minimum Q2 value
*
*     (Creation Date and Author)
*       2007.01.08 ; G.Mitsuka -- first creation
*       2007.11.05 ; G.Mitsuka -- split this function to one program
*
      implicit none

#include "nework.h"
#include "necard.h"
      real s,Wmin,Mn,Mlep
      real buf1,buf2,buf3
      real cohkartlamda,Mpr,Mnt
      external cohkartlamda

      CALL MCMASSGV(2212,Mpr)
      CALL MCMASSGV(2112,Mnt)
      Mn = NUMBNDP*Mpr + NUMBNDN*Mnt

      CALL MCMASSGV(IPNE(3),Mlep)
      buf1 = (s-Mn**2)/2
C     In NC, Mlep = 0
      if(abs(modene).eq.36) then
         if(s.ge.Wmin) then
            cohkartq2min = 0.
         else
            cohkartq2min = -1.
         endif
C     In CC, lepton mass is considered
      else 
         if(cohkartlamda(1.,Mlep**2/s,Wmin/s).ge.0.) then
            buf2 = 1.-sqrt(cohkartlamda(1.,Mlep**2/s,Wmin/s))
            buf3 = Wmin+Mlep**2-Mn**2/s*(Wmin-Mlep**2)
            cohkartq2min = buf1*buf2-(1./2.)*buf3
         else
            cohkartq2min = -1.
         endif
      endif
      
      return
      end
