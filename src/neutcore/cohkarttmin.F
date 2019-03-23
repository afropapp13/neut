******************************************
      real function cohkarttmin(Q2,nu)
*
*     (Purpose)
*       Return the kinematically allowd minimum t value
*
*     (Input)
*       Q2    : Momentum transfer from leptonic sector
*       nu    : Energy transfer
*
*     (Output)
*       The kinematically allowd minimum t value
*
*     (Creation Date and Author)
*       2007.01.08 ; G.Mitsuka -- first creation
*       2007.11.05 ; G.Mitsuka -- split this function to one program
*
      implicit none

#include "nework.h"
      real Q2,nu,Mn,Mpi,W2,cohkartlamda
      external cohkartlamda
      real buf1,buf2,Mpr,Mnt

      CALL MCMASSGV(IPNE(5),Mpi)
      CALL MCMASSGV(2212,Mpr)
      CALL MCMASSGV(2112,Mnt)
      Mn = 8*(Mpr+Mnt)
      W2 = Mn**2-Q2+2*Mn*nu
      if(cohkartlamda(W2,-Q2,Mn**2).ge.0.
     &     .and.cohkartlamda(W2,Mpi**2,Mn**2).ge.0.) then
         buf1 = (Q2+Mpi**2)**2
C         buf2 = sqrt(cohkartlamda(W2,-Q2,Mn**2))
C     &        -cohkartlamda(cohkartlamda(W2,Mpi**2,Mn**2))
         buf2 = sqrt(cohkartlamda(W2,-Q2,Mn**2))
     &        -sqrt(cohkartlamda(W2,Mpi**2,Mn**2))
         cohkarttmin = (buf1-buf2**2)/(4*W2)
      else
         cohkarttmin = 1000000.
      endif

      return
      end
