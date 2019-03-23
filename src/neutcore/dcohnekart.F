************************************************************************
*     ----------------------------
      Real Function dcohnekart(E,Q2,nu,t,x,y,z)
*     ----------------------------
*
*     (Purpose)
*       Return differential corss-section depending on E,Q2,nu,t
*
*     (Input)
*       E     : Incident Neutrino Energy ( GeV )
*       Q2    : Momentum Transfer from leptonic sector (GeV^2)
*       nu    : Energy Transfer (GeV)
*       t     : Momentum Transfer in hadronic system (GeV^2)
*
*     (Output)
*       dcohnekart : differential corss-section(dsigma/(dQ2*dnu*dt))
*     
*     (Creation Date and Author)
*       2007.01.08 ; G.Mitsuka
*       2007.11.05 ; G.Mitsuka -- function and subroutine are replaced 
*                                 to each fortran file
*
*     (Comment)
*       Kinematics of this routine is based on
*       Physical Revew D 74, 054007(2006) A.Kartavtsev et al.
*
************************************************************************
      implicit none

#include "nework.h"
#include "necard.h"
#include "pdf804/w50516.inc"

      real E,Q2,nu,t,Q2max,Mn,xi,tmax
      real Gf,Vud,fpi,L00,Lll,Ll0,Mpi
      real Mlep,q32,Pi,Mpr,Mnt
      real buf1,buf2,buf3,x,y,z,BRNTC2
      real*8 fndsigmapidt
      external fndsigmapidt
      
      PARAMETER (     xi  = 3.            ) ! Not used
      PARAMETER (     Pi  = 3.14159265    )
      PARAMETER (     Gf  = 1.16637E-5    )
      PARAMETER (     fpi = 0.092         )
      PARAMETER (     Vud = 0.9740        )
      parameter (     Q2max = 0.2         )
      parameter (     tmax = 0.05         )
      PARAMETER (   BRNTC2 = 1.E-24       )

C     Initialize
      CALL MCMASSGV(IPNE(5),Mpi)
      CALL MCMASSGV(2212,Mpr)
      CALL MCMASSGV(2112,Mnt)
      Mn = NUMBNDP*Mpr + NUMBNDN*Mnt
      dcohnekart = 0.

      buf2 = fpi**2/Q2      
      q32 = Q2+nu**2

C     In NC, lepton(neutrino) mass is negligible
      if(abs(modene).eq.36) then
         buf1 = Gf**2*nu/(4*(2*Pi)**2*E**2)

         L00 = 2*(Q2*(2*E-nu))**2/(Q2*q32)-2*Q2
         buf3 = L00
C     In CC, lepton(e,mu,tau) mass correction should be applied
      else
         CALL MCMASSGV(IPNE(3),Mlep)
         buf1 = Gf**2*Vud**2*nu/(2*(2*Pi)**2*E**2)

C        density matrix elements
         L00 = 2*(Q2*(2*E-nu)-nu*Mlep**2)**2/(Q2*q32)-2*(Q2+Mlep**2)
         Lll = 2*Mlep**2*(Mlep**2/Q2+1.)
         Ll0 = 2*Mlep**2*(Q2*(2*E-nu)-nu*Mlep**2)/(Q2*sqrt(q32))
         if(L00.lt.0..or.Lll.lt.0..or.Ll0.lt.0.) goto 100
         buf3 = L00 + Lll*(Mpi**2/(Q2+Mpi**2))**2
     &        + 2*Ll0*Mpi**2/(Q2+Mpi**2)
      endif

C     dcohnekart is the differential cross section off Carbon nucleus,
C     where fndsigmapidt is dsigma(pi+12C->pi+12C)/dt
      dcohnekart = buf1*buf2*buf3*fndsigmapidt(dble(nu),dble(t))
     &     *BRNTC2*1.D+40
 100  return
      end
