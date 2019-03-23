c-------------------------------------------------------------

      real*8 function abspidpp(x)

      implicit none

      real*8 a,b,f
      real*8 c1,g1,e1
      real*8 c2,g2,e2
      real*8 mpi,md

      real*8 aa2,aa4
      real*8 bb2,bb4
      real*8 cc2,cc4,cc6
      real*8 ga2,ga4,ga6
      real*8 ep2,ep4,ep6
      real*8 etat,eta0, eta
      real*8 tau2,tau4

      real*8 x,x2,x4
      real*8 alf0,alf2,alf4
      real*8 baseline,e
      real*8 g12,g22,ga22,ga42,ga62
      real*8 p2,p4
      real*8 ppi,s,ss,wpi
      real*8 sig
      real*8 w1,w2,w4
      

#include "abspidpp.h"
C      common /pi_energy/ tpid

      a   = -0.577d0
      b   = 3.09d0
      c1  = 12.50d0
      g1  = 69.0d0
      e1  = 2133.9d0
      c2  = 0.075d0
      g2  = 54.d0
      e2  = 2644.d0
      f   = 0.38d-3
      mpi = 139.57d0
      md  = 1875.63d0
      
      g12 = g1*g1
      g22 = g2*g2

      aa2 = -2.d-3
      aa4 = -40.d-3
      bb2 = 0.6d-3
      bb4 = 2.1d-3
      cc2 = 2.8d0
      cc4 = -0.255d0
      cc6 = -0.06d0
      ga4 = 44.d0
      ga6 = 245.d0
      ep2 = 2126.d0
      ep4 = 2184.d0
      ep6 = 2100.d0
      etat = 0.5d0
      eta0 = 1.398d0

      ga42 = ga4*ga4
      ga62 = ga6*ga6

C      write(*,*)'mpi=',mpi,'md=',md,'tpid=',tpid

      s   = mpi**2 + md**2 + 2.*(tpid+mpi)*md
     
C      write(*,*)'s=',s

      ss  = sqrt(s)

C      write(*,*)'ss=',ss
      
      wpi = (s+mpi**2-md**2)/(2.*ss)

C      write(*,*)'wpi=',wpi

      ppi = sqrt(wpi*wpi-mpi*mpi)

C      write(*,*)'ppi=',ppi

      eta = ppi/mpi
      e   = sqrt((mpi+md)**2+2.*md*tpid)
      w1  = g12/((e-e1)**2+g12)
      w2  = g22/((e-e2)**2+g22)
      sig = a + b/sqrt(tpid) + c1*w1 + c2*w2 + f*tpid

      alf0 = sig/(2.*acos(-1.))

      ga2  = 67.d0 - 14.d0*(eta-eta0)
      ga22 = ga2*ga2
      w2   = ga22/((e-ep2)**2+ga22)
      tau2 = eta/(eta+etat)
      alf2 = tau2*(aa2 + bb2*eta*eta + cc2*w2)

      w4   = ga42/((e-ep4)**2+ga42)
      tau4 = eta**3/(eta**3+etat**3)
      if (tpid.le.750.) then
         alf4 = tau4*(aa4+bb4*eta*eta+cc4*w4)
      else if (tpid.gt.750.) then
         alf4 = 0.0
      end if

      baseline = -alf2/2.d0 + 3.d0*alf4/8.d0
      x2 = x*x
      x4 = x2*x2
      p2 = (3.d0*x2 - 1.)/2.d0
      p4 = (35.d0*x4-30.*x2+3.d0)/8.d0
      abspidpp = alf2*p2 + alf4*p4 - baseline
      abspidpp = abspidpp*sqrt(1.d0-x2)

      end
