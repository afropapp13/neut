*     -----------------------------------------------------
      subroutine efdsg(ireaction,tpi,thet,dsg,ierr)
*     -----------------------------------------------------
c
*     ( purpose )
c     Routine that calculates the differential cross section for
c     all possible pion-nucleon reactions, at a given incident pion
c     lab energy (tpi) and cm scattering angle (thet).
c     The possible reactions are:
c       1) pi+ p -> pi+ p , pi- n -> pi- n
c       2) pi+ n -> pi+ n , pi- p -> pi- p
c       3) pi+ n -> pi0 p , pi- p -> pi0 n
c       4) pi0 p -> pi0 p , pi0 n -> pi0 n
c
*     ( input ) 
*       ireaction : see above  
*       tpi       : pion energy in lab frame (MeV)
*       thet      : pion/nucleon CM angle (degrees)
*
*     ( output ) 
*       dsg       : differential cross section
*       ierr      : error flag
*
*     ( history )
*       Creation:     R. Tacik
*       2010-06: P.de Perio - Modify to use TLAB up to 1995 MeV
*                             in steps of 15 MeV
*                           - Return dsg=0 upon failure instead of 'stop'
*                             by adding IERR flag

      implicit none
      
      integer ireaction, j, ierr
      
      real mp/938.28/,mpi/139.57/
      real hbarc/197.33/
      real tpi,thet,dsg
      
      real epi,ppi,bet,gam,kpi
      real p1,p2,p11,p21
 
      real f(2),g(2),f2,g2
 
      real sind,cosd
      external sind, cosd

C      real as(2),app(2),apm(2),adp(2),adm(2)      
C      common /amps/   as,app,apm,adp,adm

#include "efdsg.h"

C  Initialize
      dsg = 0.
      ierr = 0
      
C  DO SOME KINEMATICS

      epi = mpi + tpi
      ppi = sqrt(epi**2 - mpi**2)
      bet = ppi/(epi+mp)
      gam = sqrt(1./(1.-bet**2))
      kpi = -gam*bet*epi + gam*ppi
      kpi = kpi/hbarc
      
C  GET THE NUCLEAR AMPLITUDES
     
      call efdsgamp(tpi,ireaction,ierr)

      if (ierr.ne.0) then
         return
      end if

C  EVALUATE POLYNOMIALS
      
      P1     = cosd(THET)
      P2     = 0.5*( 3.*cosd(THET)**2 - 1. )
      P11    = sind(THET)
      P21    = 3.*sind(THET)*cosd(thet)

c      write (6,*) ' p ',p1,p2,p11,p21
      
C  FIRST DO L=2 PART

      do J = 1, 2
         F(J) = (3.*ADP(J)+2.*ADM(J))*P2
         G(J) = (ADP(J)-ADM(J))*P21
      end do
         
c      write (6,*) f(1),f(2),g(1),g(2)
      
C  ADD ON THE DO L=1 PART

      do J = 1, 2
         F(J) = F(J) + (2.*APP(J)+APM(J))*P1
         G(J) = G(J) + (APP(J)-APM(J))*P11
      end do

c      write (6,*) f(1),f(2),g(1),g(2)      
      
C  AND FINALLY THE L=0 PART

      DO  J = 1, 2
         F(J) = F(J) + AS(J)
      end do

c      write (6,*) f(1),f(2),g(1),g(2)      
      
C  NORMALIZE 

      DO J = 1, 2
         F(J) = F(J) / kpi
         G(J) = G(J) / kpi
      end do
      
c      write (6,*) f(1),f(2),g(1),g(2)      
      
C  FINAL RESULT

      F2 = F(1)*F(1) + F(2)*F(2)
      G2 = G(1)*G(1) + G(2)*G(2)

      DSG = 10.*( F2 + G2 )
     
      RETURN
      END

