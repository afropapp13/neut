C********************************************************************

      subroutine efdsgamp(tpi,ireaction,ierr)
*     -----------------------------------------------------
c
*     ( purpose )
C     Called by efdsg.F
C
c     The possible reactions are:
c       1) pi+ p -> pi+ p , pi- n -> pi- n
c       2) pi+ n -> pi+ n , pi- p -> pi- p
c       3) pi+ n -> pi0 p , pi- p -> pi0 n
c       4) pi0 p -> pi0 p , pi0 n -> pi0 n
c
*     ( input ) 
*       ireaction : see above  
*       tpi       : pion energy in lab frame (MeV)
*
*     ( output ) 
*       dsg       : differential cross section
*       ierr      : error flag
*
*     ( history )
*       Creation:     R. Tacik
*       2010-06: P.de Perio - Add IERR flag

      implicit none
      
      integer ireaction, j, ierr

      real sind, cosd
      external sind, cosd
      
      real tpi
      real as11(2),as31(2)
      real ap11(2),ap13(2),ap31(2),ap33(2)
      real ad13(2),ad33(2),ad15(2),ad35(2)

C      real xs11,xs31,xp11,xp31,xp13,xp33,xd13,xd33,xd15,xd35   
C      common /phases/ xs11,xs31,xp11,xp31,xp13,xp33,xd13,xd33,xd15,xd35
C      real as(2),app(2),apm(2),adp(2),adm(2)
C      common /amps/   as,app,apm,adp,adm

#include "efdsg.h"

      ierr = 0

      call efdsgps(tpi,ierr)
      
      if(ierr.ne.0) then
         return
      end if

      AS11(1) = sind(xs11)*cosd(xs11)
      AS11(2) = sind(xs11)**2
      
      AS31(1) = sind(xs31)*cosd(xs31)
      AS31(2) = sind(xs31)**2
      
      AP11(1) = sind(xp11)*cosd(xp11)
      AP11(2) = sind(xp11)**2

      AP13(1) = sind(xp13)*cosd(xp13)
      AP13(2) = sind(xp13)**2

      AP31(1) = sind(xp31)*cosd(xp31)
      AP31(2) = sind(xp31)**2

      AP33(1) = sind(xp33)*cosd(xp33)
      AP33(2) = sind(xp33)**2
      
      AD13(1) = sind(xd13)*cosd(xd13)
      AD13(2) = sind(xd13)**2
      
      AD33(1) = sind(xd33)*cosd(xd33)
      AD33(2) = sind(xd33)**2
      
      AD15(1) = sind(xd15)*cosd(xd15)
      AD15(2) = sind(xd15)**2
      
      AD35(1) = sind(xd35)*cosd(xd35)
      AD35(2) = sind(xd35)**2
      
      if (ireaction.eq.1) then
         do j = 1, 2
            AS(J)  = AS31(J)
            APP(J) = AP33(J)
            APM(J) = AP31(J)
            ADP(J) = AD35(J)
            ADM(J) = AD33(J)
         end do
      else if (ireaction.eq.2) then
         do j = 1, 2
            AS(J)  = ( AS31(J) + 2.*AS11(J) ) / 3.
            APP(J) = ( AP33(J) + 2.*AP13(J) ) / 3.
            APM(J) = ( AP31(J) + 2.*AP11(J) ) / 3.
            ADP(J) = ( AD35(J) + 2.*AD15(J) ) / 3.
            ADM(J) = ( AD33(J) + 2.*AD13(J) ) / 3.
         end do
      else if (ireaction.eq.3) then
         do j = 1, 2
            AS(J)  = sqrt(2.) * ( AS31(J) - AS11(J) ) / 3.
            APP(J) = sqrt(2.) * ( AP33(J) - AP13(J) ) / 3.
            APM(J) = sqrt(2.) * ( AP31(J) - AP11(J) ) / 3.
            ADP(J) = sqrt(2.) * ( AD35(J) - AD15(J) ) / 3.
            ADM(J) = sqrt(2.) * ( AD33(J) - AD13(J) ) / 3.
         end do
      else if (ireaction.eq.4) then
         do j = 1, 2
            AS(J)  = ( 2.*AS31(J) + AS11(J) ) / 3.
            APP(J) = ( 2.*AP33(J) + AP13(J) ) / 3.
            APM(J) = ( 2.*AP31(J) + AP11(J) ) / 3.
            ADP(J) = ( 2.*AD35(J) + AD15(J) ) / 3.
            ADM(J) = ( 2.*AD33(J) + AD13(J) ) / 3.
         end do
      end if
     
      return
      end
