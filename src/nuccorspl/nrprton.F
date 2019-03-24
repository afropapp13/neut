      subroutine nrprton ( anuc, !A of the nucleus =16 for Oxygen
     $     stpi,                !current point inside nucleus
     $     pi,idmc,             !4momentum in GeV,id (should be p or n)
     $     po,ido,stpo,no,      !the same for outgoing, creat point,#of track
     $     imode,icont)

c
c     ( Creation Date and Author )
c       2007.11.05 ; G.Mitsuka - add support other target than Oxygen
c                                other target is treated by scaling to Oxygen
c
c       2014.10.22 ; M.Hierholzer - add nucleon FSI history block to calculate systematics;
c                                   also add some documentation of the variables used in this file
c
c
c     icont =  0        generate interaction point
c              1        start from stpt  (set start(i)=0)
c     imode =  1        start on the surface
c              2        start in the volume  (probmi)
c              others   Dover mode       (probdi)
c
c     When NEUT is run by neutgeom/genev, we always have icont=1 and imode=2.
c

      implicit none
      real*4 stpi(3),pi(4),anuc
      real*4 po(4,20),stpo(3,20)
      integer ido(20)
      integer imode,icont,idmc,no

#include "vcwork.h"
#include "nework.h"
#include "neutparams.h"
#include "nrnuclparam.h"
#include "nrint.h"
#include "nrcard.h"
C     common /nucint/uin,chint,ein,pin(3),upin(3)
C     real*4 uin,chint,ein,pin,upin
#include "nrcms.h"
C     common /nuccms/rhon,pfermi,unucl,ecms2,up4cms(3),beta(3),prot
C     real*4 rhon,pfermi,unucl,ecms2,up4cms,beta
C     logical prot
#include "nrnewprd.h"
C     common /nucnewprd/piso(4,4),chiso(4)
C     real*4 piso,chiso
#include "nrstack.h"
C     common /nucstack/ntrack,
C     1             pptrack(4,20),chtrack(20),amass(20),stpt(3,20)
C     integer ntrack
C     real*4  pptrack,chtrack,amass,stpt

#include "nrhad.h"
C     common /nuchad/ihap,nel(3)
C     integer ihap,nel
#include "nrhadr.h"
C     common /hadron/inot,itrans,idec,ine1,ine2,iel,iabs,icx
C     integer inot,itrans,idec,ine1,ine2,iel,iabs,icx

#include "nucleonfsihist.h"

!
!     Used external inputs from common blocks:
!       nrrmsrad from nrnuclparam.h - unknown meaning, depends on target, set in neutcore/nrsettarg.F
!

!
!     Output variables to common blocks:
!       chtrack(*) to nrstack.h - charge of nucleon in track
!       amass(*) to nrstack.h - mass of nucleon in track
!       stpt(3,*) to nrstack.h - starting point of track inside nucleus
!       pptrack(4,*) to nrstack.h - initial 4-momentum of nucleon in track
!       ntrack to nrstack.h - number of tracks
!       ihap to nrhad.h - ??
!       nel(3) to nrhad.h - ?? 
!       itrans to nrhadr.h - ??
!

!
!     Information passed on to (and back from) subroutines called in this routine via common blocks:
!       uin to nrint.h - mass of current track
!       chint to nrint.h - charge of current track
!       ein to nrint.h - energy of current track
!       pin(3) to nrint.h - momentum of current track
!       upin(3) to nrint.h - direction of current track (never actually used outside this function)
!       rhon to nrcms.h, set in nrfermi - nuclear density at current r
!       pfermi to nrcms.h, set in nrfermi - fermi momentum
!       unucl to nrcms.h , set in nrfermi- nucleon mass minus binding energy and fermi energy (?)
!       ecms2 to nrcms.h, set in nrfermi - center-of-mass energy squared (with fermi motion smearing?)
!       up4cms(3) to nrcms.h, set in nrfermi - direction of track (with fermi motion smearing?)
!       beta(3) to nrcms.h, set in nrfermi - unused?
!       prot to nrcms.h, set in nrfermi - decide if interaction is on proton (.true.) or neutron (.false.)
!

!
!     Purely internal variables:
!       rmax - nuclus radius, calculated from nrrmsrad in dependence of target
!       pptot - absolute initial nucleon momentum
!       start(3) - offset (?) of starting point for tracking, w.r.t. stpt()
!       itdone - track currently being processed
!       xtr(3) - current position (of current step) inside nucleus
!       r - current position's radius, calculated from xtr()
!       loop - some counter used in surface mode
!       iin - flag if inside nucleus, used in surface mode
!


      real*4 xtr(3),start(3)

      real step,rmax,amn,pptot
      integer j,itdone
      real*4 pinss,pinm,r,az
      integer loop,iin,ides,ia,l1,l

      parameter (step=.2, amn=939. )

C      integer*4 nclux,ncint,nck1,nck2

      REAL*4 RMSRADOXY
      parameter (RMSRADOXY = 2.730)
      REAL effrmgas


      real*4 ptot, A,B,C
      common /martindebug/  ptot,A,B,C



C     rmax is set by scaling to 16O (rmax=6. for Oxygen)
      rmax = 6.*(NRRMSRAD/RMSRADOXY)




      if(idmc.eq.2212) then     !put the input track on stack
         chtrack(1)=1
      else if(idmc.eq.2112) then
         chtrack(1)=0
      else
         print *,'*** NUCPRTON ERROR call for wrong track id=',idmc
         return
      end if

      amass(1) = amn
      pptot=0
      r=0
      do j=1,3
         stpt(j,1)=stpi(j)
         r=r+stpi(j)**2
      end do
      r=sqrt(r)
      do j=1,4
         pptrack(j,1)=pi(j)
         if(j.lt.4) pptot = pptot + pptrack(j,1)**2
      end do
      pptot = sqrt(pptot)

C     Asmita R.:For MEC events fermi surface momentum is function of radial distance r (22-07-2013)
C     while for other modes its read from card 
C Change to local Fermi Gas ( 2017/04/03 )
      if (nucresflg.eq.1) then
         if (abs(MODENE).EQ.2) then 
            if (r.ge.rmax) then
               pfermi=EFFRMGAS(0,0,rmax)
            else
               pfermi=EFFRMGAS(0,0,r)
            endif
         else
            pfermi=NRPFSURF*1000. ! NRPFSURF is GeV units
         endif
      else if (nucresflg.eq.2) then
         if (r.ge.rmax) then
            pfermi=EFFRMGAS(0,0,rmax)
         else
            pfermi=EFFRMGAS(0,0,r)
         endif
      else
         write(*,*) 'Unknown NUCRESFLG ( expected 1 or 2 )',NUCRESFLG
         STOP
      endif

      if(pptot.lt.pfermi) then  !if it is below fermi sphere get out
         do j=1,4
            po(j,1)=pi(j)
            if(j.lt.4) stpo(j,1)=stpi(j)
         end do
         ido(1)=idmc
         no=1
         return
      end if


      call nrstp(start,icont,imode) !fix the starting point

      itdone = 0                !current track
      ntrack = 1                !number of tracks to go
      no = 0                    !number of output tracks
      itrans = 0
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 300  continue                  ! loop thru the input tracks
      itdone=itdone+1
      if(itdone.gt.ntrack) goto 100 !just return
      if(itdone.le.1)then       ! reset interaction markers
         ihap=-1
         nel(1)=0
         nel(2)=0
         nel(3)=0
      end if

      uin=amass(itdone)         !put a new track to follow
      chint=chtrack(itdone)
      ein = pptrack(4,itdone)
      pinss=0
      do j=1,3
         if(itdone.gt.1)start(j)=0.
         xtr(j)=stpt(j,itdone)+start(j)
         pin(j)=pptrack(j,itdone)
         pinss=pinss+pin(j)**2
      end do
      pinm=sqrt(pinss)
      do j=1,3
         upin(j)=pin(j)/pinm
      end do
      if(abs(uin-amn).gt.5.) goto 99 !not a nucleon

      loop=0
      iin=0

!
!     fill nucleon fsi history stack
!
      if(NFnvert.lt.maxnucleonvert.and.NFnstep.lt.maxnucleonstep) then
        NFnvert = NFnvert+1                                   ! increase number of vertices
        NFiflag(NFnvert) = 0 + chint*100                      ! interaction type and charge of propagated nucleon
        NFx(NFnvert) = xtr(1)                                 ! vertex position inside nucleus
        NFy(NFnvert) = xtr(2)
        NFz(NFnvert) = xtr(3)
        NFpx(NFnvert) = pin(1)                                ! momentum
        NFpy(NFnvert) = pin(2)
        NFpz(NFnvert) = pin(3)
        NFe(NFnvert) = ein                                    ! energy
        NFfirststep(NFnvert) = NFnstep+1                      ! index of first step of track segment

        NFnstep = NFnstep+1                                   ! increase number of steps
        NFecms2(NFnstep) = 0.                                 ! CMS energy squared (not a collision, so it's zero)
        NFptot(NFnstep) = 0.                                  ! total probability (to be removed)
      else
        print *,'nrprton.F: too many nucleon FSI vertices (',NFnvert,') or steps (',NFnstep,'), not saving all of them!'
      end if

!
!     track nuclon through nucleus
!
 200  continue

      r=0
      do j=1,3
         xtr(j)=xtr(j)+upin(j)*step
         r=r+xtr(j)**2
      end do
      r=sqrt(r)

      if(imode.eq.1)then
         if(iin.eq.0.and.r.gt.rmax)then
            loop=loop+1
            if(loop.gt.3) goto 99
            az=rmax**2-xtr(1)**2-xtr(2)**2
            if(az.lt.0.) goto 99
            xtr(3)=-sqrt(az)-0.18+0.05*loop
            goto 200
         end if
         if(r.le.rmax)iin=1
      end if

      if(r.gt.rmax) then        !particle got out of the nucleus
         if(itdone.le.ntrack) then !count #of interactions
            itrans=itrans+1
            if(ihap.lt.0)inot=inot+1
            if(nel(1).gt.0)iel=iel+nel(1)
            if(nel(2).gt.0)ine1=ine1+nel(2)
            if(nel(3).gt.0)ine2=ine2+nel(3)
         end if
         goto 99               ! ------>  put on output stack
      end if

      call nrfermi(r)           !set target nucleon
C      if(ecms2.lt.(2*amn)**2) goto 200  ! kinematically impossible / test is done inside nrnuc!

      call nrnuc(ides)          ! make interaction
      
      if(ides.ge.10) then       ! Paulu blocked interaction: save to fsi history stack
        if(NFnvert.lt.maxnucleonvert.and.NFnstep.lt.maxnucleonstep) then
          NFnvert = NFnvert+1                                   ! increase number of vertices
          NFiflag(NFnvert) = ides-10 + chint*100 + 1000         ! Pauli-blocked interaction type and charge of propagated nucleon
          if(prot) NFiflag(NFnvert) = NFiflag(NFnvert) + 10     ! charge of target nucleon
          NFx(NFnvert) = xtr(1)                                 ! vertex position inside nucleus
          NFy(NFnvert) = xtr(2)
          NFz(NFnvert) = xtr(3)
          NFpx(NFnvert) = pin(1)                                ! momentum
          NFpy(NFnvert) = pin(2)
          NFpz(NFnvert) = pin(3)
          NFe(NFnvert) = ein                                    ! energy
          NFfirststep(NFnvert) = NFnstep+1                      ! index of first step of track segment

          NFnstep = NFnstep+1                                   ! increase number of steps
          if(prot) then
            NFecms2(NFnstep) = ecms2                            ! CMS energy squared (positive: proton)
          else
            NFecms2(NFnstep) = -ecms2                           ! CMS energy squared (negative: neutron)
          end if
          NFptot(NFnstep) = ptot                                ! total probability (to be removed)
        else
          print *,'nrprton.F: too many nucleon FSI vertices (',NFnvert,') or steps (',NFnstep,'), not saving all of them!'
        end if
        goto 200                ! do the next step
      end if

      if(ides.gt.0) then        ! no interaction: do next step
        if(NFnstep.lt.maxnucleonstep) then       ! save per-step information
          NFnstep = NFnstep+1                                   ! increase number of steps
          if(prot) then
            NFecms2(NFnstep) = ecms2                            ! CMS energy squared (positive: proton)
          else
            NFecms2(NFnstep) = -ecms2                           ! CMS energy squared (negative: neutron)
          end if
          NFptot(NFnstep) = ptot                                ! total probability (to be removed)
        else
          print *,'nrprton.F: too many nucleon FSI steps, not saving all of them!'
        end if
        goto 200
      end if

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ia=abs(ides)-1         !number of pions produced
       do l=1,ia              !put pions on output stack
         if(no.lt.20) then
            no = no + 1
            if(chiso(2*l).ne.0) then !pions are on even indexes
               ido(no) = sign(211.01,chiso(2*l))
            else
               ido(no) = 111
            end if
            do l1=1,4
               po(l1,no)=piso(l1,2*l)
               if(l1.lt.4) stpo(l1,no)=xtr(l1)
            end do
         end if
      end do

c     continue with the first nucleon, put second on stack

      uin=amass(itdone)
      chint=chiso(1)
      ein = piso(4,1)
      pinss=0
      do j=1,3
         pin(j)=piso(j,1)
         pinss=pinss+pin(j)**2
      end do
      pinm=sqrt(pinss)
      do j=1,3
         upin(j)=pin(j)/pinm
      end do

      if(ntrack.lt.20) then   !put second nucleon on stack
         ntrack = ntrack+1
         chtrack(ntrack)=chiso(3)
         amass(ntrack) = amn
         do j=1,4
            pptrack(j,ntrack)=piso(j,3)
            if(j.lt.4) stpt(j,ntrack)=xtr(j)
         end do
      end if

!
!     fill nucleon fsi history stack
!
      if(NFnvert.lt.maxnucleonvert.and.NFnstep.lt.maxnucleonstep) then
        NFnvert = NFnvert+1                                   ! increase number of vertices
        NFiflag(NFnvert) = abs(ides) + chint*100              ! interaction type and charge of propagated nucleon
        if(prot) NFiflag(NFnvert) = NFiflag(NFnvert) + 10     ! charge of target nucleon
        NFx(NFnvert) = xtr(1)                                 ! vertex position inside nucleus
        NFy(NFnvert) = xtr(2)
        NFz(NFnvert) = xtr(3)
        NFpx(NFnvert) = pin(1)                                ! momentum
        NFpy(NFnvert) = pin(2)
        NFpz(NFnvert) = pin(3)
        NFe(NFnvert) = ein                                    ! energy
        NFfirststep(NFnvert) = NFnstep+1                      ! index of first step of track segment

        NFnstep = NFnstep+1                                   ! increase number of steps
        if(prot) then
          NFecms2(NFnstep) = ecms2                            ! CMS energy squared (positive: proton)
        else
          NFecms2(NFnstep) = -ecms2                           ! CMS energy squared (negative: neutron)
        end if
        NFptot(NFnstep) = ptot                                ! total probability (to be removed)
      else
        print *,'nrprton.F: too many nucleon FSI vertices (',NFnvert,') or steps (',NFnstep,'), not saving all of them!'
      end if

      goto 200              !continue stepping

c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 99   continue                  ! done with the current nucleon

!
!     fill nucleon fsi history stack
!
      if(abs(uin-amn).le.5.) then                               ! only fill fsi stack if dealing with a nucleon
        if(NFnvert.lt.maxnucleonvert.and.NFnstep.lt.maxnucleonstep) then
          NFnvert = NFnvert+1                                   ! increase number of vertices
          NFiflag(NFnvert) = 4 + chint*100                      ! interaction type and charge of propagated nucleon
          NFx(NFnvert) = xtr(1)                                 ! vertex position inside nucleus
          NFy(NFnvert) = xtr(2)
          NFz(NFnvert) = xtr(3)
          NFpx(NFnvert) = pin(1)                                ! momentum
          NFpy(NFnvert) = pin(2)
          NFpz(NFnvert) = pin(3)
          NFe(NFnvert) = ein                                    ! energy
          NFfirststep(NFnvert) = -1                             ! index of first step of track segment (no new segment starts here...)
        else
          print *,'nrprton.F: too many nucleon FSI vertices (',NFnvert,') or steps (',NFnstep,'), not saving all of them!'
        end if
      end if


C
C-- This will create NULL ( code=0 ) particles
C                           ( only the index is incremented )
C      no = no + 1
      if(no.lt.20) then
         no = no + 1
         do j=1,3               ! put him on output stack
            stpo(j,no)=xtr(j)
            po(j,no)=pin(j)
         end do
         po(4,no)=ein
         if(chint.eq.0) then
            ido(no) = 2112
         else
            ido(no) = 2212
         end if
      end if

      if(no.lt.20) goto 300    !-------^ take next from the stack

 100  continue
      
      return
      end
