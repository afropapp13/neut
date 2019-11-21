      subroutine nrprton ( anuc, !A of the nucleus =16 for Oxygen
     $     stpi,                !current point inside nucleus
     $     pi,idmc,             !4momentum in GeV,id (should be p or n)
     $     po,ido,stpo,no,      !the same for outgoing, creat point,#of track
     $     imode,icont)

c     
c     ( Creation Date and Author )
c     2007.11.05 ; G.Mitsuka - add support other target than Oxygen
c     other target is treated by scaling to Oxygen
c     
c     2014.10.22 ; M.Hierholzer - add nucleon FSI history block to calculate systematics;
c     also add some documentation of the variables used in this file
c     
c     
c     icont =  0        generate interaction point
c     1        start from stpt  (set start(i)=0)
c     imode =  1        start on the surface
c     2        start in the volume  (probmi)
c     others   Dover mode       (probdi)
c     
c     When NEUT is run by neutgeom/genev, we always have icont=1 and imode=2.
c     
      
c     Arguments of nrprton:
c     stpi = input position
c     pi = input momentum
c     anuc = "A" value for nucleus
c     po = output momentum
c     stpo = output position
c     ido = output particle ID
c     imode = (used in nrstp.F) 1 = surface mode, 2 = volume mode, other = Dover mode(like volume mode but probdi vs probmi) 
c     icont = (used in nrstp.F/nrranve.F) 0 = generate interaction point, 1 = start from 0,0,0 
c     idmc = input paticle ID
c     no = number of output tracks

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
#include "fsinucleonreweight.h"
#include "necard.h"
!     #include "ntest.h"
!     Used external inputs from common blocks:
!     nrrmsrad from nrnuclparam.h - unknown meaning, depends on target, set in neutcore/nrsettarg.F
!     


!     
!     Output variables to common blocks:
!     chtrack(*) to nrstack.h - charge of nucleon in track
!     amass(*) to nrstack.h - mass of nucleon in track
!     stpt(3,*) to nrstack.h - starting point of track inside nucleus
!     pptrack(4,*) to nrstack.h - initial 4-momentum of nucleon in track
!     ntrack to nrstack.h - number of tracks
!     ihap to nrhad.h - ??
!     nel(3) to nrhad.h - ??
!     itrans to nrhadr.h - ??
!     

!     
!     Information passed on to (and back from) subroutines called in this routine via common blocks:
!     uin to nrint.h - mass of current track
!     chint to nrint.h - charge of current track
!     ein to nrint.h - energy of current track
!     pin(3) to nrint.h - momentum of current track
!     upin(3) to nrint.h - direction of current track (never actually used outside this function)
!     rhon to nrcms.h, set in nrfermi - nuclear density at current r
!     pfermi to nrcms.h, set in nrfermi - fermi momentum
!     unucl to nrcms.h , set in nrfermi- nucleon mass minus binding energy and fermi energy (?)
!     ecms2 to nrcms.h, set in nrfermi - center-of-mass energy squared (with fermi motion smearing?)
!     up4cms(3) to nrcms.h, set in nrfermi - direction of track (with fermi motion smearing?)
!     beta(3) to nrcms.h, set in nrfermi - unused?
!     prot to nrcms.h, set in nrfermi - decide if interaction is on proton (.true.) or neutron (.false.)
!     

!     
!     Purely internal variables:
!     rmax - nuclus radius, calculated from nrrmsrad in dependence of target
!     pptot - absolute initial nucleon momentum
!     start(3) - offset (?) of starting point for tracking, w.r.t. stpt()
!     itdone - track currently being processed
!     xtr(3) - current position (of current step) inside nucleus
!     r - current position's radius, calculated from xtr()
!     loop - some counter used in surface mode
!     iin - flag if inside nucleus, used in surface mode
!     

c     xtr - current position, used to find r, start - starting point with respect to stpt?
      real*4 xtr(3),start(3)

c     step - step size through nucleus, rmax - nucleus radius (from nrrmsrad)
c     amn - mass of nucleon, pptot - absolute initial nucleon momentum
      real step,rmax,amn,pptot
c     j - counter, itdone - track being processed
      integer j,itdone
c     pinss - squared sum of pin, pinm - root pinss, 
c     r - radius from xtr, az - square of rmax - sum of squares of x and y
      real*4 pinss,pinm,r,az
c     loop - counter in surface mode, iin - inside nucleus flag (surface mode),
c     ides - something to do with interaction type, ia = abs(ides)-1 number of pions produced
c     l1 - counter, l - counter
      integer loop,iin,ides,ia,l1,l
c     set step size and nucleon mass (constant - can't be changed further down)
      parameter (step=.2 )
c     not sure if 0.2 is optimal? probably is
      
C     integer*4 nclux,ncint,nck1,nck2

c     rmsradius for oxygen 
      REAL*4 RMSRADOXY
      parameter (RMSRADOXY = 2.730)
c     something for fermi gas
      REAL effrmgas


C     integer RW_step_i

c     ptot - total probability, ABC seemingly unused
      real*4 ptot, A,B,C
      common /martindebug/  ptot,A,B,C

c      pcascprob = 1
C      pelnuc =1
C      pspnuc =1 
C      pdpnuc =1
c     initialise before loop
c      write (*,*) "nfnstep start of nrprton = " , NFNStep

c      write (*,*) "Start of nrprton"
c      write (*,*) "NReWeightFlag = ", NReweightFlag
c      write (*,*) "Step numberA nfnstep ", NFnstep

      if(NReweightFlag.eq.1) then

c     NFnstep = 0
         RW_step_i =0

c         write (*,*) "RW StepA  ", RW_step_i
c         write (*,*) "Step numberA nfnstep ", NFnstep
         nucresflg=1
         

c     NFnstep = 0
         goto 200
      endif

C     are we reweighting? flag 1 if we are! (this might need to go higher up)                                                                                       C     I think we need to this right at the beginning - might be possible to kove it a bit further down later 
      if(NReweightFlag.eq.0) then
C     if not, then set fsiNFnucresflg from nucresflg                                                                                                                 
         NFnucresflg = 0;
c     if so, then other way around                                                                                                                                   
      else if(NReweightFlag.eq.1) then
c     ie. nucresflg should equal 1

         nucresflg=1;
      endif
c     this is a bit silly, just keeping it here for now





C     rmax is set by scaling to 16 O (rmax=6. for Oxygen)
c     nrrmsrad from nrnuclparam
c     For Reweight: Should be ok?
      rmax = 6.*(NRRMSRAD/RMSRADOXY)





c     For Reweight: idmc came from inputs, so should be set ok
c     if it's a proton set charge 1
      if(idmc.eq.2212) then     !put the input track on stack
         chtrack(1)=1
         amn = 938.28
c     or if it's a neutron set charge 0
      else if(idmc.eq.2112) then
         chtrack(1)=0
         amn = 939.57
c     otherwise have an existential crisis (and return)
      else
         print *,'*** NUCPRTON ERROR call for wrong track id=',idmc
         return
      end if

c     For Reweight: amn const
c     initialise mass to amn (set above in parameter), set r to 0 and pptot to 0, these will change soon
      amass(1) = amn
      pptot=0
      r=0
c     from 1 to 3
c     For Reweight: these come from input so should be ok
      do j=1,3
c     set starting point of this track (1) to input starting points
         stpt(j,1)=stpi(j)
c     calculate square of radius
         r=r+stpi(j)**2
      end do
c     then square root
      r=sqrt(r)
c     from 1 to 4
      do j=1,4
c     set initial nucleon 4 momentum of this track (1) from input momentum 
         pptrack(j,1)=pi(j)
c     for 3 momenta, take square sum and sqrt it below
         if(j.lt.4) pptot = pptot + pptrack(j,1)**2
      end do
      pptot = sqrt(pptot)
      

C     Asmita R.:For MEC events fermi surface momentum is function of radial distance r (22-07-2013)
C     while for other modes it's read from card
C     Change to local Fermi Gas ( 2017/04/03 )


c     not sure what nucresflg 1 vs 2 means
c     if  the nucres flag is 1, let's go:
      if (nucresflg.eq.1) then
c     modene from nework, neutmode
         if (abs(MODENE).EQ.2) then
c     can't go beyond rmax, set fermi momentum
            if (r.ge.rmax) then
               pfermi=EFFRMGAS(0,0,rmax)
            else
               pfermi=EFFRMGAS(0,0,r)
            endif
c     not sure what nrpfsurf is, some alternative to effrmgas
         else 
            pfermi=NRPFSURF*1000. ! NRPFSURF is GeV units
         endif
c     if nucresflg is 2
      else if (nucresflg.eq.2) then
         if (r.ge.rmax) then
            pfermi=EFFRMGAS(0,0,rmax)
         else
            pfermi=EFFRMGAS(0,0,r)
         endif
c     so nucresflg 1 & 2 treated the same, except for 1 check if modene is not 2
c     if not 1 or 2, kick off
      else
         write(*,*) 'Unknown NUCRESFLG ( expected 1 or 2 )',NUCRESFLG
         STOP
      endif

c     now we've set pfermi
c     is our initial total momentum less than fermi momentum?
      if(pptot.lt.pfermi) then  !if it is below fermi sphere get out
c     if so, set ouput momentum/position to input momentum/position
c     ido set to nucleon ID and number of tracks is 1. QUIT
         do j=1,4
            po(j,1)=pi(j)
            if(j.lt.4) stpo(j,1)=stpi(j)
         end do
         ido(1)=idmc
         no=1
         return
      end if

c     icont and imode are inputs
C     call nrstp which fixes the starting point
      call nrstp(start,icont,imode) !fix the starting point
      
c     currently on track 0
      itdone = 0                !current track
c     number of tracks left 1
      ntrack = 1                !number of tracks to go
c     number of output tracks
      no = 0                    !number of output tracks
c     itrans iterates with every new track (itdone) 
      itrans = 0
c!!!!!OD!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c     loop over input tracks
      NFnstep = NFnstep+1       ! increase number of steps                                                                                                                                 
c      write (*,*) "nfstep increased line 308", NFnstep

c      write (*,*) "before 300  "
 300  continue                  ! loop thru the input tracks
c      write (*,*) "after 300  "
c     loop counter itdone
      itdone=itdone+1
c     if the counter goes beyond ntracks, QUIT
      if(itdone.gt.ntrack) goto 100 !just return
c     if this is the first track, reset nel and ihap interaction marker
      if(itdone.le.1)then       ! reset interaction markers
         ihap=-1
         nel(1)=0
         nel(2)=0
         nel(3)=0
      end if
C     write (*,*) "RW Step B  ", RW_step_i
C     write (*,*) "Step numberB nfnstep ", NFnstep

c     For Reweight: amass/chtrack/stpt should be set in nrstack.h common block (check this)
C     For Reweight: pptrack set from pi (input) above
c     set uin to mass of current track
      uin=amass(itdone)         !put a new track to follow
c     set chint to charge of current track
      chint=chtrack(itdone)
c     set ein to energy of current track
      ein = pptrack(4,itdone)
      pinss=0
c     set xtr/pinss to position and momentum squared of current track
c     For NOT Reweight: only want to set position and momentum this way if not reweighting:
      if(NReweightFlag.eq.0) then
         do j=1,3
            if(itdone.gt.1)start(j)=0.
            xtr(j)=stpt(j,itdone)+start(j)
            pin(j)=pptrack(j,itdone)
            pinss=pinss+pin(j)**2
         end do

c     For Reweight: if we are reweighting, set x,mom, and ein from common block
C     else  if(NReweightFlag.eq.1) then
C     NFnvert = NFnvert+1                                   ! increase number of vertices                                                                      
c     Set NFiflag to 100*chint (presumably 100 for proton 000 for pion)                                                                                            
C     NFiflag(NFnvert) = 0 + chint*100                      ! interaction type and charge of propagated nucleon                                                 
         
c     set xyx position, momentum, and then energy                                                                                                                   
C     NFx(NFnvert) = xtr(1)  ! vertex position inside nucleus                                                                       
C     NFy(NFnvert) = xtr(2)
C     NFz(NFnvert) = xtr(3)
C     NFpx(NFnvert) = pin(1) ! momentum                                                                                             
C     NFpy(NFnvert) = pin(2)
C     NFpz(NFnvert) = pin(3)
C     NFe(NFnvert) = ein     ! energy                                                                                               
c     set first step of track to nfnstep+1                                                                                                                           
C     NFfirststep(NFnvert) = NFnstep+1 ! index of first step of track segment                                                                 


C     NFnstep = NFnstep+1                                   ! increase number of steps                                                                             
c     set cms2 to 0 because it's not a collision                                                                                                                     
C     NFecms2(NFnstep) = 0.                                 ! CMS energy squared (not a collision, so it's zero)                                                   
c     conservation of momentum when there's no interaction                                                                                                           
C     NFptot(NFnstep) = 0.
      endif
      pinm=sqrt(pinss)
c     normalise momentum for this track
      do j=1,3
         upin(j)=pin(j)/pinm
      end do
c     if mass this track - mass of first track is greater than 5, done for this nucleon
c     check for pions mass
      if(abs(uin-amn).gt.5.) goto 99 !not a nucleon

c     initialsie loop and iin
      loop=0
      iin=0

!     
!     fill nucleon fsi history stack
!     
c     still outside event loop, but if NFnvert and NFnstep are below maximum, start filling common
      if(NFnvert.lt.maxnucleonvert.and.NFnstep.lt.maxnucleonstep) then
c     increase NFNvert and NFNstep by 1
         NFnvert = NFnvert+1    ! increase number of vertices
c     write (*,*) "nfnvert1 is here to stay", NFnvert

c     Set NFiflag to 100*chint (presumably 100 for proton 000 for pion)
         NFiflag(NFnvert) = 0 + chint*100 ! interaction type and charge of propagated nucleon
c     set xyx position, momentum, and then energy
         NFx(NFnvert) = xtr(1)  ! vertex position inside nucleus
         NFy(NFnvert) = xtr(2)
         NFz(NFnvert) = xtr(3)
         NFpx(NFnvert) = pin(1) ! momentum
         NFpy(NFnvert) = pin(2)
         NFpz(NFnvert) = pin(3)
         NFe(NFnvert) = ein     ! energy
c     set first step of track to nfnstep+1
         NFfirststep(NFnvert) = NFnstep ! index of first step of track segment

c         NFnstep = NFnstep+1    ! increase number of steps
c         write (*,*) "nfstep increased line 401", NFnstep

c     write (*,*) "nfnstep1 is here to stay", NFnstep
c     write (*,*) "NfiFlagStep = ", NfiFlagStep(NFnstep)


         NFxStep(NFnstep) = xtr(1)   
         NFyStep(NFnstep) = xtr(2)
         NFzStep(NFnstep) = xtr(3)
         NFpxStep(NFnstep) = pin(1)  
         NFpyStep(NFnstep) = pin(2)
         NFpzStep(NFnstep) = pin(3)
         NFeStep(NFnstep) = ein 
C     NFiflagStep(NFnstep)=ides
C     write (*,*) "NFiFlagStep (generated)3= ", NFiflagStep(NFnstep)

c     set cms2 to 0 because it's not a collision
         NFecms2(NFnstep) = 0.  ! CMS energy squared (not a collision, so it's zero)
c     conservation of momentum when there's no interaction
         NFptot(NFnstep) = 0.   ! total probability (to be removed)
      else
c     if we're over the limit of maxnucleonvert/step, kick off and throw a message
         print *,'nrprton.F: too many nucleon FSI vertices (',NFnvert,')
     &or steps (',NFnstep,'), not saving all of them!'
      end if

c     write (*,*) "RW StepC  ", RW_step_i
c     write (*,*) "Step numberC nfnstep ", NFnstep


!     
!     track nuclon through nucleus
!     
      RW_step_i =0
c     Event loop time. Loop over steps through the nucleus
c      write (*,*) "before 200  "
 200  continue
c      write (*,*) "after 200  "
c     set r from xtr and step*direction. Does upin stay the same throughout?
      if(NReweightFlag.eq.1) then
         
         if(RW_step_i.lt.NFnstep-1) then
            RW_step_i = RW_step_i+1
            


            xtr(1) = NFxStep(RW_step_i)
            xtr(2) = NFyStep(RW_step_i)
            xtr(3) = NFzStep(RW_step_i)
            pin(1) = NFpxStep(RW_step_i)
            pin(2) = NFpyStep(RW_step_i)
            pin(3) = NFpzStep(RW_step_i)
            ein = NFeStep(RW_step_i)
c     possibly needs a +/- distinction for protons and pions                                                                                                         
            ecms2 = NFecms2(RW_step_i)
c     ides = NFiflagstep(RW_step_i)
            
C            do j=1,3
C               if(itdone.gt.1)start(j)=0.
C               xtr(j)=stpt(j,itdone)+start(j)
C               pin(j)=pptrack(j,itdone)
C               pinss=pinss+pin(j)**2
C            end do


C            pinm=sqrt(pinss)
c     normalise momentum for this track                                                                                                                                                      
C            do j=1,3
C               upin(j)=pin(j)/pinm
C            end do
      
c     write (*,*) "RW Step  ", RW_step_i
c     write (*,*) "Step number nfnstep ", NFnstep
c     write (*,*) "xtr1 = ", xtr(1)
c     write (*,*) "xtr2 = ", xtr(2)
c     write (*,*) "xtr3 = ", xtr(3)
c     write (*,*) "pin1 = ", pin(1)
c     write (*,*) "pin2 = ", pin(2)
c     write (*,*) "pin3 = ", pin(3)
c     write (*,*) "ein = ", ein
c     write (*,*) "ecms2 = ", ecms2
c     write (*,*) "IDES = ", ides



         else
c            write (*,*) "finished number of steps for event", NFnstep
            goto 300
         end if
         
      end if





      r=0
      do j=1,3
         xtr(j)=xtr(j)+upin(j)*step
         r=r+xtr(j)**2
c         write (*,*) "xtr(j) = ", xtr(j)   
c         write (*,*) "upin(j) = ", upin(j)
c         write (*,*) "step = ", step
      end do
      r=sqrt(r)



      if(NReweightFlag.eq.0) then


c     if surface mode
         if(imode.eq.1)then
c     and if in nucleus flag=0 and r bigger than rmax
            if(iin.eq.0.and.r.gt.rmax)then
c     count on loop, and if bigger than 3, we're done with this nucleon 
c     ie, if we're outside the nucleus 3 times in a row, we're done?
               loop=loop+1
               if(loop.gt.3) goto 99
c     calculate if position is within rmax
               az=rmax**2-xtr(1)**2-xtr(2)**2
c     if we're outside the circle, then we're done with this nucleon
               if(az.lt.0.) goto 99
c     each time we go along the loop, make xtr slightly bigger?
c     not sure why it's 0.18 or 0.05 (numerical precision?)
               xtr(3)=-sqrt(az)-0.18+0.05*loop
c     go to next step
               goto 300
            end if
c     if r is less than rmax iin=1 back in nucleon
            if(r.le.rmax)iin=1
         end if


         if(r.gt.rmax) then     !particle got out of the nucleus
c     uf the track number is still less than the total number of tracks
            if(itdone.le.ntrack) then !count #of interactions
c     coubt itrans
               itrans=itrans+1
C     count inot if ihap(?) less than 0 NB counted negative in contrast to nel(j) 
               if(ihap.lt.0)inot=inot+1
c     count iel if nel(1) greater than 0
               if(nel(1).gt.0)iel=iel+nel(1)
c     count ine1 if nel(2) greater than 0
               if(nel(2).gt.0)ine1=ine1+nel(2)
c     count ine2 if nel(3) greater than 0
               if(nel(3).gt.0)ine2=ine2+nel(3)
c     inot, iel, ine1, ine2 must be different interaction types. Need to sort this out
c     pretty sure nel(1/2/3)= pel/psp/pdp. maybe ihap is survival?
            end if
            goto 99             ! ------>  put on output stack
         end if
      end if

C     if(NReweightFlag.eq.1) then
C     xtr(1) = NFxStep(RW_step_i)
C     xtr(2) = NFyStep(RW_step_i) 
C     xtr(3) = NFzStep(RW_step_i) 
C     pin(1) = NFpxStep(RW_step_i) 
C     pin(2) = NFpyStep(RW_step_i)
C     pin(3) = NFpzStep(RW_step_i)
C     ein = NFeStep(RW_step_i)
c     possibly needs a +/- distinction for protons and pions
C     ecms2 = NFecms2(RW_step_i) 
      
C     end if
c     write (*,*) "RW Step  ", RW_step_i
c     write (*,*) "Step number nfnstep ", NFnstep
c     write (*,*) "xtr1 = ", xtr(1)
c     write (*,*) "xtr2 = ", xtr(2)
c     write (*,*) "xtr3 = ", xtr(3)
c     write (*,*) "pin1 = ", pin(1)
c     write (*,*) "pin2 = ", pin(2)
c     write (*,*) "pin3 = ", pin(3)
c     write (*,*) "ein = ", ein
c     write (*,*) "ecms2 = ", ecms2          
c     write (*,*) "r = ", r
C     if(NReweightFlag.eq.1) then
C     write (*,*) "NFIflag = ", ides
C     endif
c      write (*,*) "r before nfermi = ", r

      call nrfermi(r)           !set target nucleon
C     if(ecms2.lt.(2*amn)**2) goto 200  ! kinematically impossible / test is done inside nrnuc!

c     NFiflagStep(NFnstep) = 99999
      if(NReweightFlag.eq.1) then
c     write (*,*) "I HAVE MADE IT THIS FAR HOORAY"
         ides = NFiflagStep(RW_step_i)
      endif
c     write (*,*) "amn before nrnc = ", amn
c     write (*,*) "itdone number before nrnuc = ", itdone    
c     write (*,*) "no number before nrnuc = ", no
c     write (*,*) "track number before nrnuc = ", ntrack

      call nrnuc(ides)          ! make interaction
c     if ides was set to 11,12,13 then it was pauli blocked
c      write (*,*) "nfnstep coming out of nrnuc = ", NFnStep

      
c     write (*,*) "Interaction flag =" , ecms2
c     write (*,*) "NFiflagStep after nrnuc" , NFiflagStep(NFnstep)
      if(NReweightFlag.eq.1) then
         
         goto 200
      endif
      if(NReweightFlag.eq.0) then
         if(ides.ge.10) then    ! Paulu blocked interaction: save to fsi history stack
            if(NFnvert.lt.maxnucleonvert.and.NFnstep.lt.maxnucleonstep) then
c     assuming we're still below max vert and step
               NFnvert = NFnvert+1 ! increase number of vertices
c               write (*,*) "nfnvert2 is here to stay", NFnvert
               NFiflag(NFnvert) = ides-10 + chint*100 + 1000 ! Pauli-blocked interaction type and charge of propagated nucleon
c     previously we set NFiFlag to chint*100.
c     now we need to set 1000 for Pauli times
C     and we set ides-1 so we know what would have happened
               if(prot) NFiflag(NFnvert) = NFiflag(NFnvert) + 10 ! charge of target nucleon
c     and we add 10 if our target was a proton
c     now set position and momentum as per usual
               NFx(NFnvert) = xtr(1) ! vertex position inside nucleus
               NFy(NFnvert) = xtr(2)
               NFz(NFnvert) = xtr(3)
               NFpx(NFnvert) = pin(1) ! momentum
               NFpy(NFnvert) = pin(2)
               NFpz(NFnvert) = pin(3)
               NFe(NFnvert) = ein ! energy
               NFfirststep(NFnvert) = NFnstep+1 ! index of first step of track segment
               
               NFnstep = NFnstep+1 ! increase number of steps
c               write (*,*) "nfstep increased line 606", NFnstep

c     write (*,*) "nfnstep2 is here to stay", NFnstep

c     write (*,*) "NfiFlagStep = ", NfiFlagStep(NFnstep)

c     NFiflagStep(NFnstep)=ides
c     write (*,*) "NFiFlagStep (generated) = ", NFiflagStep(NFnstep)
               NFxStep(NFnstep) = xtr(1)
               NFyStep(NFnstep) = xtr(2)
               NFzStep(NFnstep) = xtr(3)
               NFpxStep(NFnstep) = pin(1)
               NFpyStep(NFnstep) = pin(2)
               NFpzStep(NFnstep) = pin(3)
               NFeStep(NFnstep) = ein
               
               

c     why is cms^2 positive for protons, and negative for pions
               if(prot) then
                  NFecms2(NFnstep) = ecms2 ! CMS energy squared (positive: proton)
               else
                  NFecms2(NFnstep) = -ecms2 ! CMS energy squared (negative: neutron)
               end if
c     I think this is just something Martin did for debugging 
               NFptot(NFnstep) = ptot ! total probability (to be removed)
c     kick off if we've gone over max vert and step
            else
               print *,'nrprton.F: too many nucleon FSI vertices (',NFnvert,')
     &or steps (',NFnstep,'), not saving all of them!'
            end if
c     head to next step
            goto 200            ! do the next step
         end if

c     if it wasn't Pauli blocked
c     ides was set to -1,-2,-3 for interactions, or 1 for survival

         if(ides.gt.0) then     ! no interaction: do next step
c     therefore if ides is +ve, there's no interaction, head to next step
            if(NFnstep.lt.maxnucleonstep) then ! save per-step information
               NFnstep = NFnstep+1 ! increase number of steps
c               write (*,*) "nfstep increased line 648", NFnstep

c     write (*,*) "nfnstep3 is here to stay", NFnstep

c     write (*,*) "NfiFlagStep = ", NfiFlagStep(NFnstep)

               NFxStep(NFnstep) = xtr(1)
               NFyStep(NFnstep) = xtr(2)
               NFzStep(NFnstep) = xtr(3)
               NFpxStep(NFnstep) = pin(1)
               NFpyStep(NFnstep) = pin(2)
               NFpzStep(NFnstep) = pin(3)
               NFeStep(NFnstep) = ein
c     NFiflagStep(NFnstep)=ides
c     write (*,*) "NFiFlagStep (generated2) = ", NFiflagStep(NFnstep)

               if(prot) then
                  NFecms2(NFnstep) = ecms2 ! CMS energy squared (positive: proton)
               else
                  NFecms2(NFnstep) = -ecms2 ! CMS energy squared (negative: neutron)
               end if
               NFptot(NFnstep) = ptot ! total probability (to be removed)
            else
               print *,'nrprton.F: too many
     &nucleon FSI steps, not saving all of them!'
            end if
            goto 200
         end if

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c     otherwise, we're not pauli blocked and we've got an interaction!
c     how many pions were produced? =ia
         ia=abs(ides)-1         !number of pions produced
c     for each pion...
         do l=1,ia              !put pions on output stack
c     for each one in the set of tracks (up to a maximum of 20)
            if(no.lt.20) then
c     increment track counter
               no = no + 1
c     for each charge of pions (even indeces apparently)
               if(chiso(2*l).ne.0) then !pions are on even indexes
c     set id for this pion      
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

         if(ntrack.lt.20) then  !put second nucleon on stack
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
            NFnvert = NFnvert+1 ! increase number of vertices
c     write (*,*) "nfnvert3 is here to stay", NFnvert

            NFiflag(NFnvert) = abs(ides) + chint*100 ! interaction type and charge of propagated nucleon
            if(prot) NFiflag(NFnvert) = NFiflag(NFnvert) + 10 ! charge of target nucleon
            NFx(NFnvert) = xtr(1) ! vertex position inside nucleus
            NFy(NFnvert) = xtr(2)
            NFz(NFnvert) = xtr(3)
            NFpx(NFnvert) = pin(1) ! momentum
            NFpy(NFnvert) = pin(2)
            NFpz(NFnvert) = pin(3)
            NFe(NFnvert) = ein  ! energy
            NFfirststep(NFnvert) = NFnstep+1 ! index of first step of track segment

            NFnstep = NFnstep+1 ! increase number of steps
c            write (*,*) "nfstep increased line 748", NFnstep

c     write (*,*) "nfnstep4 is here to stay", NFnstep

c     write (*,*) "NfiFlagStep = ", NfiFlagStep(NFnstep)

            NFxStep(NFnstep) = xtr(1)
            NFyStep(NFnstep) = xtr(2)
            NFzStep(NFnstep) = xtr(3)
            NFpxStep(NFnstep) = pin(1)
            NFpyStep(NFnstep) = pin(2)
            NFpzStep(NFnstep) = pin(3)
            NFeStep(NFnstep) = ein
c     NFiflagStep(NFnstep)=ides
c     write (*,*) "NFiFlagStep (generated) = ", NFiflagStep(NFnstep)


            if(prot) then
               NFecms2(NFnstep) = ecms2 ! CMS energy squared (positive: proton)
            else
               NFecms2(NFnstep) = -ecms2 ! CMS energy squared (negative: neutron)
            end if
            NFptot(NFnstep) = ptot ! total probability (to be removed)
         else
c     print *,'nrprton.F: too many nucleon FSI vertices 
c     & (',NFnvert,') or steps (',NFnstep,'), not saving all of them!'
         end if
      endif
      goto 200                  !continue stepping

c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 99   continue                  ! done with the current nucleon

!     
!     fill nucleon fsi history stack
!     
      if(abs(uin-amn).le.5.) then ! only fill fsi stack if dealing with a nucleon
         if(NFnvert.lt.maxnucleonvert.and.NFnstep.lt.maxnucleonstep) then
            NFnvert = NFnvert+1 ! increase number of vertices
c     write (*,*) "nfnvert4 is here to stay", NFnvert
            NFiflag(NFnvert) = 4 + chint*100 ! interaction type and charge of propagated nucleon
            NFx(NFnvert) = xtr(1) ! vertex position inside nucleus
            NFy(NFnvert) = xtr(2)
            NFz(NFnvert) = xtr(3)
            NFpx(NFnvert) = pin(1) ! momentum
            NFpy(NFnvert) = pin(2)
            NFpz(NFnvert) = pin(3)
            NFe(NFnvert) = ein  ! energy
            NFfirststep(NFnvert) = -1 ! index of first step of track segment (no new segment starts here...)
         else
c     print *,'nrprton.F: too many nucleon FSI vertices 
c     & (',NFnvert,') or steps (',NFnstep,'), not saving all of them!'
         end if
      end if


C     
C--   This will create NULL ( code=0 ) particles
C     ( only the index is incremented )
C     no = no + 1
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
c         write (*,*) "heading back to 300  "
      end if

      if(no.lt.20) goto 300     !-------^ take next from the stack

 100  continue

      return
      end
