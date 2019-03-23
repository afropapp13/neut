      subroutine nrprton ( anuc, !A of the nucleus =16 for Oxygen
     $     stpi,                !current point inside nucleus
     $     pi,idmc,             !4momentum in GeV,id (should be p or n)
     $     po,ido,stpo,no,      !the same for outgoing, creat point,#of track
     $     imode,icont)

*
*     ( Creation Date and Author )
*       2007.11.05 ; G.Mitsuka - add support other target than 16O
*                                other target is treated by scaling
*                                to 16O
*

      implicit none
c     icont =  0  generate interaction point
c     1  start from stpt  (set start(i)=0)
c     imode  = 1  start on the surface
c     2  start in the volume  (probmi)
c     others Dover mode       (probdi)
c     anuc   = 16.

      real*4 stpi(3),pi(4),anuc
      real*4 po(4,20),stpo(3,20)
      integer ido(20)
      integer imode,icont,idmc,no

#include "vcwork.h"
#include "nework.h"
#include "nrnuclparam.h"
#include "nrint.h"
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
#include "neutparams.h"
      real*4 xtr(3),start(3)

      real step,rmax,amn,pptot
      integer iraz,j,itdone
      real*4 pinss,pinm,r,az
      integer loop,iin,ides,ia,l1,l

      parameter (step=.2, amn=939. )
      data iraz/0/

      integer*4 nclux,ncint,nck1,nck2

      REAL*4 RMSRADOXY
      parameter (RMSRADOXY = 2.730)
      REAL effrmgas


C     rmax is set by scaling to 16O(rmax=6. for 16O)
      rmax = 6.*(NRRMSRAD/RMSRADOXY)
      if(iraz.eq.0)then
         iraz=1
         call nrinit(anuc,icont)      !set nuclear matter distribution
      end if

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
        pfermi= NRPFSURF*1000. ! NRPFSURF is GeV uniti

       if(abs(MODENE).EQ.2) then 
 	pfermi=EFFRMGAS(0,0,r)
       endif
        print*, 'MODE = ',MODENE, pfermi

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
      if(itdone.gt.ntrack)go to 100 !just return
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
      if(abs(uin-amn).gt.5.) go to 99 !not a nucleon

      loop=0
      iin=0

 200  continue                  !****start tracking***********

      r=0
      do j=1,3
         xtr(j)=xtr(j)+upin(j)*step
         r=r+xtr(j)**2
      end do
      r=sqrt(r)

      if(imode.eq.1)then
         if(iin.eq.0.and.r.gt.rmax)then
            loop=loop+1
            if(loop.gt.3)go to 99
            az=rmax**2-xtr(1)**2-xtr(2)**2
            if(az.lt.0.)go to 99
            xtr(3)=-sqrt(az)-0.18+0.05*loop
            go to 200
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
         go to 99               ! ------>  put on output stack
      end if

      call RLUXAT(nclux,ncint,nck1,nck2)
C      write(*,*) 'RLUX before nucfermi: (nclux,ncint,nck1,nck2)=(',
C     $     nclux,',',ncint,',',nck1,',',nck2,')'

      call nrfermi(r)           !set target nucleon
      call RLUXAT(nclux,ncint,nck1,nck2)
C      write(*,*) 'RLUX after nucfermi : (nclux,ncint,nck1,nck2)=(',
C     $     nclux,',',ncint,',',nck1,',',nck2,')'

      if(ecms2.lt.(2*amn)**2) go to 200
      ides=itdone*10
      call nrnuc(ides)          !make interaction

      if(ides.eq.1)go to 200    !no interaction make next step

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(ides.lt.0)then         !there was an interaction
         ia=abs(ides)-1         !number of pions produced
         do l=1,ia              !put pions on output stack
C
C-- This will create NULL ( code=0 ) particles
C                           ( only the index is incremented )
C           no=no+1
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

C
C-- This will create NULL ( code=0 ) particles
C                           ( only the index is incremented )
C         ntrack = ntrack +1     !second nucleon on stack
         if(ntrack.lt.20) then
            ntrack = ntrack +1   !second nucleon on stack
            chtrack(ntrack)=chiso(3)
            amass(ntrack) = amn
            do j=1,4
               pptrack(j,ntrack)=piso(j,3)
               if(j.lt.4) stpt(j,ntrack)=xtr(j)
            end do
         end if

         go to 200              !continue stepping
      end if

c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 99   continue                  ! done with the current nucleon
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

      if(no.lt.20) go to 300    !-------^ take next from the stack

 100  continue
      
      return
      end
