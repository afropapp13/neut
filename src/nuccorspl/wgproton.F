      subroutine wgproton(itrack,xpos,anuc,nvc2)

c     interface routine to deal with inter-nuclear interactions of nucleon
c     written by W. Gajewski (3/25/2001)
c     It connects NEUT subroutine  NEFFECT  with W.G's programs NUCCOR
c     NUCCOR uses Bertini's database from MECC7 (1971) which is been
c     currently used by GCALOR (in the file chetc.dat)
c     itrack - track number in the VCWORK stack
c     xpos   - coordinates of the interaction point in nucleus (in fm)
c     anuc   - A of the nucleus (I do not know if it works for A.ne.16)
c     nvc2   - position of the last track in VCWORK stack

      implicit none

      integer itrack,nvc2,isuffix
      real*4  xpos(3)

#include "vcwork.h"
#include "vcvrtx.h"
#include "posinnuc.h"

      real*4 stpi(3),pi(4),anuc
C      real*4 po(4,20),stpo(3,20)
C      integer ido(20)
      real*4 po(4,MAXVC),stpo(3,MAXVC)
      integer ido(MAXVC)
      integer imode,icont,idmc,no,itpmax,itnmax
      real*4 amn,p4pmax,p4nmax,pcut
      parameter (amn=939.,pcut=amn+8.)
      data imode/2/

      integer i

      icont=1
      idmc = ipvc(itrack)
      call ucopy(xpos,stpi,3)
      call ucopy(pvc(1,itrack),pi,3)
      pi(4)=sqrt(amn**2+pi(1)**2+pi(2)**2+pi(3)**2)
      no=0

      call nrprton ( anuc,     !A of the nucleus =16 for Oxygen
     $     stpi,		!current point inside nucleus
     $     pi,idmc,	        !4momentum in GeV,id (should be p or n)
     $     po,ido,stpo,no,      !the same for outgoing, creat point,#of track
     $     imode,icont)

c     replace the incoming track with the most energetic proton (or neutron)
      p4pmax=0.
      p4nmax=0.
      itpmax=0
      itnmax=0
      do i=1,no
         if(ido(i).eq.2212.and.po(4,i).gt.p4pmax) then
            p4pmax=po(4,i)
            itpmax=i
         end if
         if(ido(i).eq.2112.and.po(4,i).gt.p4nmax) then
            p4nmax=po(4,i)
            itnmax=i
         end if
      end do

C      if(itpmax.ne.0) then     ! put fastest proton
C         call ucopy(po(1,itpmax),PVC(1,itrack),3)
C         ipvc(itrack)=ido(itpmax)
C         call ucopy(stpo(1,itpmax),POSFVC(1,itrack),3)
C      else if(itnmax.ne.0) then  !if no protons put fastest neutron
C         call ucopy(po(1,itnmax),PVC(1,itrack),3)
C         IPVC(itrack)=ido(itnmax)
C         call ucopy(stpo(1,itnmax),POSFVC(1,itrack),3)
C         itpmax=itnmax
C      else                        !if no nucleons absorb this track

      if (no.eq.1) then
         do 80 isuffix=1,3
            if (po(isuffix,1).ne.pvc(isuffix,itrack)) goto 85
 80      continue
C         write(*,*) "wg_proton:Nothing happened."
         return

 85      continue
      endif

      if ((itpmax.eq.0).and.(itnmax.eq.0)) then
C--- absorption
         IFLGVC(itrack)=3         ! I do not know what I am doing here !!!!!
         ICRNVC(itrack)=0
         IVTFVC(itrack)=1
      else
C--- interaction
         IFLGVC(itrack)=7
         ICRNVC(itrack)=0
         IVTFVC(itrack)=1
      end if
      
      if (no.eq.1) then
         if (ido(1).eq.ipvc(itrack)) then
            IFLGVC(itrack)=8
         else
            IFLGVC(itrack)=4
         endif
      endif

c     put the rest of outgoing tracks on stack (pions first)

      do i=1,no
         if(i.ne.itpmax .and. ido(i).lt.1000) then  !pions first
            IF (nvc2.ge.MAXVC-1) then
               write(*,*) 'wgproton: Too many particles..(1)'
               write(*,*) 'wgproton: Stop storing some of them.'
               goto 1000 
            endif
            nvc2 = nvc2 +1
            IPVC(NVC2)=ido(i)
            CALL UCOPY(PO(1,i),PVC(1,NVC2),3)
            IORGVC(NVC2)=itrack
            IFLGVC(NVC2)=0
            ICRNVC(NVC2)=1
C            CALL UCOPY(stpo(1,i),POSIVC(1,NVC2),3)
            CALL UCOPY(stpo(1,i),POSNUC(1,NVC2),3)
            IVTIVC(NVC2)=1
            CALL MCMASS(IPVC(NVC2),AMASVC(NVC2))
         end if
      end do

      do i=1,no
C         if(  i.ne.itpmax .and.         !track not used 
C     &        ido(i).gt.1000 .and.      !now nucleons
C     &        po(4,i).gt.pcut) then     !can exit nucleus ?
         if(  ido(i).gt.1000 .and.      !now nucleons
     &        po(4,i).gt.pcut) then     !can exit nucleus ?
            IF (nvc2.ge.MAXVC-1) then
               write(*,*) 'wgproton: Too many particles..(2)'
               write(*,*) 'wgproton: Stop storing some of them.'
               goto 1000 
            endif
            nvc2 = nvc2 +1
            IPVC(NVC2)=ido(i)           !track id
            CALL UCOPY(po(1,i),PVC(1,NVC2),3)   !momentum
            IORGVC(NVC2)=itrack                 !mother particle
            IFLGVC(NVC2)=0
            ICRNVC(NVC2)=1
C            CALL UCOPY(stpo(1,i),POSIVC(1,NVC2),3)  !production point
            CALL UCOPY(stpo(1,i),POSNUC(1,NVC2),3)  !production point
            IVTIVC(NVC2)=1
            CALL MCMASS(IPVC(NVC2),AMASVC(NVC2))
         end if
      end do
 1000 return
      end
