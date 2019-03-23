c**********************************************************************
c     ------------------------------
      subroutine abspi(ivc)
c     ------------------------------
c
c     ( purpose )
c       Generate nucleons following Delta abbsorption
c
c     ( input )
c       via vcwork common block
c
c     ( output )
c       via vcwork common block
c
c     ( creation date and author )
c       2009.11.04 ; R. Tacik
C
C       2009.12.09 ; P. de Perio 
C                    - Added argument to accept slot of pion in vcwork
C                    - Added flags for outgoing nucleons
C
C       2010.06.10 ; P. de perio
C                    - Bug fix ipvc(5) -> ipvc(ivc)
c
c       2012.11.26 ; R. Tacik
c                    - initialize variable "pi"
c                    - Added new common block /labpion/ in absneject.h
c
c**********************************************************************

      implicit none

#include "vcwork.h"
#include "posinnuc.h"
#include "neutparams.h"

      real mn/0.94/
      real mmu/0.1056/
      real me/0.000511/
      real mpi/0.13957/
      real ml

      real be/0.008/

      real pi

      real enu,en,el,ed
      real pdx,pdy,pdz
      real s,tpi

      real rndm(12)

C      real EandP(4)
C      real p1cm(4),p2cm(4),p3cm(4),p4cm(4)
C      real p1lab(4),p2lab(4),p3lab(4),p4lab(4)

      real ppix,ppiy,ppiz,epi

      real p1,th1,ph1,p1x,p1y,p1z
      real p2,th2,ph2,p2x,p2y,p2z
      real p3,th3,ph3,p3x,p3y,p3z
      real p4,th4,ph4,p4x,p4y,p4z

      real qx,qy,qz,qe

      real pmu,ct,emu,rnu,q2
      real mp/0.93828/
      real mnv/0.93687/

      real p13,p23,p14,p34
      real p37,p47,p18,p12,p78

      integer absprob,numn,ivc

#include "absneject.h"

C     integer absprob,lpdg,npdg(4),numn,ivc
C      common /labEandP/ lpdg,EandP
C      common /cmnucleons/ p1cm,p2cm,p3cm,p4cm
C      common /labnucleons/ npdg,p1lab,p2lab,p3lab,p4lab

 
c  ivc'th particle got absorbed
      
C      write(*,*) 'PVC(i,ivc)=(', pvc(1,ivc), pvc(2,ivc), pvc(3,ivc),ivc, ')'

      pi = acos(-1.)

      ppix = pvc(1,ivc)/1000.
      ppiy = pvc(2,ivc)/1000.
      ppiz = pvc(3,ivc)/1000.
      epi  = sqrt(ppix**2+ppiy**2+ppiz**2+mpi**2)
      tpi  = epi - mpi

c  Fill elements of /labpion/ common block

      pilab(1) = epi
      pilab(2) = ppix
      pilab(3) = ppiy
      pilab(4) = ppiz


c  Determine 2-, 3-, and 4-body absorption
c  probabilities based on pion energy

      numn = absprob(tpi)

      if (numn.eq.3) go to 1000
      if (numn.eq.4) go to 2000

c  This section for 2-body absorption
c------------------------------------

c  Choose momenta for absorbing nucleons

      call ranlux(rndm,6)

      p1  = (rndm(1)*pfmax**3)**(1./3.)
      th1 = acos(1.-2.*rndm(2))
      ph1 = 2.*pi*rndm(3)
      p1x = p1*cos(ph1)*sin(th1)
      p1y = p1*sin(ph1)*sin(th1)
      p1z = p1*cos(th1)

      p2  = (rndm(4)*pfmax**3)**(1./3.)
      th2 = acos(1.-2.*rndm(5))
      ph2 = 2.*pi*rndm(6)
      p2x = p2*cos(ph2)*sin(th2)
      p2y = p2*sin(ph2)*sin(th2)
      p2z = p2*cos(th2)

c  Find total lab energy and momentum
c  Stored in 4-vector plab

      eandp(1) = epi + 2.*mn - 2.*be
      eandp(2) = ppix + p1x + p2x
      eandp(3) = ppiy + p1y + p2y
      eandp(4) = ppiz + p1z + p2z

      call abstwobody(tpi)

c  Decide of charge states of the outgoing nucleons

      call ranlux(rndm,1)

c     pi+
      if (ipvc(ivc).eq.211) then
         if (rndm(1).lt.0.70) then
            npdg(1) = 2212
            npdg(2) = 2212
         else
            npdg(1) = 2112
            npdg(2) = 2212
         end if

c     pi0
      else if (ipvc(ivc).eq.111) then
         npdg(1) = 2212
         npdg(2) = 2112


c     pi-
      else if (ipvc(ivc).eq.-211) then
         if (rndm(1).lt.0.70) then
            npdg(1) = 2112
            npdg(2) = 2112
         else
            npdg(1) = 2212
            npdg(2) = 2112
         end if
      end if

c  Fill vcwork
      if (nvc.GT.(maxvc-2)) then
         nvc = maxvc-2
      end if

      ipvc(nvc+1)  = npdg(1)
      icrnvc(nvc+1) = 1 
      ivtivc(nvc+1) = 1
      iorgvc(nvc+1) = ivc
      pvc(1,nvc+1) = p1lab(2)*1000.
      pvc(2,nvc+1) = p1lab(3)*1000.
      pvc(3,nvc+1) = p1lab(4)*1000.
      posnuc(1,nvc+1) = posnuc(1,ivc)
      posnuc(2,nvc+1) = posnuc(2,ivc)
      posnuc(3,nvc+1) = posnuc(3,ivc)
      iflgvc(nvc+1) = 0 
      CALL MCMASS(IPVC(NVC+1),AMASVC(NVC+1))

      ipvc(nvc+2)  = npdg(2)
      icrnvc(nvc+2) = 1
      ivtivc(nvc+2) = 1
      iorgvc(nvc+2) = ivc
      pvc(1,nvc+2) = p2lab(2)*1000.
      pvc(2,nvc+2) = p2lab(3)*1000.
      pvc(3,nvc+2) = p2lab(4)*1000.
      posnuc(1,nvc+2) = posnuc(1,ivc)
      posnuc(2,nvc+2) = posnuc(2,ivc)
      posnuc(3,nvc+2) = posnuc(3,ivc)
      iflgvc(nvc+2) = 0 
      CALL MCMASS(IPVC(NVC+2),AMASVC(NVC+2))

      nvc = nvc + 2

      return

c  This section for 3-body absorption
c------------------------------------

 1000 continue

c  Choose momentum for absorbing nucleons

      call ranlux(rndm,9)

      p1  = (rndm(1)*pfmax**3)**(1./3.)
      th1 = acos(1.-2.*rndm(2))
      ph1 = 2.*pi*rndm(3)
      p1x = p1*cos(ph1)*sin(th1)
      p1y = p1*sin(ph1)*sin(th1)
      p1z = p1*cos(th1)

      p2  = (rndm(4)*pfmax**3)**(1./3.)
      th2 = acos(1.-2.*rndm(5))
      ph2 = 2.*pi*rndm(6)
      p2x = p2*cos(ph2)*sin(th2)
      p2y = p2*sin(ph2)*sin(th2)
      p2z = p2*cos(th2)

      p3  = (rndm(7)*pfmax**3)**(1./3.)
      th3 = acos(1.-2.*rndm(8))
      ph3 = 2.*pi*rndm(9)
      p3x = p3*cos(ph3)*sin(th3)
      p3y = p3*sin(ph3)*sin(th3)
      p3z = p3*cos(th3)


c  Find total lab energy and momentum
c  Stored in 4-vector plab

      eandp(1) = epi + 3.*mn - 3.*be
      eandp(2) = ppix + p1x + p2x + p3x
      eandp(3) = ppiy + p1y + p2y + p3y
      eandp(4) = ppiz + p1z + p2z + p3z

      call absthreebody

c  Decide of charge states of the outgoing nucleons

      call ranlux(rndm,1)

c     pi+
      if (ipvc(ivc).eq.211) then
         if (rndm(1).lt.0.25) then
            npdg(1) = 2212
            npdg(2) = 2212
            npdg(3) = 2212
         else if (rndm(1).ge.0.85) then
            npdg(1) = 2212
            npdg(2) = 2112
            npdg(3) = 2112
         else
            npdg(1) = 2212
            npdg(2) = 2212
            npdg(3) = 2112
         end if

c     pi0
      else if (ipvc(ivc).eq.111) then
         if (rndm(1).lt.0.5) then
            npdg(1) = 2212
            npdg(2) = 2212
            npdg(3) = 2112
         else
            npdg(1) = 2212
            npdg(2) = 2112
            npdg(3) = 2112
         end if

c     pi-
      else if (ipvc(ivc).eq.-211) then 
         if (rndm(1).lt.0.25) then
            npdg(1) = 2112
            npdg(2) = 2112
            npdg(3) = 2112
         else if (rndm(1).ge.0.85) then
            npdg(1) = 2212
            npdg(2) = 2212
            npdg(3) = 2112
         else
            npdg(1) = 2212
            npdg(2) = 2112
            npdg(3) = 2112
         end if
      end if

c  Fill vcwork
      if (nvc.GT.(maxvc-3)) then
         nvc = maxvc - 3
      end if
      
      ipvc(nvc+1)  = npdg(1)
      icrnvc(nvc+1) = 1 
      ivtivc(nvc+1) = 1
      iorgvc(nvc+1) = ivc
      pvc(1,nvc+1) = p1lab(2)*1000.
      pvc(2,nvc+1) = p1lab(3)*1000.
      pvc(3,nvc+1) = p1lab(4)*1000.
      posnuc(1,nvc+1) = posnuc(1,ivc)
      posnuc(2,nvc+1) = posnuc(2,ivc)
      posnuc(3,nvc+1) = posnuc(3,ivc)
      iflgvc(nvc+1) = 0 
      CALL MCMASS(IPVC(NVC+1),AMASVC(NVC+1))

      ipvc(nvc+2)  = npdg(2)
      icrnvc(nvc+2) = 1 
      ivtivc(nvc+2) = 1
      iorgvc(nvc+2) = ivc
      pvc(1,nvc+2) = p2lab(2)*1000.
      pvc(2,nvc+2) = p2lab(3)*1000.
      pvc(3,nvc+2) = p2lab(4)*1000.
      posnuc(1,nvc+2) = posnuc(1,ivc)
      posnuc(2,nvc+2) = posnuc(2,ivc)
      posnuc(3,nvc+2) = posnuc(3,ivc)
      iflgvc(nvc+2) = 0 
      CALL MCMASS(IPVC(NVC+2),AMASVC(NVC+2))

      ipvc(nvc+3)  = npdg(3)
      icrnvc(nvc+3) = 1 
      ivtivc(nvc+3) = 1
      iorgvc(nvc+3) = ivc
      pvc(1,nvc+3) = p3lab(2)*1000.
      pvc(2,nvc+3) = p3lab(3)*1000.
      pvc(3,nvc+3) = p3lab(4)*1000.
      posnuc(1,nvc+3) = posnuc(1,ivc)
      posnuc(2,nvc+3) = posnuc(2,ivc)
      posnuc(3,nvc+3) = posnuc(3,ivc)
      iflgvc(nvc+3) = 0 
      CALL MCMASS(IPVC(NVC+3),AMASVC(NVC+3))

      nvc = nvc + 3

      return

c  This section for 4-body absorption
c------------------------------------

 2000 continue

c  Choose momentum for interacting nucleons

      call ranlux(rndm,12)

      p1  = (rndm(1)*pfmax**3)**(1./3.)
      th1 = acos(1.-2.*rndm(2))
      ph1 = 2.*pi*rndm(3)
      p1x = p1*cos(ph1)*sin(th1)
      p1y = p1*sin(ph1)*sin(th1)
      p1z = p1*cos(th1)

      p2  = (rndm(4)*pfmax**3)**(1./3.)
      th2 = acos(1.-2.*rndm(5))
      ph2 = 2.*pi*rndm(6)
      p2x = p2*cos(ph2)*sin(th2)
      p2y = p2*sin(ph2)*sin(th2)
      p2z = p2*cos(th2)

      p3  = (rndm(7)*pfmax**3)**(1./3.)
      th3 = acos(1.-2.*rndm(8))
      ph3 = 2.*pi*rndm(9)
      p3x = p3*cos(ph3)*sin(th3)
      p3y = p3*sin(ph3)*sin(th3)
      p3z = p3*cos(th3)

      p4  = (rndm(10)*pfmax**3)**(1./3.)
      th4 = acos(1.-2.*rndm(11))
      ph4 = 2.*pi*rndm(12)
      p4x = p4*cos(ph4)*sin(th4)
      p4y = p4*sin(ph4)*sin(th4)
      p4z = p4*cos(th4)

c  Find total lab energy and momentum
c  Stored in 4-vector plab

      eandp(1) = epi + 4.*mn - 4.*be
      eandp(2) = ppix + p1x + p2x + p3x + p4x
      eandp(3) = ppiy + p1y + p2y + p3y + p4y
      eandp(4) = ppiz + p1z + p2z + p3z + p4z

      call absfourbody

c  Decide of charge states of the outgoing nucleons

      call ranlux(rndm,1)

c     pi++
      if (ipvc(ivc).eq.211) then
         if (rndm(1).lt.0.05) then
            npdg(1) = 2212
            npdg(2) = 2212
            npdg(3) = 2212
            npdg(4) = 2212
         else if ((rndm(1).ge.0.05).and.(rndm(1).lt.0.55)) then
            npdg(1) = 2212
            npdg(2) = 2212
            npdg(3) = 2212
            npdg(4) = 2112
         else if ((rndm(1).ge.0.55).and.(rndm(1).lt.0.95)) then
            npdg(1) = 2212
            npdg(2) = 2212
            npdg(3) = 2112
            npdg(4) = 2112
         else
            npdg(1) = 2212
            npdg(2) = 2112
            npdg(3) = 2112
            npdg(4) = 2112
         end if

c     pi0
      else if (ipvc(ivc).eq.111) then
         npdg(1) = 2212
         npdg(2) = 2212
         npdg(3) = 2112
         npdg(4) = 2112

c     pi-
      else if (ipvc(ivc).eq.-211) then
         if (rndm(1).lt.0.05) then
            npdg(1) = 2112
            npdg(2) = 2112
            npdg(3) = 2112
            npdg(4) = 2112
         else if ((rndm(1).ge.0.05).and.(rndm(1).lt.0.55)) then
            npdg(1) = 2212
            npdg(2) = 2112
            npdg(3) = 2112
            npdg(4) = 2112
         else if ((rndm(1).ge.0.55).and.(rndm(1).lt.0.95)) then
            npdg(1) = 2212
            npdg(2) = 2212
            npdg(3) = 2112
            npdg(4) = 2112
         else
            npdg(1) = 2212
            npdg(2) = 2212
            npdg(3) = 2212
            npdg(4) = 2112
         end if
      end if

c  Fill vcwork
      if (nvc.GT.(maxvc-4)) then
         nvc = maxvc - 4
      end if

      ipvc(nvc+1)  = npdg(1)
      icrnvc(nvc+1) = 1 
      ivtivc(nvc+1) = 1
      iorgvc(nvc+1) = ivc
      pvc(1,nvc+1) = p1lab(2)*1000.
      pvc(2,nvc+1) = p1lab(3)*1000.
      pvc(3,nvc+1) = p1lab(4)*1000.
      posnuc(1,nvc+1) = posnuc(1,ivc)
      posnuc(2,nvc+1) = posnuc(2,ivc)
      posnuc(3,nvc+1) = posnuc(3,ivc)
      iflgvc(nvc+1) = 0 
      CALL MCMASS(IPVC(NVC+1),AMASVC(NVC+1))

      ipvc(nvc+2)  = npdg(2)
      icrnvc(nvc+2) = 1 
      ivtivc(nvc+2) = 1
      iorgvc(nvc+2) = ivc
      pvc(1,nvc+2) = p2lab(2)*1000.
      pvc(2,nvc+2) = p2lab(3)*1000.
      pvc(3,nvc+2) = p2lab(4)*1000.
      posnuc(1,nvc+2) = posnuc(1,ivc)
      posnuc(2,nvc+2) = posnuc(2,ivc)
      posnuc(3,nvc+2) = posnuc(3,ivc)
      iflgvc(nvc+2) = 0 
      CALL MCMASS(IPVC(NVC+2),AMASVC(NVC+2))

      ipvc(nvc+3)  = npdg(3)
      icrnvc(nvc+3) = 1 
      ivtivc(nvc+3) = 1
      iorgvc(nvc+3) = ivc
      pvc(1,nvc+3) = p3lab(2)*1000.
      pvc(2,nvc+3) = p3lab(3)*1000.
      pvc(3,nvc+3) = p3lab(4)*1000.
      posnuc(1,nvc+3) = posnuc(1,ivc)
      posnuc(2,nvc+3) = posnuc(2,ivc)
      posnuc(3,nvc+3) = posnuc(3,ivc)
      iflgvc(nvc+3) = 0 
      CALL MCMASS(IPVC(NVC+3),AMASVC(NVC+3))

      ipvc(nvc+4)  = npdg(4)
      icrnvc(nvc+4) = 1 
      ivtivc(nvc+4) = 1
      iorgvc(nvc+4) = ivc
      pvc(1,nvc+4) = p4lab(2)*1000.
      pvc(2,nvc+4) = p4lab(3)*1000.
      pvc(3,nvc+4) = p4lab(4)*1000.
      posnuc(1,nvc+4) = posnuc(1,ivc)
      posnuc(2,nvc+4) = posnuc(2,ivc)
      posnuc(3,nvc+4) = posnuc(3,ivc)
      iflgvc(nvc+4) = 0 
      CALL MCMASS(IPVC(NVC+4),AMASVC(NVC+4))

      nvc = nvc + 4

      return

      end
