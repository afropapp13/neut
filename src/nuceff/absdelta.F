c**********************************************************************
c     ------------------------------
      subroutine absdelta
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
C                    - Added flags for outgoing nucleons
c
c       2012.11.26 ; R. Tacik
c                    - Added new common block /labpion/ in absneject.h
c
c**********************************************************************

      implicit none

#include "vcwork.h"
#include "neutparams.h"
#include "posinnuc.h"

      real mn/0.94/
      real mmu/0.1056/
      real me/0.000511/
      real mtau/1.77684/
      real mpi/0.13957/
      real ml

      real be/0.008/

      real pi

      real enu,en,el,ed
      real pdx,pdy,pdz
      real s,tpi

      real rndm(9)

C      real EandP(4)
C      real p1cm(4),p2cm(4),p3cm(4),p4cm(4)
C      real p1lab(4),p2lab(4),p3lab(4),p4lab(4)

      real p1,th1,ph1,p1x,p1y,p1z
      real p2,th2,ph2,p2x,p2y,p2z
      real p3,th3,ph3,p3x,p3y,p3z

      real dummy
      real qx,qy,qz,qe

      real pmu,ct,emu,rnu,q2
      real mp/0.93828/
      real mnv/0.93687/

      real p13,p23,p14,p34
      real p37,p47,p18,p12,p78

      real pd,epi,ppi

C      integer absprob,lpdg,npdg(4),numn
      integer absprob,numn

#include "absneject.h"
C      common /labEandP/ lpdg,EandP
C      common /cmnucleons/ p1cm,p2cm,p3cm,p4cm
C      common /labnucleons/ npdg,p1lab,p2lab,p3lab,p4lab

C      common /abs_constants/ pi,p13,p23,p14,p34,
C     &                       p37,p47,p18,p12,p78

      pi = acos(-1.)
      p13 = 1./3.
      p23 = 2./3.
      p14 = 1./4.
      p34 = 3./4.
      p37 = 3./7.
      p47 = 4./7.
      p18 = 1./8.
      p12 = 1./2.
      p78 = 7./8.

c  Expect to deal with nu N -> l Delta events
c  Start by determining Delta energy and momentum

      if (abs(ipvc(3)).eq.11) then
         ml = me
      else if (abs(ipvc(3)).eq.13) then
         ml = mmu
      else if (abs(ipvc(3)).eq.15) then
         ml = mtau
      else if (abs(ipvc(3)).eq.12 .or. abs(ipvc(3)).eq.14
     $         .or. abs(ipvc(3)).eq.16) then
         ml = 0
      else
         write(*,*) "ABSDELTA Warning: Unknown outgoing lepton = ",
     $                ipvc(3)
         return
      end if

      enu = sqrt(pvc(1,1)**2+pvc(2,1)**2+pvc(3,1)**2)/1000.
      en  = sqrt((pvc(1,2)/1000.)**2+(pvc(2,2)/1000.)**2+
     &     (pvc(3,2)/1000.)**2+mn**2)
      el  = sqrt((pvc(1,3)/1000.)**2+(pvc(2,3)/1000.)**2+
     &     (pvc(3,3)/1000.)**2+ml**2)

      ed  = enu + en - el

      pdx = (pvc(1,1) + pvc(1,2) - pvc(1,3))/1000.
      pdy = (pvc(2,1) + pvc(2,2) - pvc(2,3))/1000.
      pdz = (pvc(3,1) + pvc(3,2) - pvc(3,3))/1000.
      pd  = sqrt(pdx**2 + pdy**2 + pdz**2)

      s = ed**2 - pdx**2 - pdy**2 -pdz**2

c  Find equivalent pion kinetic energy

      epi = (s-mpi**2-mn**2)/(2.*mn)
      tpi = epi - mpi

c  Fill elements of /labpion/ common block
c  assuming pion direction is same as Delta

      if (pd.le.0.) return
      ppi      = sqrt(epi**2 + mpi**2)
      pilab(1) = epi
      pilab(2) = ppi*pdx/pd
      pilab(3) = ppi*pdy/pd
      pilab(4) = ppi*pdz/pd


C  This happens very rarely, and the problem is somewhere upstream
C  in determining the event kinematics, but I don't know exactly (-Patrick)
      if (tpi.le.0) then
         return
      end if

C      write(*,*) "ipvc = ",ipvc(1), ipvc(2), ipvc(3), ipvc(4)
C      write(*,*) "pvc(1) = ",pvc(1,1), pvc(2,1), pvc(3,1)
C      write(*,*) "pvc(2) = ",pvc(1,2), pvc(2,2), pvc(3,2)
C      write(*,*) "pvc(3) = ",pvc(1,3), pvc(2,3), pvc(3,3)
C      write(*,*) "ml = ",ml
C      write(*,*) "mn = ",mn
C      write(*,*) "enu = ",enu
C      write(*,*) "en = ",en
C      write(*,*) "el = ",el
C      write(*,*) "ed = ",ed
C      write(*,*) "pdx = ",pdx
C      write(*,*) "pdy = ",pdy
C      write(*,*) "pdz = ",pdz
C      write(*,*) "s = ",s
C      write (6,*) "TPI = ",tpi

c  Determine 2-, 3-, and 4-body absorption
c  probabilities based on pion energy

      numn = absprob(tpi)

      if (numn.eq.3) go to 1000
      if (numn.eq.4) go to 2000

c  This section for 2-body absorption
c------------------------------------

c  Choose momentum for 2nd nucleon

      call ranlux(rndm,3)

      p1  = (rndm(1)*pfmax**3)**(1./3.)
      th1 = acos(1.-2.*rndm(2))
      ph1 = 2.*pi*rndm(3)
      p1x = p1*cos(ph1)*sin(th1)
      p1y = p1*sin(ph1)*sin(th1)
      p1z = p1*cos(th1)

c  Find total lab energy and momentum
c  Stored in 4-vector plab

      eandp(1) = ed + mn - be
      eandp(2) = pdx + p1x
      eandp(3) = pdy + p1y
      eandp(4) = pdz + p1z

      call abstwobody(tpi)

c  Decide of charge states of the outgoing nucleons

      call ranlux(rndm,1)

c     Delta++
      if (ipvc(4).eq.2224) then
         npdg(1) = 2212
         npdg(2) = 2212

c     Delta+
      else if (ipvc(4).eq.2214) then
         if (rndm(1).lt.p12) then
            npdg(1) = 2212
            npdg(2) = 2212
         else if ((rndm(1).ge.p12).and.(rndm(1).lt.p34)) then
            npdg(1) = 2212
            npdg(2) = 2112
         else
            npdg(1) = 2112
            npdg(2) = 2212
         end if

c     Delta0
      else if (ipvc(4).eq.2114) then
         if (rndm(1).lt.p12) then
            npdg(1) = 2112
            npdg(2) = 2112
         else if ((rndm(1).ge.p12).and.(rndm(1).lt.p34)) then
            npdg(1) = 2212
            npdg(2) = 2112
         else
            npdg(1) = 2112
            npdg(2) = 2212
         end if

c     Delta-
      else if (ipvc(4).eq.1114) then
         npdg(1) = 2112
         npdg(2) = 2112

      end if

c  Fill vcwork
      if (nvc.GT.(maxvc-2)) then
         nvc = maxvc - 2
      end if

      ipvc(nvc+1)  = npdg(1)
      icrnvc(nvc+1) = 1
      ivtivc(nvc+1) = 1
      iorgvc(nvc+1) = 4
      pvc(1,nvc+1) = p1lab(2)*1000.
      pvc(2,nvc+1) = p1lab(3)*1000.
      pvc(3,nvc+1) = p1lab(4)*1000.
      posnuc(1,nvc+1) = posnuc(1,4)
      posnuc(2,nvc+1) = posnuc(2,4)
      posnuc(3,nvc+1) = posnuc(3,4)
      iflgvc(nvc+1) = 0
      CALL MCMASS(IPVC(NVC+1),AMASVC(NVC+1))

      ipvc(nvc+2)  = npdg(2)
      icrnvc(nvc+2) = 1
      ivtivc(nvc+2) = 1
      iorgvc(nvc+2) = 4
      pvc(1,nvc+2) = p2lab(2)*1000.
      pvc(2,nvc+2) = p2lab(3)*1000.
      pvc(3,nvc+2) = p2lab(4)*1000.
      posnuc(1,nvc+2) = posnuc(1,4)
      posnuc(2,nvc+2) = posnuc(2,4)
      posnuc(3,nvc+2) = posnuc(3,4)
      iflgvc(nvc+2) = 0
      CALL MCMASS(IPVC(NVC+2),AMASVC(NVC+2))

      nvc = nvc + 2



      return

c  This section for 3-body absorption
c------------------------------------

 1000 continue

c  Choose momentum for interacting nucleons

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

      eandp(1) = ed + 2.*mn - 2.*be
      eandp(2) = pdx + p1x + p2x
      eandp(3) = pdy + p1y + p2y
      eandp(4) = pdz + p1z + p2z

      call absthreebody

c  Decide of charge states of the outgoing nucleons

      call ranlux(rndm,1)

c     Delta++
      if (ipvc(4).eq.2224) then
         if (rndm(1).lt.p13) then
            npdg(1) = 2212
            npdg(2) = 2212
            npdg(3) = 2112
         else
            npdg(1) = 2212
            npdg(2) = 2212
            npdg(3) = 2212
         end if

c     Delta+
      else if (ipvc(4).eq.2214) then
         if (rndm(1).lt.p14) then
            npdg(1) = 2212
            npdg(2) = 2212
            npdg(3) = 2212
         else if (rndm(1).ge.p34) then
            npdg(1) = 2212
            npdg(2) = 2112
            npdg(3) = 2112
         else
            npdg(1) = 2212
            npdg(2) = 2212
            npdg(3) = 2112
         end if

c     Delta0
      else if (ipvc(4).eq.2114) then
         if (rndm(1).lt.p14) then
            npdg(1) = 2112
            npdg(2) = 2112
            npdg(3) = 2112
         else if (rndm(1).ge.p34) then
            npdg(1) = 2212
            npdg(2) = 2212
            npdg(3) = 2112
         else
            npdg(1) = 2212
            npdg(2) = 2112
            npdg(3) = 2112
         end if

c     Delta-
      else if (ipvc(4).eq.1114) then
         if (rndm(1).lt.p13) then
            npdg(1) = 2212
            npdg(2) = 2112
            npdg(3) = 2112
         else
            npdg(1) = 2112
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
      iorgvc(nvc+1) = 4
      iflgvc(nvc+1) = 0
      pvc(1,nvc+1) = p1lab(2)*1000.
      pvc(2,nvc+1) = p1lab(3)*1000.
      pvc(3,nvc+1) = p1lab(4)*1000.
      posnuc(1,nvc+1) = posnuc(1,4)
      posnuc(2,nvc+1) = posnuc(2,4)
      posnuc(3,nvc+1) = posnuc(3,4)
      CALL MCMASS(IPVC(NVC+1),AMASVC(NVC+1))

      ipvc(nvc+2)  = npdg(2)
      icrnvc(nvc+2) = 1
      ivtivc(nvc+2) = 1
      iorgvc(nvc+2) = 4
      iflgvc(nvc+2) = 0
      pvc(1,nvc+2) = p2lab(2)*1000.
      pvc(2,nvc+2) = p2lab(3)*1000.
      pvc(3,nvc+2) = p2lab(4)*1000.
      posnuc(1,nvc+2) = posnuc(1,4)
      posnuc(2,nvc+2) = posnuc(2,4)
      posnuc(3,nvc+2) = posnuc(3,4)
      CALL MCMASS(IPVC(NVC+2),AMASVC(NVC+2))

      ipvc(nvc+3)  = npdg(3)
      icrnvc(nvc+3) = 1
      ivtivc(nvc+3) = 1
      iorgvc(nvc+3) = 4
      iflgvc(nvc+3) = 0
      pvc(1,nvc+3) = p3lab(2)*1000.
      pvc(2,nvc+3) = p3lab(3)*1000.
      pvc(3,nvc+3) = p3lab(4)*1000.
      posnuc(1,nvc+3) = posnuc(1,4)
      posnuc(2,nvc+3) = posnuc(2,4)
      posnuc(3,nvc+3) = posnuc(3,4)
      CALL MCMASS(IPVC(NVC+3),AMASVC(NVC+3))

      nvc = nvc + 3

      return

c  This section for 4-body absorption
c------------------------------------

 2000 continue

c  Choose momentum for interacting nucleons

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

      eandp(1) = ed + 3.*mn - 3.*be
      eandp(2) = pdx + p1x + p2x + p3x
      eandp(3) = pdy + p1y + p2y + p3y
      eandp(4) = pdz + p1z + p2z + p3z

      call absfourbody

c  Decide of charge states of the outgoing nucleons

      call ranlux(rndm,1)

c     Delta++
      if (ipvc(4).eq.2224) then
         if (rndm(1).lt.p37) then
            npdg(1) = 2212
            npdg(2) = 2212
            npdg(3) = 2212
            npdg(4) = 2212
         else if (rndm(1).gt.p47) then
            npdg(1) = 2212
            npdg(2) = 2212
            npdg(3) = 2212
            npdg(4) = 2112
         else
            npdg(1) = 2212
            npdg(2) = 2212
            npdg(3) = 2112
            npdg(4) = 2112
         end if

c     Delta+
      else if (ipvc(4).eq.2214) then
         if (rndm(1).lt.p18) then
            npdg(1) = 2212
            npdg(2) = 2212
            npdg(3) = 2212
            npdg(4) = 2212
         else if ((rndm(1).ge.p18).and.(rndm(1).lt.p12)) then
            npdg(1) = 2212
            npdg(2) = 2212
            npdg(3) = 2212
            npdg(4) = 2112
         else if ((rndm(1).ge.p12).and.(rndm(1).lt.p78)) then
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

c     Delta0
      else if (ipvc(4).eq.2114) then
         if (rndm(1).lt.p18) then
            npdg(1) = 2212
            npdg(2) = 2212
            npdg(3) = 2212
            npdg(4) = 2112
         else if ((rndm(1).ge.p18).and.(rndm(1).lt.p12)) then
            npdg(1) = 2212
            npdg(2) = 2212
            npdg(3) = 2112
            npdg(4) = 2112
         else if ((rndm(1).ge.p12).and.(rndm(1).lt.p78)) then
            npdg(1) = 2212
            npdg(2) = 2112
            npdg(3) = 2112
            npdg(4) = 2112
         else
            npdg(1) = 2112
            npdg(2) = 2112
            npdg(3) = 2112
            npdg(4) = 2112
         end if

c     Delta-
      else if (ipvc(4).eq.1114) then
         if (rndm(1).lt.p37) then
            npdg(1) = 2112
            npdg(2) = 2112
            npdg(3) = 2112
            npdg(4) = 2112
         else if (rndm(1).gt.p47) then
            npdg(1) = 2212
            npdg(2) = 2112
            npdg(3) = 2112
            npdg(4) = 2112
         else
            npdg(1) = 2112
            npdg(2) = 2112
            npdg(3) = 2112
            npdg(4) = 2112
         end if
      end if

c     Fill vcwork
      if (nvc.GT.(maxvc-4)) then
         nvc = maxvc - 4
      end if

      ipvc(nvc+1)  = npdg(1)
      icrnvc(nvc+1) = 1
      ivtivc(nvc+1) = 1
      iorgvc(nvc+1) = 4
      pvc(1,nvc+1) = p1lab(2)*1000.
      pvc(2,nvc+1) = p1lab(3)*1000.
      pvc(3,nvc+1) = p1lab(4)*1000.
      posnuc(1,nvc+1) = posnuc(1,4)
      posnuc(2,nvc+1) = posnuc(2,4)
      posnuc(3,nvc+1) = posnuc(3,4)
      iflgvc(nvc+1) = 0
      CALL MCMASS(IPVC(NVC+1),AMASVC(NVC+1))

      ipvc(nvc+2)  = npdg(2)
      icrnvc(nvc+2) = 1
      ivtivc(nvc+2) = 1
      iorgvc(nvc+2) = 4
      pvc(1,nvc+2) = p2lab(2)*1000.
      pvc(2,nvc+2) = p2lab(3)*1000.
      pvc(3,nvc+2) = p2lab(4)*1000.
      posnuc(1,nvc+2) = posnuc(1,4)
      posnuc(2,nvc+2) = posnuc(2,4)
      posnuc(3,nvc+2) = posnuc(3,4)
      iflgvc(nvc+2) = 0
      CALL MCMASS(IPVC(NVC+2),AMASVC(NVC+2))

      ipvc(nvc+3)  = npdg(3)
      icrnvc(nvc+3) = 1
      ivtivc(nvc+3) = 1
      iorgvc(nvc+3) = 4
      pvc(1,nvc+3) = p3lab(2)*1000.
      pvc(2,nvc+3) = p3lab(3)*1000.
      pvc(3,nvc+3) = p3lab(4)*1000.
      posnuc(1,nvc+3) = posnuc(1,4)
      posnuc(2,nvc+3) = posnuc(2,4)
      posnuc(3,nvc+3) = posnuc(3,4)
      iflgvc(nvc+3) = 0
      CALL MCMASS(IPVC(NVC+3),AMASVC(NVC+3))

      ipvc(nvc+4)  = npdg(4)
      icrnvc(nvc+4) = 1
      ivtivc(nvc+4) = 1
      iorgvc(nvc+4) = 4
      pvc(1,nvc+4) = p4lab(2)*1000.
      pvc(2,nvc+4) = p4lab(3)*1000.
      pvc(3,nvc+4) = p4lab(4)*1000.
      posnuc(1,nvc+4) = posnuc(1,4)
      posnuc(2,nvc+4) = posnuc(2,4)
      posnuc(3,nvc+4) = posnuc(3,4)
      iflgvc(nvc+4) = 0
      CALL MCMASS(IPVC(NVC+4),AMASVC(NVC+4))

      nvc = nvc + 4

      return

      end
