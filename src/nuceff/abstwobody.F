c**********************************************************************
c     ------------------------------
      subroutine abstwobody(tpi)
c     ------------------------------

      implicit none

      real tpi

      real mn/0.93893/
      real pi

C      real EandP(4)
C      real p1cm(4),p2cm(4),p3cm(4),p4cm(4)
C      real p1lab(4),p2lab(4),p3lab(4),p4lab(4)

      real esys,psys,s,w,k
      real ct,st,ph

      real beta(3),mbeta(3)

      real rndm

      real pimom,thpi,phpi,tmp(4),picm(4)

C      integer lpdg,npdg(4)
C
#include "absneject.h"
C      common /labEandP/ lpdg,EandP
C      common /cmnucleons/ p1cm,p2cm,p3cm,p4cm
C      common /labnucleons/ npdg,p1lab,p2lab,p3lab,p4lab

      pi = acos(-1.)
      call ranlux(rndm,1)

c  Find parameters relating lab and cm frames

      esys  = eandp(1)
      psys  = sqrt(eandp(2)**2 + eandp(3)**2 + eandp(4)**2)

      beta(1) = eandp(2)/eandp(1)
      beta(2) = eandp(3)/eandp(1)
      beta(3) = eandp(4)/eandp(1)

      mbeta(1) = -beta(1)
      mbeta(2) = -beta(2)
      mbeta(3) = -beta(3)

      s = esys**2 - psys**2

      w = sqrt(s)/2.
      k = sqrt(w**2 - mn**2)

c  Find 4-momentum of absorbed pion in cm frame

      call absboost(mbeta,pilab,picm)

c  Find rotation angles which put pion cm
c  3-momentum along z-axis

      pimom  = sqrt(picm(2)**2+picm(3)**2+picm(4)**2)
      if (pimom.lt.0.000001) then
        thpi = 0.
        phpi = 0.
      else
        thpi   = acos(picm(4)/pimom)
        phpi   = atan(picm(3)/picm(2))
      end if

c  Find angle of one nucleon in the cm frame
c  in which absorbed pion momentum is along z-axis
c  --> absgen_two returns cosine (ct) of forward
c      going nucleon, i.e. 0 < ct <1

      call absgen_two(tpi,ct)
      st = sqrt(1.-ct**2)
      ph = 2.*pi*rndm

c  Choose whether nucleon 1 or 2 goes forward

      call ranlux(rndm,1)

      if (rndm.lt.0.5) then
         p1cm(1) = w
         p1cm(2) = k*cos(ph)*st
         p1cm(3) = k*sin(ph)*st
         p1cm(4) = k*ct

         p2cm(1) = w
         p2cm(2) = -p1cm(2)
         p2cm(3) = -p1cm(3)
         p2cm(4) = -p1cm(4)
      else
         p2cm(1) = w
         p2cm(2) = k*cos(ph)*st
         p2cm(3) = k*sin(ph)*st
         p2cm(4) = k*ct

         p1cm(1) = w
         p1cm(2) = -p2cm(2)
         p1cm(3) = -p2cm(3)
         p1cm(4) = -p2cm(4)
      end if

c  Now rotate the two nucleon 3-momenta from the
c  frame in which the pion momentum is along z-axis
c  to the actual cm frame

c  First p1
c  rotate around y by thpi
      tmp(2) =  cos(thpi)*p1cm(2) + sin(thpi)*p1cm(4)
      tmp(3) =  p1cm(3)
      tmp(4) = -sin(thpi)*p1cm(2) + cos(thpi)*p1cm(4)
c  rotate around z by -phipi
      p1cm(2) =  cos(-phpi)*tmp(2) + sin(-phpi)*tmp(3)
      p1cm(3) = -sin(-phpi)*tmp(2) + cos(-phpi)*tmp(3)
      p1cm(4) =  tmp(4)

c  Then p2
c  rotate around y by thpi
      tmp(2) =  cos(thpi)*p2cm(2) + sin(thpi)*p2cm(4)
      tmp(3) =  p2cm(3)
      tmp(4) = -sin(thpi)*p2cm(2) + cos(thpi)*p2cm(4)
c  rotate around z by -phipi
      p2cm(2) =  cos(-phpi)*tmp(2) + sin(-phpi)*tmp(3)
      p2cm(3) = -sin(-phpi)*tmp(2) + cos(-phpi)*tmp(3)
      p2cm(4) =  tmp(4)


c  Now boost to lab frame

      call absboost(beta,p1cm,p1lab)
      call absboost(beta,p2cm,p2lab)

      return

      end
