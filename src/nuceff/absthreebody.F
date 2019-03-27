c**********************************************************************
c     ------------------------------
      subroutine absthreebody
c     ------------------------------

      implicit none

C      real EandP(4)
C      real p1cm(4),p2cm(4),p3cm(4),p4cm(4)
C      real p1lab(4),p2lab(4),p3lab(4),p4lab(4)

      real esys,psys,s

      real beta(3)

C      integer lpdg,npdg(4)

#include "absneject.h"

C      common /labEandP/ lpdg,EandP
C      common /cmnucleons/ p1cm,p2cm,p3cm,p4cm
C      common /labnucleons/ npdg,p1lab,p2lab,p3lab,p4lab


c  Find parameters relating lab and cm frames

      esys  = eandp(1)
      psys  = sqrt(eandp(2)**2 + eandp(3)**2 + eandp(4)**2)

      beta(1) = eandp(2)/eandp(1)
      beta(2) = eandp(3)/eandp(1)
      beta(3) = eandp(4)/eandp(1)

      s = esys**2 - psys**2

c  Generate nucleons in cm frame

      call absgen_three(s)

c  Boost to lab frame

      call absboost(beta,p1cm,p1lab)
      call absboost(beta,p2cm,p2lab)
      call absboost(beta,p3cm,p3lab)

      return

      end

