      integer function absprob(tpigev)
c
c     ( purpose )
c       Choose nucleon multiplicity following pion absorption
c
c     ( input )
c       tpigev  incident pion kinetic energy in lab frame (in GeV)
c
c     ( output )
c       absprob = 2, 3, or 4
c
c     algorithm:
c        2-body absorption is 80% at tpi=0, 0% at tpi=400 MeV
c        4-body is 0% at tpi=50 MeV, 80% at tpi=450 MeV
c        otherwise is is 3-body
c

      implicit none

      real tpigev,tpi
      real a2,b2,a4,b4
      real p2,p3,p4,rndm, p24

      a2 = 0.8
      b2 = -0.8/400.
      b4 = 0.8/450.
      a4 = -50.*b4

      tpi = 1000.*tpigev

      call ranlux(rndm,1)

      if (tpi.lt.50.) then
         p2 = a2 + b2*tpi
         if (rndm.gt.p2) then
            absprob = 3
         else
            absprob = 2
         end if
         return
      else if (tpi.ge.612.5) then
         absprob = 4
         return
      else if ((tpi.gt.400.).and.(tpi.lt.612.5)) then
         p4 = a4 + b4*tpi
         if (rndm.gt.p4) then
            absprob = 3
         else
            absprob = 4
         end if
         return
      else
         p2 = a2 + b2*tpi
         p4 = a4 + b4*tpi
         p3 = 1. - p2 - p4
         p24 = p2 + p4
         if (rndm.gt.p24) then
            absprob = 3
         else if (rndm.gt.p2) then
            absprob = 4
         else
            absprob = 2
         end if
         return
      end if

      end
