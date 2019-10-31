* Adapted from the SF code - C. Wilkinson 2014
*
***************************************
*
* inputs:
*   ipar - neutrino type
*   E - neutrino energy
*   DIRNEU - neutrino direction
*   IERR - error code
*   mode - NEUT reaction code (should be 1(CC) or 51/52(NC p/n))
*   single - if set to false, code is used to find the total xsec and max weight
*   first - simple flag to do some things on the first pass only
*
* variables
*   k - neutrino 4-vector
*   kPrime - lepton 4-vector
*   p - incoming nucleon 4-vector
*   pPrime - outgoing nucleon 4-vector
*   q - 4-momentum transfer
*
*   see fourVector.F for 4-vector implementations and subroutines
*
******************************************
      subroutine effsfevent(IPAR, E, DIRNEU, IERR, single, mode, first)
        implicit none

#include "nework.h"
c ----- temporary values to copy to output common block
        integer tmp_MODENE
        integer tmp_NUMNE
        integer tmp_IPNE(5)
        integer tmp_IORGNE(5)
        integer tmp_IFLGNE(5)
        integer tmp_ICRNNE(5)


#include "necard.h"
#include "neutmodel.h"

        real E, DIRNEU(3)
        integer IPAR, IERR, mode
        logical single, nc, first, isCorr, isTEM

c ----- pdg code of target
        integer targ

C ----- dummy indices when needed
        integer i, j
        real dummy
        double precision dummyVect(4), testVect(4)
C ----- masses
        double precision leptonMass
        double precision inNucleonMass
        double precision outNucleonMass
        common /masses/ leptonMass, inNucleonMass, outNucleonMass

C        double precision dQQ
C        REAL dQQ

C ----- 4-vectors
        double precision k(4), kPrime(4), p(4), pPrime(4)
        double precision q(4), qTilde(4)

C     CW: new variables required for the effective SF calculation
        double precision removal, testVal, pauli

C ----- centre of mass (com) related variables
C ----- temporary vector for use when needed
        double precision kPrimeCom(4), comVel(4)
        double precision temp_vector(4)

C ----- common block for dealing with spectral function selection
#include "effspecfunc.h"

C ----- terms in the cross-section
        double precision initialTerm, volumeTerm, crossSection
        double precision tableMax

C ----- when calculating tables or total cross-sections
        double precision mean, maxValue
        double precision maxValueNC(2), meanNC(2)
        integer N_tries
        common /makeEffSFTable/ mean,maxValue,N_tries,meanNC,maxValueNC

C ----- external functions
        double precision modMomentum, LH, dotProduct
        double precision Jacobian, effSFMaxDiff

C ----- Functions defined in this file
        logical isCorrEff
        integer getTarget
        double precision pauliBlocking

        real RLU
        external modMomentum, LH, dotProduct, RLU, effSFMaxDiff
C ----- parameters
        double precision Gf, cosThetaC, Pi

        parameter(Gf = 1.16639E-11,cosThetaC = 0.97418
     &  ,Pi = 3.14159265358)

C        write(*,*) "MAKING CCQE EVENT"
c ----- E is in GeV here - change it to MeV for the calculation

        targ = 0
        targ = getTarget(NUMATOM)

C     Check if TEM is used
        if (mod(MDLQE, 1000)/100 .EQ. 7) then
           isTEM = .TRUE.
        else
           isTEM = .FALSE.
        endif

        if (targ .EQ. 0) then
         write (*,*) "using RFG for nucleus with atomic mass", NUMATOM
         if (abs(mode) .eq.1) then
           call NEELSVCP(IPAR,E,DIRNEU,IERR)
           return
         else if ((abs(mode) .eq. 51) .or. (abs(mode).eq.52)) then
           call NEELSVNP(IPAR,E,DIRNEU,IERR,MODE)
           return
         endif
       endif

        nc = (abs(mode) .NE.1)

        E = E * 1000

C ----- set initial neutrino four-vector components
        k(1) = dble(E)
        do i = 2, 4
          k(i) = dble(E*DIRNEU(i-1))
        enddo

C ----- set masses of particles
        if (nc .EQV. .FALSE.) then
          call MCMASS(abs(IPAR)-1, dummy)
          leptonMass = dble(dummy)
        else
          leptonMass = 0
        endif

c nc
        if (nc .EQV. .TRUE.) then
          if ((single .EQV..TRUE. .and. abs(MODE) .EQ. 52)
     &.OR. ((single .EQV. .FALSE. .AND. RLU(DUMMY)
     & .LT. dble(NUMBNDN/dble(NUMATOM))))) then
            call MCMASS(2112, dummy)
            inNucleonMass = dummy
            outNucleonMass = dummy
            tmp_IPNE(2) = 2112
            tmp_IPNE(4) = 2112
            IPNE(2) = 2112
            IPNE(4) = 2112
            tmp_MODENE = mode
          else if ((single .EQV. .TRUE. .AND. abs(MODE) .EQ. 51)
     & .OR. (single .EQV. .FALSE.)) then
            call MCMASS(2212, dummy)
            inNucleonMass = dummy
            outNucleonMass = dummy
            tmp_IPNE(2) = 2212
            tmp_IPNE(4) = 2212
            IPNE(2) = 2212
            IPNE(4) = 2212
            tmp_MODENE = mode
          endif
c cc
        else if (IPAR .LT. 0) then
          call MCMASS(2212, dummy)
          inNucleonMass = dummy
          call MCMASS(2112, dummy)
          outNucleonMass = dummy
          tmp_IPNE(2) = 2212
          tmp_IPNE(4) = 2112
          tmp_MODENE = -1
        else
          call MCMASS(2212, dummy)
          outNucleonMass = dummy
          call MCMASS(2112, dummy)
          inNucleonMass = dummy
          tmp_IPNE(2) = 2112
          tmp_IPNE(4) = 2212
          tmp_MODENE = 1
        endif

C ----- build the map of spectral function values that can be picked from
        call buildEFFSF(XMAQE)
C ----- find the maximum allowed value to throw against - note energy in GeV in tables, MeV in code.

        if (first .EQV. .TRUE.) then
           if (single.EQV..TRUE.)then
             tableMax=effSFMaxDiff(dble(E)/1000.0,IPAR, targ, nc, MDLQE)
           endif
           first = .FALSE.
        endif

c ----- setup temp values for output common block
        tmp_NUMNE = 4
        tmp_IPNE(1) = IPAR
        if (nc .EQV. .FALSE.)then
          tmp_IPNE(3) = sign(abs(IPAR)-1, IPAR)
        else
          tmp_IPNE(3) = IPAR
        endif
        data tmp_IORGNE / 0, 0, 1, 2, 0 /
        data tmp_IFLGNE / -1, -1, 0, 0, 0 /
        data tmp_ICRNNE / 0, 0, 1, 1, 1 / ! note that an additional nucleon is allowed for

C ******************************************************************** C
C ----- Begin event generation -----
c       write(*,*) "Beginning event generation Effective SF"
        N_tries = 0
10      continue

        if (tmp_IPNE(2) < 0) then
          write(*,*) "in nucleon pdg < 0!"
        endif

c ----- For table generation only, if nc, we select proton or neutron again.
c ----- For generating actual events we have already decided.
        if ( single .EQV. .FALSE.) then
            if ((nc .EQV. .TRUE.) .AND. (RLU(dummy) .LT.
     &                dble(NUMBNDN)/dble(NUMATOM))) then
              call MCMASS(2212, dummy)
              inNucleonMass = dummy
              outNucleonMass = dummy
              tmp_IPNE(2) = 2212
              tmp_IPNE(4) = 2212
              IPNE(2) = 2212
              IPNE(4) = 2212
              tmp_MODENE = 51
            else if(nc .EQV. .TRUE.) then
              call MCMASS(2212, dummy)
              inNucleonMass = dummy
              outNucleonMass = dummy
              tmp_IPNE(2) = 2112
              tmp_IPNE(4) = 2112
              IPNE(2) = 2112
              IPNE(4) = 2112
              tmp_MODENE = 52
            endif
        endif

        if (NUMATOM .LT. 3) then
          p(2) = 0
          p(3) = 0
          p(4) = 0
          ETilde = 0
          fermiMomentum = 0
        else
          call selectEffSfValues(p)
        endif
        isCorr = isCorrEff(NUMATOM, isTEM)

        if (isCorr .EQV. .TRUE.) then
C     Deuterium mass
           call MCMASS(1000010020, dummy)
           testVect(1) = dummy - 2*ETilde -
     +          sqrt(outNucleonMass**2 + modMomentum(p)**2)
        else
C     Target nucleus mass
           call MCMASS(targ, dummy)
           testVect(1) = inNucleonMass - ETilde - (modMomentum(p)**2/
     +          (2*(dummy - inNucleonMass + ETilde)))
        endif

C     CW: the cross section has to be calculated with the nucleon on mass shell
        p(1) = sqrt(inNucleonMass**2 + modMomentum(p)**2)
        removal = inNucleonMass - testVect(1)

C ----- temp_vector is the c.o.m 4-momentum
        do i = 2, 4
           testVect(i) = p(i)
           temp_vector(i) = k(i) + p(i)
        enddo
c ----- energy conservation from spectral function definition
C     CW: this is the neutrino-nucleon centre of mass energy
        temp_vector(1) = k(1) + inNucleonMass - removal

        if (dotProduct(temp_vector, temp_vector) .LT.
     +       (leptonMass + outNucleonMass)**2 ) then
c ----- not enough energy - skip straight to end of event
          crossSection = 0
          goto 20
        endif

C ----- calculate outgoing momenta, and pick a direction
C     CW: modify. Can I pick a direction like this. Yes.
C     CW: this needs to be modified.
        kPrime(1) = (dotProduct(temp_vector, temp_vector)
     &  - outNucleonMass**2 + leptonMass**2)
     &  / (2 * sqrt(dotProduct(temp_vector, temp_vector)) )
        call pickDirection(kPrime, leptonMass)

C     CW: momentum conservation!
        do i = 2, 4
          pPrime(i) = - kPrime(i)
        enddo
        pPrime(1) = sqrt(modMomentum(pPrime)**2 + outNucleonMass**2)

C     CW: need to be careful to check that the method used when boosting into COM frame is still valid
C ----- Boost
c temp_vector contains c.o.m momentum and energy
        do i = 2, 4
          comVel(i) = -temp_vector(i) / temp_vector(1)
          kPrimeCom(i) = kPrime(i)
        enddo

        call boost(kPrime, comVel)
        call boost(pPrime, comVel)

C ----- calculate q and qTilde
        do i = 1, 4
           q(i) = k(i) - kPrime(i)
           dummyVect(i) = p(i) + q(i)
           qTilde(i) = q(i)
        enddo

C     CW: Now apply Pauli suppression if q3 < 2*fermiMomentum
        pauli = pauliBlocking(q, fermiMomentum, NUMATOM, nc)
        if (RLU(dummy) .GT. pauli ) then
           crossSection = 0
           goto 20
        endif

***********************************************************************
C ----- Calculate xsec
        initialTerm = Gf**2 * cosThetaC**2
     &   / ( 8*Pi*Pi*k(1)*kPrime(1)*p(1)*pPrime(1) )

c ----- for nc events there is no cabibbo mixing to worry about
        if (nc.EQV..TRUE.)then
          initialTerm = initialTerm / (cosThetaC**2)
        endif

        volumeTerm = 1
     &               * 4*Pi*(modMomentum(kPrimeCom)**2)
     &               * Jacobian(kPrime, kPrimeCom, pPrime, comVel)

        crossSection = initialTerm*volumeTerm *
     &  LH(k, kPrime, p, qTilde, IPAR, nc,tmp_IPNE(2), NUMATOM)

c ----- turn cross-sections into per-nucleon cross-sections in cm^2
        crossSection = crossSection * 3.89E-22

c sent straight here if event is impossible
20      continue
        if (nc .EQV. .FALSE.)then
          mean = mean + crossSection
          maxValue = max(maxValue, crossSection)
        else if (abs(tmp_MODENE) .EQ. 51) then
          meanNC(1) = meanNC(1) + crossSection
          maxValueNC(1) = max(maxValueNC(1), crossSection)
        else
          meanNC(2) = meanNC(2) + crossSection
          maxValueNC(2) = max(maxValueNC(2), crossSection)
        endif
        N_tries = N_tries + 1

**********  For table generation only!
**********  Once enough events generated, stop
**********  This value can be changed to increase precision/speed
        if (single .EQV. .FALSE.) then
           if (N_tries .GT. 1000000) then
              mean = mean / N_tries
              meanNC(1) = meanNC(1) *dble(NUMATOM) /
     &             (dble(N_tries)*dble(NUMBNDP))
              meanNC(2) = meanNC(2) *dble(NUMATOM) /
     &             (dble(N_tries)*dble(NUMBNDN))
              E = E / 1000.0    ! turn back into GeV for the rest of neut
              write(*,*) E, meanNC(1), meanNC(2)
              return
           endif
c     event complete - go back to start to generate more
           goto 10
C     ----- for event generation -  compare to value loaded from xsec table
c     ----- 1% factor reduces efficiency but avoids errors from interpolation
        else if (crossSection*1E38 .LT. RLU(dummy)*tableMax*1.05) then
           goto 10
        endif

c ----- if event has got to here - event is saved
c ----- need to change everything into GeV
        MODENE = tmp_MODENE
        NUMNE = tmp_NUMNE
        do i = 1, 4
          IPNE(i) = tmp_IPNE(i)
          IORGNE(i) = tmp_IORGNE(i)
          IFLGNE(i) = tmp_IFLGNE(i)
          ICRNNE(i) = tmp_ICRNNE(i)
        enddo
        do i = 1, 3
          PNE(i, 1) = real(k(i+1)) / 1000
          PNE(i, 2) = real(p(i+1)) / 1000
          PNE(i, 3) = real(kPrime(i+1)) / 1000
          PNE(i, 4) = real(pPrime(i+1)) / 1000
        enddo

        E = E / 1000 ! turn back into GeV for the rest of neut

C     Simulate second nucleon for a correlated event
        if (isCorr .eqv..true.) then
          NUMNE = 5
          if (IPNE(2) .eq. 2212) then ! all correlated pairs treated as n-p
            IPNE(5) = 2112
          else
            IPNE(5) = 2212
          endif
          IORGNE(5) = tmp_IORGNE(5)
          IFLGNE(5) = tmp_IFLGNE(5)
          ICRNNE(5) = tmp_ICRNNE(5)
          do i = 1,3
            PNE(i,5) = -1*real(p(i+1)/1000)
          enddo
        endif

      return
      end

***********************************************************************
* function returns true if initial nucleon comes from correlated tail
***********************************************************************
        logical function isCorrEff(NUMATOM, isTEM)
          double precision val, f1p1h
          logical isTEM
          real dummy, RLU
          external RLU

          isCorrEff = .FALSE.
          f1p1h     = 1.
          val = RLU(dummy)

          if (NUMATOM .eq. 2) then ! H2
             f1p1h = 0
          else if (NUMATOM .eq. 3) then ! He3
             f1p1h = 0.312
          else if (NUMATOM .eq. 4) then ! He4
             f1p1h = 0.791
          else if (NUMATOM .eq. 12 .OR. NUMATOM .eq. 16) then ! C12 or O16
             f1p1h = 0.808
          else if (NUMATOM .eq. 20) then ! Ne20
             f1p1h = 0.765
          else if (NUMATOM .eq. 27) then ! Al27
             f1p1h = 0.774
          else if (NUMATOM .eq. 40) then ! Ar40
             f1p1h = 0.809
          else if (NUMATOM .eq. 56) then ! Fe56
             f1p1h = 0.822
          else if (NUMATOM .eq. 63 .OR. NUMATOM .eq. 64 .OR.
     &            NUMATOM .eq. 208) then ! Pb208
             f1p1h = 0.896
          endif

C     Modify if TEM is turned on and nucleus is heavier than carbon
          if (isTEM .EQV. .TRUE. .AND. NUMATOM .GT. 11) then
             f1p1h = f1p1h/1.18
          endif

          if (val > f1p1h) then
             isCorrEff = .TRUE.
          endif

          return
          end

      integer function getTarget(NUMATOM)

          if (NUMATOM .eq. 2) then  ! H2
             getTarget = 1000010020
          else if (NUMATOM .eq. 3) then ! He3
             getTarget = 1000020030
          else if (NUMATOM .eq. 4) then ! He4
             getTarget = 1000020040
          else if (NUMATOM .eq. 12) then ! C12
             getTarget = 1000060120
          else if ( NUMATOM .eq. 16) then ! O16
             getTarget = 1000080160
          else if (NUMATOM .eq. 20) then ! Ne20
             getTarget = 1000100200
          else if (NUMATOM .eq. 27) then ! Al27
             getTarget = 1000130270
          else if (NUMATOM .eq. 40) then ! Ar40
             getTarget = 1000180400
          else if (NUMATOM .eq. 56) then ! Fe56
             getTarget = 1000260560
          else if (NUMATOM .eq. 63) then ! Cu63
             getTarget = 1000290630
          else if (NUMATOM .eq. 64) then ! Zn64
             getTarget = 1000300640
          else if (NUMATOM .eq. 208) then ! Pb208
             getTarget = 1000822080
          endif
          return
          end

      double precision function pauliBlocking(q, fermiMomentum, NUMATOM,
     &     nc)
          double precision dQQ, dotProduct, modMomentum, fermiMomentum
          logical nc
          double precision q(4)

          external dotProduct, modMomentum

          pauliBlocking = 1.
          dQQ = REAL(dotProduct(q,q)/(-1.e6))

C     Deuterium is a special case
          if (NUMATOM .EQ. 2) then
C     Pauli blocking only applies for CC events
             if (nc .EQV. .FALSE.) then
                pauliBlocking = 1 - 0.588918*EXP(-17.2306*dQQ)**0.749157
             endif
C     For all other targets except deuterium
          else
             if (modMomentum(q) .LT. 2*fermiMomentum) then
                pauliBlocking = (3./4.)*(modMomentum(q)/fermiMomentum)*
     +               (1 - ((modMomentum(q)/fermiMomentum)**2)/12)
             endif
          endif
          return
          end
