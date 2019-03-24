* Author - A. Furmanski 2013
*
* generates elastic events for nuclei with spectral functions
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
*   ETilde - sf parameter related to binding energy
*   pPrime - outgoing nucleon 4-vector
*   q - 4-momentum transfer
*   qTilde - reduced 4-momentum transfer
*
*   see fourVector.F for 4-vector implementations and subroutines
*
******************************************
      subroutine sfevent(IPAR, E, DIRNEU, IERR, single, mode, first)
        implicit none

#include "nework.h"
c ----- temporary values to copy to output common block
        integer tmp_MODENE
        integer tmp_NUMNE
        integer tmp_IPNE(5)
        integer tmp_IORGNE(5)
        integer tmp_IFLGNE(5)
        integer tmp_ICRNNE(5)

        
#include "../neutcore/necard.h"
#include "../neutcore/neutmodel.h"

        real E, DIRNEU(3)
        integer IPAR, IERR, mode
        logical single, nc, first

c ----- pdg code of target
        integer targetType

C ----- dummy indices when needed
        integer i, j
        real dummy

C ----- masses
        double precision leptonMass
        double precision inNucleonMass
        double precision outNucleonMass
        common /masses/ leptonMass, inNucleonMass, outNucleonMass

C        double precision dQQ
        REAL dQQ

C ----- 4-vectors
        double precision k(4), kPrime(4), p(4), pPrime(4)
        double precision q(4), qTilde(4)

C ----- centre of mass (com) related variables
C ----- temporary vector for use when needed
        double precision kPrimeCom(4), comVel(4)
        double precision temp_vector(4)

C ----- common block for dealing with spectral function selection
# include "../specfunc/specfunc.h"

C ----- terms in the cross-section
        double precision initialTerm, volumeTerm, jacobian, crossSection
        double precision tableMax
C ----- when calculating tables or total cross-sections
        double precision mean, maxValue
        double precision maxValueNC(2), meanNC(2)
        integer N_tries
        common /makeTable/ mean, maxValue, N_tries, meanNC, maxValueNC

C ----- external functions
        double precision modMomentum, LH, dotProduct, maxDiff
        logical isCorr

        real RLU, FNELS_RPASCL
        external modMomentum, LH, dotProduct, RLU, maxDiff, 
     &  FNELS_RPASCL
C ----- parameters        
        double precision Gf, cosThetaC, Pi

        parameter(Gf = 1.16639E-11,cosThetaC = 0.97418
     &  ,Pi = 3.14159265358)
C ******************************************************** C
c        write (*,*) "MAKING CCQE EVENT"
c        write (*,*) numbndp,numatom, mode, E
c ----- E is in GeV here - change it to MeV for the calculation

        if ((NUMBNDP .ne. 6) .AND. (NUMBNDP .ne. 8) .AND. 
     &      (NUMBNDP .ne. 26)) then
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
c        MODENE = mode

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
        call buildSF
C ----- find the maximum allowed value to throw against - note energy in GeV in tables, MeV in code.
c ----- Carbon-12

        if (NUMBNDP .EQ. 6)then
          targetType = 1000060120
          if (first .EQV. .TRUE.) then
c            call buildSF_C12
            if (single.EQV..TRUE.)then
              tableMax = maxDiff(dble(E)/1000.0, IPAR, targetType, nc)
            endif
            first = .FALSE.
          endif
c ----- Oxygen-16
        else if (NUMBNDP .EQ. 8) then
          targetType = 1000080160
          if (first .EQV. .TRUE.) then
c            call buildSF_O16
            if (single .EQV..TRUE.)then
              tableMax = maxDiff(dble(E)/1000.0, IPAR, targetType, nc)
            endif
            first = .FALSE.
          endif
c ----- Iron-56
        else if (NUMBNDP .EQ. 26) then
          targetType = 1000260560
          if (first .EQV. .TRUE.) then
c            call buildSF_Fe56
            if (single .EQV..TRUE.)then
              tableMax = maxDiff(dble(E)/1000.0, IPAR, targetType, nc)
            endif
            first = .FALSE.
          endif
c ----- Deuterium - currently not used - smith-moniz used instead
        else if ((NUMBNDP .EQ. 1) .AND. (NUMBNDN .EQ. 1)) then
          targetType = 1000010020
          if (first .EQV. .TRUE.) then
            if (single .EQV..TRUE.)then
              tableMax = maxDiff(dble(E)/1000.0, IPAR, targetType, nc)
            endif
            first = .FALSE.
          endif
        endif
c        write (*,*) "tablemax = ", tableMax

        
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

C ----- select values from the spectral function array
        if (NUMATOM .LT. 3) then
          p(2) = 0
          p(3) = 0
          p(4) = 0
          ETilde = 0
          fermiMomentum = 0
        else
          call selectSfValues(p)
        endif


        p(1) = sqrt(inNucleonMass**2 + modMomentum(p)**2)

C ----- temp_vector is the c.o.m 4-momentum
        do i = 2, 4
          temp_vector(i) = k(i) + p(i)
        enddo
c ----- energy conservation from spectral function definition
        temp_vector(1) = k(1) + inNucleonMass - ETilde

        if (dotProduct(temp_vector, temp_vector) .LT.
     +  (leptonMass + outNucleonMass)**2 ) then
c ----- not enough energy - skip straight to end of event
          crossSection = 0
          goto 20
        endif
C ----- calculate outgoing momenta, and pick a direction
        kPrime(1) = (dotProduct(temp_vector, temp_vector) 
     &  - outNucleonMass**2 + leptonMass**2)
     &  / (2 * sqrt(dotProduct(temp_vector, temp_vector)) )
        call pickDirection(kPrime, leptonMass)

        do i = 2, 4
          pPrime(i) = - kPrime(i)
        enddo
        pPrime(1) = sqrt(modMomentum(pPrime)**2 + outNucleonMass**2)

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
          qTilde(i) = q(i)
        enddo
        qTilde(1) = pPrime(1) - p(1)

C ----- pauli blocking
        if (modMomentum(pPrime) .LT. fermiMomentum) then
c event impossible due to pauli blocking - go straight to end of event          
          crossSection = 0
          goto 20
c          write(*,*) "pauliBlocking"
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


c ----- RPA correction to dsigma/dQ2 spectrum.
c ----- Not fully clear whether this is double counting effects.
c ----- Not used for total cross-section calculation, only vector generation
        if ((mod(MDLQE,10000)/1000 .eq. 1)
     &     .and.(single .eqv. .TRUE.)) then
           dQQ = REAL(dotProduct(q,q)/(-1.e6))
           crossSection=crossSection*FNELS_RPASCL(dQQ,IPAR)
        endif

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
     &                  (dble(N_tries)*dble(NUMBNDP))
            meanNC(2) = meanNC(2) *dble(NUMATOM) /
     &                  (dble(N_tries)*dble(NUMBNDN))
            E = E / 1000.0 ! turn back into GeV for the rest of neut
            return
          endif
c event complete - go back to start to generate more
          goto 10
        


C ----- for event generation -  compare to value loaded from xsec table
c ----- 1% factor reduces efficiency but avoids errors from interpolation
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
      
      
c ----- Add second nucleon, if interaction in correlated tail
        if (isCorr(p, Etilde, targetType) .eqv..true.) then
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

***********************************************************************

C ----- Jacobian converts from c.o.m integration to lab integration
        double precision function Jacobian(kPrime, kPrimeCom, pPrime
     &   , comVel)
          implicit none
          double precision kPrime(4), kPrimeCom(4), pPrime(4), comVel(4)
        double precision vkPrime(4), vpPrime(4), velDiff(4)
        double precision gam, cosTheta0, sinTheta0
        double precision modMomentum
        integer i

        external modMomentum
        
        gam = 1. / sqrt(1. - modMomentum(comVel)**2)

        cosTheta0 = (kPrimeCom(2)*comVel(2) + kPrimeCom(3)*comVel(3)
     +               + kPrimeCom(4)*comVel(4)) 
     + /( modMomentum(kPrimeCom) * modMomentum(comVel) )
        
        sinTheta0 = sqrt( 1. - cosTheta0**2 )

        do i = 2, 4
          vkPrime(i) = kPrime(i) / kPrime(1)
          vpPrime(i) = pPrime(i) / pPrime(1)
        enddo
        do i = 2, 4
          velDiff(i) = vkPrime(i) - vpPrime(i)
        enddo
        
        
c        Jacobian = sqrt(sinTheta0**2 + gam*(cosTheta0**2) )
        Jacobian = sqrt(1 + (1 - cosTheta0**2)*(gam**2 - 1))
     +             / modMomentum(velDiff)


        return
        end


***********************************************************************
* function returns true if initial nucleon comes from correlated tail
***********************************************************************
        logical function isCorr(p, Etilde, targetType)
          double precision p(4), Etilde
          integer targetType

          double precision pMax, eMax

          double precision modMomentum
          external modMomentum

          pMax = 0
          eMax = 0

          if (targetType .eq. 1000060120) then ! C12
            pMax = 300 ! return to values later
            eMax = 100 ! return to values later
          endif

          if (targetType .eq. 1000080160) then ! O16
            pMax = 300 ! return to values later
            eMax = 100 ! return to values later
          endif

          if (targetType .eq. 1000260560) then ! Fe56
            pMax = 300 ! return to values later
            eMax = 100 ! return to values later
          endif
          
C         write(*,*) modMomentum(p), Etilde
          if ((modMomentum(p) .gt. pMax) .or. (Etilde .gt. eMax))then
            isCorr = .TRUE.
            return 
          endif

          isCorr = .FALSE.
          
          return
          end
