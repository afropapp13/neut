c subroutine that contracts the leptonic and free hadronic tensors

c      include "fourVector.F"

      double precision function LH(k_MeV,kPrime_MeV,p_MeV,
     $                             qTilde_MeV,IPAR,nc,nucPdg, NUMATOM)
        implicit none
#include "neutmodel.h"
#include "rscons.h"
c ----- neutrino pdg, nucleon pdg, nc logical flag
        integer IPAR
        integer nucPdg, NUMATOM
        logical nc
        double precision k_MeV(4), kPrime_MeV(4), p_MeV(4), qTilde_MeV(4)
        double precision leptonMass_MeV, inNucleonMass_MeV, outNucleonMass_MeV
        common /Masses/ leptonMass_MeV, inNucleonMass_MeV, outNucleonMass_MeV
        double precision GEp_di, GEn_di, GMp_di, GMn_di
        logical BBBA5, BBBA7
        real*8 gep, gmp, gen, gmn
        double precision Q2_GeV
        common /ff/ gep, gmp, gen, gmn, BBBA5, BBBA7
        double precision h1, h2, h3, h4, h5
        double precision tau, temVal
        double precision f1, f2, fa, fp, fv3, fa3
        double precision CalcFa, CalcFp, GM, GE

        double precision term1, term2, term3, term4, term5

        double precision dotProduct
        external dotProduct
c ----- sgn is negative for an incoming neutron - applies to NC events only
        integer sgn

        if (nucPdg .EQ. 2112) then
          sgn = -1
        else if (nucPdg .EQ. 2212) then
          sgn = +1
        else
          write (*,*) "incoming nucleon not a proton or a neutron!"
          write (*,*) nucPdg
        endif


c        WRITE(*,*) "In LH: Q2_GeV", Q2_GeV, ", k_MeV[4]: ", k_MeV(4), ", p_MeV[4]: ",p_MeV(4), "qTilde_MeV[4]: ", qTilde_MeV(4)

cRT: set masses for NC & CC interactions
C cwret: but the Rein-Sehgal values are in GeV
C        spectral function calculations assumes MeV
        if (nc .eqv. .FALSE.) then
           if (abs(IPAR) .eq. 12) then
              leptonMass_MeV = XME*1000.
           else if (abs(IPAR) .eq. 14) then
              leptonMass_MeV = XMMU*1000.
           else if (abs(IPAR) .eq. 16) then
              leptonMass_MeV = XMTAU*1000.
           else
              leptonMass_MeV = 0
           endif
           if (IPAR .gt. 0) then
              inNucleonMass_MeV = XMNE*1000.
              outNucleonMass_MeV = XMP*1000.
           else
              inNucleonMass_MeV = XMP*1000.
              outNucleonMass_MeV = XMNE*1000.
           endif
        else
           leptonMass_MeV = 0
           if (nucPDG.eq.2212) then
              inNucleonMass_MeV = XMP*1000.
              outNucleonMass_MeV = XMP*1000.
           else
              inNucleonMass_MeV = XMNE*1000.
              outNucleonMass_MeV = XMNE*1000.
           endif
        endif

c ----- decide here whether to use dipole, or BBBA fits for form factors
        BBBA5 = .FALSE.
        BBBA7 = .FALSE.
        if (mod(mod(MDLQE, 1000),100) .EQ. 2) then
          BBBA5 = .TRUE.
        else if (mod(mod(MDLQE, 1000),100) .EQ. 3) then
          BBBA7 = .TRUE.
        else
          BBBA5 = .FALSE.
          BBBA7 = .FALSE.
        endif

        Q2_GeV = real(-dotProduct(qTilde_MeV, qTilde_MeV)/1E6,8)

        if (BBBA5 .EQV. .TRUE.) then
          call BBBA05(Q2_GeV, gep, gmp, gen, gmn)
        else if (BBBA7 .EQV. .TRUE.) then
          call BBBA07(Q2_GeV, gep, gmp, gen, gmn)
        else
          gep = GEp_di(qTilde_MeV)
          gmp = GMp_di(qTilde_MeV)
          gen = GEn_di(qTilde_MeV)
          gmn = GMn_di(qTilde_MeV)
        endif

        tau = -dotProduct( qTilde_MeV, qTilde_MeV )
     $        / (4 * outNucleonMass_MeV**2)

C     Transverse enhancement -- Modify the magnetic FFs
        if (mod(MDLQE, 1000)/100 .EQ. 7 .AND. NUMATOM .GT. 11) then
           temVal = sqrt(1 + 5.19*Q2_GeV*EXP(-Q2_GeV/0.376))
           gmp = gmp*temVal
           gmn = gmn*temVal
        endif

c form factors
        f1 = (GE(qTilde_MeV,nc,sgn) + tau * GM(qTilde_MeV,nc,sgn))/(1 + tau)

        f2 = ( GM(qTilde_MeV,nc,sgn) - GE(qTilde_MeV,nc,sgn) ) / (1 + tau)

        fa = CalcFa(Q2_GeV,nc,sgn)

c RT:Allow for error term on fp as well for this form factor
        fp = FPQE * CalcFp(Q2_GeV,nc,sgn)

c second-class form factors - off by default.  Use card parameter SCC to turn on
c fv3 is treated as a dipole (using same vector mass)
c fa3 is treated as a scaling of fa
        fv3 = SCCFV*  (-1.) /( (1. + Q2_GeV/0.71)**2)
        fa3 = SCCFA * fa

c  terms from free nucleon tensor
        h1 = fa**2 *(1.0 + tau) + tau*(f1 + f2)*(f1 + f2)
        h2 = fa**2 + f1**2 + tau*(f2**2 + 4*fa3**2)
        h3 = 2.0*fa*(f1 + f2)
        h4 = 0.25*(f2**2 -tau*(f2 - 2*fv3)**2 - 4*tau*(fp+fa3)**2) -
     $  fv3**2 - 0.5*(f1*(2*fv3 - f2) - 2*fa*(fp + fa3))
        h5 = h2 + 2*(fv3*(f1 - tau*f2) - fa3*(fa - 2*tau*fp))

c calculate 5 individual terms for LH
        term1 = 2.0*dotProduct(k_MeV, kPrime_MeV)*(inNucleonMass_MeV**2)

        term2 = 2.0*dotProduct(k_MeV, p_MeV)*dotProduct(kPrime_MeV, p_MeV)
     +          - dotProduct(k_MeV, kPrime_MeV)*dotProduct(p_MeV, p_MeV)

        term3 = dotProduct(k_MeV, qTilde_MeV)*dotProduct(kPrime_MeV, p_MeV)
     $           -  dotProduct(k_MeV, p_MeV)*dotProduct(kPrime_MeV, qTilde_MeV)

        term4 = dotProduct(k_MeV, kPrime_MeV)*dotProduct(qTilde_MeV, qTilde_MeV)
     $  - 2.0*dotProduct(k_MeV, qTilde_MeV)*dotProduct(kPrime_MeV, qTilde_MeV)
     $

        term5 = dotProduct(k_MeV, p_MeV)*dotProduct(kPrime_MeV, qTilde_MeV)
     $          + dotProduct(kPrime_MeV, p_MeV)*dotProduct(k_MeV, qTilde_MeV)
     $          - dotProduct(k_MeV, kPrime_MeV)*dotProduct(p_MeV, qTilde_MeV)

c term3 becomes negative for antineutrinos
        if (IPAR .LT. 0) then
          term3 = -term3
          h5 = h2 - 2*(fv3*(f1 - tau*f2) - fa3*(fa - 2*tau*fp))
        endif

c LH calculated here
        LH = 2 * ( term1*h1 + term2*h2 + term3*h3 + term4*h4 +
     $  term5 * h5)

      if (LH .LT.0) then
        write (*,*) "oddly, LH has been evaluated to be negative..."
        write (*,*) LH
c        LH = 0
      endif

      return
      end

c function to calculate Axial form factor
      double precision function CalcFA(Q2_GeV,nc,sgn)
        implicit none
        logical nc
        integer sgn
C
        double precision Q2_GeV
        double precision ELAXFF
        external ELAXFF

        CalcFA = ELAXFF(Q2_GeV,nc,sgn)

        return
        end

c pseudo-scalar form factor
      double precision function CalcFP(Q2_GeV,nc,sgn)
        implicit none
        real partMass_MeV
        logical nc
        integer sgn
        double precision protonMass2_GeV, PiMass2_GeV, Q2_GeV
        double precision CalcFA

        external dotProduct

        call MCMASS(2212, partMass_MeV)
        protonMass2_GeV = real(partMass_MeV/1E3,8)**2
        call MCMASS(211, partMass_MeV)
        piMass2_GeV = real(partMass_MeV/1E3,8)**2

        CalcFP = CalcFA(Q2_GeV,nc,sgn)*protonMass2_GeV
     $           * (2.0/(PiMass2_GeV + Q2_GeV))

        if (nc .EQV. .TRUE.)then
          CalcFP = 0
        endif

        return
        end

c functions calculate electric and magnetic form factors for F1 and F2
      double precision function GE(q,nc,sgn)
        implicit none
        double precision q(4)
        logical nc
        integer sgn
        logical BBBA5, BBBA7
        real*8 gep, gmp, gen, gmn
        common /ff/ gep, gmp, gen, gmn, BBBA5, BBBA7
        double precision sin2ThetaW

c ----- for CC interactions:
        GE = dble(gep - gen)
C ----- or for NC interactions:
        if (nc.EQV..TRUE.) then
          sin2ThetaW = 0.2312
          if (sgn .LT. 0) then
            GE = sgn*0.5*GE - 2*sin2ThetaW*gen
          else
            GE = sgn*0.5*GE - 2*sin2ThetaW*gep
          endif
        endif
        return
        end

c GM = GMp - GMn
      double precision function GM(q,nc,sgn)
        implicit none
        double precision q(4)
        logical nc
        integer sgn
        logical BBBA5, BBBA7
        real*8 gep, gmp, gen, gmn
        common /ff/ gep, gmp, gen, gmn, BBBA5, BBBA7
        double precision sin2ThetaW

c ----- for CC interactions:
        GM = dble(gmp - gmn)

c ----- if nc, use nc form factors - nucleon dependent
        if (nc.EQV..TRUE.) then
          sin2ThetaW = 0.2312
          if (sgn .LT. 0) then
            GM = -0.5*GM - 2*sin2ThetaW*gmn
          else
            GM = 0.5*GM - 2*sin2ThetaW*gmp
          endif
        endif
        return
        end

c GEp_di - proton electric FF - dipole
      double precision function GEp_di(q)
        implicit none
        double precision q(4)
        double precision a, MV2
        double precision dotProduct

        external dotProduct
        MV2 = 710000
        a = 1.0 - (dotProduct(q, q) / MV2)
        GEp_di = 1. / (a**2)
      return
      end

c GEn_di - neutron electric FF - dipole
      double precision function GEn_di(q)
        implicit none
        double precision q(4)
        GEn_di = 0
        return
        end

c GMp_di - proton magnetic FF - dipole
      double precision function GMp_di(q)
        implicit none
        double precision q(4)
        double precision GEp_di
        GMp_di = 2.793 * GEp_di(q)
        return
        end

c GMn_di - neutron magnetic FF - dipole
      double precision function GMn_di(q)
        implicit none
        double precision dotProduct
        double precision q(4)
        double precision mu_n
        double precision GEp_di

        external dotProduct

        mu_n = -1.9130
        GMn_di = mu_n * GEp_di(q)
        return
        end

      double precision function  LH_GeV(k_GeV,kPrime_GeV,p_GeV,
     $                                  qTilde_GeV,IPAR,nc,nucPdg, NUMATOM)
        implicit none
        integer IPAR
        integer nucPdg, NUMATOM
        logical nc
        integer i
        double precision k_GeV(4), kPrime_GeV(4), p_GeV(4), qTilde_GeV(4)
        double precision k_MeV(4), kPrime_MeV(4), p_MeV(4), qTilde_MeV(4)
        double precision LH

        do 2207 i = 1, 4
          k_MeV(i) = k_GeV(i) * 1.0E3
          kPrime_MeV(i) = kPrime_GeV(i) * 1.0E3
          p_MeV(i) = p_GeV(i) * 1.0E3
          qTilde_MeV(i) = qTilde_GeV(i) * 1.0E3
2207    continue

        LH_GeV = LH(k_MeV,kPrime_MeV,p_MeV,qTilde_MeV,IPAR,nc,nucPdg, NUMATOM)
        return
        end
