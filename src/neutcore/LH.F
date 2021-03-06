c subroutine that contracts the leptonic and free hadronic tensors
      
c      include "fourVector.F"

      double precision function LH(k,kPrime,p,qTilde,IPAR,nc,nucPdg, 
     $     NUMATOM)
        implicit none
#include "../neutcore/neutmodel.h"
#include "../neutcore/rscons.h"
c#define __GFORTRAN__
c ----- neutrino pdg, nucleon pdg, nc logical flag        
        integer IPAR
        integer nucPdg, NUMATOM
        logical nc
        double precision k(4), kPrime(4), p(4), qTilde(4)
        double precision leptonMass, inNucleonMass, outNucleonMass
        common /Masses/ leptonMass, inNucleonMass, outNucleonMass
        double precision GEp_di, GEn_di, GMp_di, GMn_di
        logical BBBA5, BBBA7
        real*8 gep, gmp, gen, gmn
        real*8 q2
        common /ff/ gep, gmp, gen, gmn, BBBA5, BBBA7
        double precision h1, h2, h3, h4, h5
        double precision tau, temVal
        double precision f1, f2, fa, fp, fv3, fa3
        double precision CalcFp, GM, GE

        double precision ELAXFF
        external ELAXFF

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

cRT: set masses for NC & CC interactions
        if (nc .eqv. .FALSE.) then
           if (abs(IPAR) .eq. 12) then
              leptonMass = XME
           else if (abs(IPAR) .eq. 14) then
              leptonMass = XMMU
           else if (abs(IPAR) .eq. 16) then
              leptonMass = XMTAU
           else
              leptonMass = 0
           endif
           if (IPAR .gt. 0) then
              inNucleonMass = XMNE
              outNucleonMass = XMP
           else 
              inNucleonMass = XMP
              outNucleonMass = XMNE
           endif
        else 
           leptonMass = 0
           if (nucPDG.eq.2212) then
              inNucleonMass = XMP
              outNucleonMass = XMP
           else
              inNucleonMass = XMNE
              outNucleonMass = XMNE
           endif
        endif

        if (BBBA5 .EQV. .TRUE.) then
#ifdef __GFORTRAN__
c          q2 = real(-dotProduct(qTilde, qTilde)/1000000.0,8)
          q2 = real(-dotProduct(qTilde, qTilde),8)
#else
c          q2 = dble(-dotProduct(qTilde, qTilde)/1000000.0)
          q2 = dble(-dotProduct(qTilde, qTilde))
#endif
          call BBBA05(q2, gep, gmp,
     &    gen, gmn)
c          write (*,*) gep, gmp, gen, gmn
        else if (BBBA7 .EQV. .TRUE.) then
#ifdef __GFORTRAN__
c          q2 = real(-dotProduct(qTilde, qTilde)/1000000.0,8)
          q2 = real(-dotProduct(qTilde, qTilde),8)
#else
c          q2 = dble(-dotProduct(qTilde, qTilde)/1000000.0)
          q2 = dble(-dotProduct(qTilde, qTilde))
#endif          
          call BBBA07(q2, gep, gmp, gen, gmn)
c          write (*,*) gep, gmp, gen, gmn
        else  
          gep = GEp_di(qTilde)
          gmp = GMp_di(qTilde)
          gen = GEn_di(qTilde)
          gmn = GMn_di(qTilde)
c          write (*,*) gep, gmp, gen, gmn
        endif
        
        tau = -1 * dotProduct( qTilde, qTilde )/(4 * outNucleonMass**2)

C     Transverse enhancement -- Modify the magnetic FFs
        if (mod(MDLQE, 1000)/100 .EQ. 7 .AND. NUMATOM .GT. 11) then
           temVal = sqrt(1 + 5.19*q2*EXP(-q2/0.376))
           gmp = gmp*temVal
           gmn = gmn*temVal
        endif

c form factors
        f1 = (GE(qTilde,nc,sgn) + tau * GM(qTilde,nc,sgn))/(1 + tau)
        
        f2 = ( GM(qTilde,nc,sgn) - GE(qTilde,nc,sgn) ) / (1 + tau)

#ifdef __GFORTRAN__
        q2 = real(-dotProduct(qTilde, qTilde)/1000000.0,8)
#else
        q2 = dble(-dotProduct(qTilde, qTilde)/1000000.0)
#endif  
        fa = ELAXFF(q2,nc,sgn)

c RT:Allow for error term on fp as well for this form factor        
        fp = FPQE * CalcFp(qTilde,nc,sgn)

c second-class form factors - off by default.  Use card parameter SCC to turn on
c fv3 is treated as a dipole (using same vector mass)
c fa3 is treated as a scaling of fa
        fv3 = SCCFV*
     +  (-1.) /( (1. - (dotProduct(qTilde, qTilde)/0.71))**2)
        fa3 = SCCFA * fa 
        
c  terms from free nucleon tensor
        h1 = fa**2 *(1.0 + tau) + tau*(f1 + f2)*(f1 + f2)
        h2 = fa**2 + f1**2 + tau*(f2**2 + 4*fa3**2)
        h3 = 2.0*fa*(f1 + f2)
        h4 = 0.25*(f2**2 -tau*(f2 - 2*fv3)**2 - 4*tau*(fp+fa3)**2) -
     +  fv3**2 - 0.5*(f1*(2*fv3 - f2) - 2*fa*(fp + fa3))
        h5 = h2 + 2*(fv3*(f1 - tau*f2) - fa3*(fa - 2*tau*fp))

c calculate 5 individual terms for LH
        term1 = 2.0*dotProduct(k, kPrime)*(inNucleonMass**2)
        
        term2 = 2.0*dotProduct(k, p)*dotProduct(kPrime, p)
     +          - dotProduct(k, kPrime)*dotProduct(p, p)
        
        term3 = dotProduct(k, qTilde)*dotProduct(kPrime, p)
     +           -  dotProduct(k, p)*dotProduct(kPrime, qTilde)
        
        term4 = dotProduct(k, kPrime)*dotProduct(qTilde, qTilde)
     +  - 2.0*dotProduct(k, qTilde)*dotProduct(kPrime, qTilde)
     +          
        
        term5 = dotProduct(k, p)*dotProduct(kPrime, qTilde)
     +          + dotProduct(kPrime, p)*dotProduct(k, qTilde)
     +          - dotProduct(k, kPrime)*dotProduct(p, qTilde)
        
c term3 becomes negative for antineutrinos        
        if (IPAR .LT. 0) then
          term3 = -term3
          h5 = h2 - 2*(fv3*(f1 - tau*f2) - fa3*(fa - 2*tau*fp))
        endif

c LH calculated here
        LH = 2 * ( term1*h1 + term2*h2 + term3*h3 + term4*h4 +
     +  term5 * h5)

cRT: possible debugging line here:
c        write (*,*) fv3,fa3

      if (LH .LT.0) then
        write (*,*) "oddly, LH has been evaluated to be negative..."
        write (*,*) LH
c        LH = 0
      endif
      
      return
      end

c pseudo-scalar form factor
      double precision function CalcFP(q,nc,sgn)
        implicit none
        real dummy
        double precision q(4)
        logical nc
        integer sgn
        double precision protonMass, PiMass2, q2
        double precision dotProduct, ELAXFF

        external dotProduct, ELAXFF

        call MCMASSGV(2212, dummy)
        protonMass = dble(dummy)
        call MCMASSGV(211, dummy)
        piMass2 = dble(dummy)**2
        q2 = dotProduct(q, q)

        CalcFP = ELAXFF(q,nc,sgn)*(protonMass**2) * 2.0/(PiMass2 - q2)

        if (nc .EQV. .TRUE.)then
          CalcFP = 0
        endif

        return
        end

c functions calculate electric and magnetic form factors for F1 and F2
c GE = GEp - GEn
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
c        MV2 = 710000
        MV2 = 0.71
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



