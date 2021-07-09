***********************************************************************
*     ---------------------------------------------------
      SUBROUTINE RSCALCY(IMODE, INLEP, E, Q2, W, theta, phi, yfn)
*     ---------------------------------------------------
*
*     ( purpose )
*       calculate kinematics of final state particles
*
*     ( input )
*       imode       : interaction mode
*       inlep       : kind of incoming neutrino
*       e           : energy of incoming neutrino
*
*     ( output )
*       yfn: angular distribution function (called W(theta,phi) in
*       Rein-Sehgal). Previously called YFN in NEUT
*
*
*     ( creation date and author )
*       2021.Feb. ; C.Wret: Add to reweight pion ejection method
*                           Carved out from rsdt33
*
*     ( comment )
*       Reference: Rein-Sehgal 1981, Appendix
*               
*
*       ----------- IFLAG means Interaction mode ----------
*
*       IFLAG < 10 -> neutrino
*       IFLAG > 10 -> anti neutrino
*
*       IFLAG = 1:'CHARGED CURRENT, P - PI+'
*       IFLAG = 2:'CHARGED CURRENT, P - PI0'
*       IFLAG = 3:'CHARGED CURRENT, N - PI+'
*       IFLAG = 4:'NEUTRAL CURRENT, P - PI0'
*       IFLAG = 5:'NEUTRAL CURRENT, N - PI+'
*       IFLAG = 6:'NEUTRAL CURRENT, N - PI0'
*       IFLAG = 7:'NEUTRAL CURRENT, P - PI-'
*
*       IFLAG =11 ... -> charge conjugate
*
***********************************************************************

      IMPLICIT NONE

C     Number of accept-reject tries in the main W, Q2 loop
C     Lives so we can avoid infinite loops (set to something like 10k)
      INTEGER  INFLPCK

      INTEGER IMODE,INLEP

      INTEGER*4 IERR

#include <rscons.h>
#include <neutmodel.h>
#include "necard.h"

      REAL*4 XMASPI2, XMVSPI2

C     Store the computed resonance amplitudes from Rein-Sehgal and FKR
C     Do not need to be 50; should be 31?
      REAL AM3(31),AM1(31),AP1(31),AP3(31),A0P(31),A0M(31)
      REAL AP1C(31),AP3C(31),AM1C(31),AM3C(31),A0PC(31),A0MC(31)
      REAL FFBW(31)

C     Stores all the amplitudes together to make it easier to pass to
C     helper functions
      REAL AM(6,31)
C     Iterator for number of resonances in making AM3C... to AM
      INTEGER NCTR
C     Number of resonances we have to loop to make outgoing particles
C     Depends on pion ejection model used.
C     DRES is just to skip certain resonances we don't need to calculate
      INTEGER MAXRES, DRES

C     These are the 2l+1 and number of resonances
      INTEGER IPP(31),NRR(31)

      INTEGER IFLAG,IB,NR,IP,IBLOCK,IQ2

C     Angles of pion in resonance rest frame
      REAL E,W,Q2
      REAL PHI, THETA

C     Setting up auxilliary kinematic variables
      real W2, ES, tempU

#include <neutcrs.h>

      REAL QV,QV2,QVS,QVS2,A,V,U,V2,U2
      REAL BR,FBWNO,XE,FFV,FFA
      REAL R,RA,RV,RM,RP,T,TA,TV,TP,TM,S,B,C
      REAL FM3C,FM1C,FP3C,FP1C,F0PC,F0MC
      REAL FM3,FM1,FP3,FP1,F0P,F0M
      REAL XXE,RAD,RADR,PQ,PQR,FBW,BRA

      REAL XNU,XNU2,XNUS,XLAM,XLAM2,XMR
      REAL XMLEP,QMIN,ELEPL,APLEPL,APLEPL2,XMLEP2

      REAL GV3,GV1,NFF,CV1,CV2,CV3A,CV3B,CV4,CV5
      REAL CA5

C     Variables for dealing with ejecting the pion in the Npi frame
      COMPLEX ETABW(31)
      COMPLEX ROP3P3,ROM3M3,ROP3P1,ROM1M3,ROP3M1,ROP1M3,ROP1P1,ROM1M1
      COMPLEX RHOTIL,ROTI33,ROTI31,ROT3M1
      REAL BETSH,ALPSH,DELSH
C     The Axy for Rein-Sehgal
      COMPLEX A00,A10,A20,A30,A11,A21,A22,A31,A32,A33
      COMPLEX B11,B21,B22,B31,B32,B33

      COMPLEX RHOCALC
      EXTERNAL RHOCALC

C     Accept-reject for the pion ejection
      REAL YFN

C     For generating outgoing particles
C     Spherical harmonics (in neutcore/rsy20.F etc)
C     Some are real, others complex!
      REAL      RSY00,RSY10,RSY20,RSY30
      EXTERNAL  RSY00,RSY10,RSY20,RSY30
      COMPLEX   RSY11,RSY21,RSY22,RSY31,RSY32,RSY33
      EXTERNAL  RSY11,RSY21,RSY22,RSY31,RSY32,RSY33
C     The normalised rho tildes in Rein-Sehgal pion ejection
      COMPLEX RHOP3P3,RHOP3P1,RHOP3M1

C     Store the signs of the resonances; here be dragons!
      REAL RESSGN(31)
C     See Rein-Sehgal paper, page 114
      DATA RESSGN/ 1., ! P33(1236) P,N (IBLOCK = 0 in rsclfm.F)
     &            -1., ! D13(1525) P
     &            -1., ! D13(1525) N
     &            -1., ! S11(1540) P
     &            -1., ! S11(1540) N (IBLOCK = 4)
     &             1., ! S31(1620) P,N
     &             1., ! S11(1640) P
     &             1., ! S11(1640) N
     &            -1., ! D13(1670) P
     &            -1., ! D13*1670) N
     &             1., ! D15(1680) P
     &             1., ! D15(1680) N
     &             1., ! D33(1730) P,N
     &             1., ! P11(1450) P (IBLOCK 13)
     &             1., ! P11(1450) N
     &             1., ! P33(1640) P,N
     &             1., ! F15(1680) P
     &             1., ! F15(1680) N
     &             1., ! P11(1710) P (IBLOCK 18)
     &             1., ! P11(1710) N
     &             1., ! P13(1740) P
     &             1., ! P13(1740) N
     &            -1., ! P31(1920) P,N
     &            -1., ! F35(1920) P,N
     &             1., ! F37(1950) P,N
     &             1., ! P33(1960) P,N
     &             1., ! F17(1970) P (IBLOCK 26)
     &             1., ! F17(1970) N
     &             1., ! Proton background
     &             1., ! Neutron background
     &             1. / ! J=3/2 background (not active)

C     Store the various constants for given Enu, W, Q2 combination for Rein-Sehgal
      REAL PAR(14)
C     Store the returned amplitudes per resonance
      REAL RETF(6)

      INTEGER I

C     Random numbers
      REAL RLU
      EXTERNAL RLU

C     These are the 2l+1 constants (IPP stores which interaction mode the resonance is relevant to)
      DATA IPP/3,5,5,1,1,1,1,1,5,5,5,5,5,3,3,3,7,7,3,3,3,3,3,7,7,3,7,7,
     1     1,1,3/
C     Harmonic oscillation of the resonance
      DATA NRR/0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
     1     0,0,0/


      IERR = 0

      YFN = 0

C     Set the constants
      XMASPI2= XMASPI * XMASPI
      XMVSPI2= XMVSPI * XMVSPI

      XMLEP = 0.00000
      if ((IMODE.lt.4).or.((IMODE.gt.10).and.(IMODE.lt.14)))
     $     then
         IF (ABS(INLEP).EQ.12) XMLEP=XME
         IF (ABS(INLEP).EQ.14) XMLEP=XMMU
         IF (ABS(INLEP).EQ.16) XMLEP=XMTAU
      endif

      XMLEP2 = XMLEP*XMLEP

C     Set iflag
      IFLAG = IMODE
      IF (IMODE.GT.10) IFLAG = IFLAG -10

      W2=W*W

C     ===================================================================
C     ===================================================================
C     We have now selected W and Q2 for the event from E_nu and the
C     differential cross-section. Now move onto generating outgoing
C     particles, as described by Rein-Sehgal, section 2.3 and appendix A.
C     Have various methods of doing this, see Rein-Sehgal (all
C     implemented here)
C     ===================================================================

C     CWRET: maybe separate subroutine?

C     CWRET: there's a copy of this in rsdcrs too
C     KINEMATICAL QUANTITIES

      XNU = (W2 - XMN2 - Q2) / (2 * XMN)
      XNU2 = XNU**2
      ES = E-XNU

C     Some Rein-Sehgal fundamental constants for a given W, Enu, Q2
C     combination
C     EQ 3.11 roughly
      QV = SQRT(XNU2-Q2)
      QV2 = QV**2
      QVS = QV*XMN/W
      QVS2 = QVS**2
      XNUS = ABS(W2-XMN2+Q2)/(2.*W)
      XLAM = QVS*SQRT(2./OMEG)
      XLAM2 = XLAM**2
      A = ((W+XMN)**2-Q2)/(2.*W)

      U = (E+ES+QV)/(2.*E)
      V = (E+ES-QV)/(2.*E)
      if (V .LT. 0) then
        YFN = -1
        GOTO 72
      endif
      U2 = U**2
      V2 = V**2

C     Moved the U and V replacement to ALPSH, BETSH and DELSH, since
C     this is more clear and what Rein-Sehgal/Ravndal actually do

C     Set the resonance rhos to zero (initialise them)
C     For explanation on what rho is, see Rein-Sehgal appendix
      ROP3P3 = 0.
      ROM3M3 = 0.
      ROP3P1 = 0.
      ROM1M3 = 0.
      ROP1M3 = 0.
      ROP3M1 = 0.
      ROP1P1 = 0.
      ROM1M1 = 0.
      A00 = 0. ! Keep it coming!

C     Although computationally unnecessary doing the larger loop adds
C     virtually no time to generation because the W Q2 accept-reject is the
C     slowest part by far
      MAXRES = 31
      DRES = 1

 200  DO 300 IB = 1, MAXRES, DRES

C     Can skip some calculations depending on MDLSPIEJ

C     If we're doing isotropic can skip all calculations as we won't use
C     them
      if (MDLSPIEJ .EQ. 0) then
        GOTO 300

C     If we're only using the Delta, only need its amplitude
      else if (MDLSPIEJ .EQ. 1 .OR. MDLSPIEJ .EQ. 3) then
        if (IB.NE.1) then
          GOTO 300
        endif

C     If we're using multiple resonances, calculate them
      else if (MDLSPIEJ .EQ. 2) then ! 
C       Pure I3/2
        if (IFLAG .EQ. 1) then
C         Only P33(1232) P33(1640)
          if (IB.NE.1 .and. IB.NE.16) then
            CYCLE
          endif
C       Mixed isospin
        else
          if (IB.NE.1 .and. IB.NE.3 .and. IB.NE.5 .AND. IB.NE.15)
     &      CYCLE
        endif
      endif

C     IBLOCK is just IB-1; has to do with how we've save NRR(IB) and
C     rsclfm.F. Could be fixed if you want :)

      IBLOCK = IB-1
      NR = NRR(IB)
      XMR = XMRR(IB)
      BR = BRR(IB)

C     SETS UP NORMALIZATION CONSTANTS
C     G&S don't use the BW normalizing factors for their model 
C     NOTE: This should only work for G&S form factors only
C     Set branching ratio of Delta to one when G&S form factors
      if(IFFSPI.ne.0) then
        FBWNO = 1
      else
        FBWNO = FBWNOO(IB)
      endif

C     BWNOO is set in rscons.h
      XE = XEE(IB)
      IP = IPP(IB)

C     **********************************
C     FORMFACTORS AND RELATED QUANTITIES

      if (IFFSPI.ne.0) then

C     K.M. Graczyk and J.T. Sobczyk, Phys.Rev.D 77,053001 (2008)
C     Equations are correction in erratum
C     CA5(0) and XMA2 are set by the user

C     Assume CA5 has dipole form
C     CWRET might want to change this!
        CA5=   RCA5ISPI/((1-Q2/XMASPI2)**2)

C       J.A. Nowak, arXiv:0909.3659
C       Vector Form Factors
C       CV1, CV2 are for simplicity, not true form factors
C       Added the proportional factors of the vector form factors GV
C       to CV3,CV4,CV5

C       (PAR) Yargh, the variable Q2 here has a minus sign relative to
C       Q2 in the Nowak/Graczyk & Sobczyk papers

C       CWRET this could be updated to arxiv 1407.5445v1 (2014)
        CV1=  1./((1.-Q2/(4.*XMVSPI2)) * ((1.-Q2/XMVSPI2)**2))
        CV2=  1./((1.-Q2/(4.*XMVSPI2)) * ((1.-Q2/(0.766*XMVSPI2))**2))
        CV3A= 2.13*CV1*((W+XMN)/XMN)
        CV3B= 2.13*CV1*((XMN2-Q2+XMN*W)/(XMN*W))
        CV4= -1.51*CV1*((W2+Q2-XMN2)/(2.*XMN2))
        CV5=  0.48*CV2*((W2-Q2-XMN2)/(2.*XMN2))

        GV3=  (.5/SQRT(3.))*(CV4+CV5+CV3A)
        GV1=  (-.5/SQRT(3.))*(CV4+CV5-CV3B)

        NFF=  SQRT((3.*(GV3)**2)+GV1**2)

C       Note, J.A. Nowak has different form factors that those quoted
C       in G&S (eq. 58). G&S seems to match previous MC better

C     Neutrino-production
      if(NRTYPESPI.eq.0) then
        FFV = 0.5 * (1 - Q2/(XMN+W)**2)**0.5
     &        * (1 - Q2/(4*XMN2))**(-1*NR)
     &        * NFF

C       Axial. Divide by Z, the arbitrary Rein-Sehgal factor, because
C       we multiply by it later

        FFA = 0.5 * SQRT(3.0) * (1 - Q2/(XMN+W)**2)**0.5
     &        * (1 - Q2/(4*XMN2))**(-1*NR)
     &        * (1 - (W2 + Q2 - XMN2)/(8*XMN2))
     &        * CA5 / Z

C     Electro-production
      else
        FFV = 0.5 * (1 - Q2/(XMN+W)**2)**0.5
     &        * (1 - Q2/(4*W2))**(-1*NR/2)
     &        * NFF

        FFA = 0.5 * SQRT(3.0) * (1 - Q2/(XMN+W)**2)**0.5
     &        * (1 - Q2/(4*W2))**(-1*NR/2)
     &        * (1 - (W2 + Q2 - XMN2)/(8*XMN2))
     &        * CA5 / Z

      endif

C     Rein-Seghal form factors
C     Difference between Rein-Seghal form factors are shown in
C     Phys Rev D 77 (053001), equation 12, 13
      else
        if(NRTYPESPI.eq.0) then
           FFV = ((1.-Q2/(4.*XMN2))**(.5-NR/2.))/(1.-Q2/XMVSPI2)**2
           FFA = ((1.-Q2/(4.*XMN2))**(.5-NR/2.))/(1.-Q2/XMASPI2)**2
        else
           FFV = ((1.-Q2/(4.*XMN2))**(.5-NR))/(1.-Q2/XMVSPI2)**2
           FFA = ((1.-Q2/(4.*XMN2))**(.5-NR))/(1.-Q2/XMASPI2)**2
        endif
      endif

C     ATTENTION: NR CHANGED INTO NR/2.
C     EQ 3.11 of Rein-Sehgal
      RV = SQRT(2.)*FFV*QVS*(W+XMN)/(2.*W*A)
      RA = SQRT(2.)*FFA*Z*(W+XMN+NR*OMEG/A)/(6.*W)
      R  = RV
      RP = -(RV+RA)
      RM = -(RV-RA)
      TV = FFV*SQRT(OMEG/2.)/(3.*W)
      TA = FFA*SQRT(OMEG/2.)*Z*QVS/(3.*W*A)
      T  = TV
      TP = -(TV+TA)
      TM = -(TV-TA)
      S  = (FFV*(3.*W*XMN-XMN2+Q2)/(6.*W2))*(-Q2/QVS2)
      B  = FFA*Z*SQRT(OMEG/2.)*(1.+XNUS/A)/(3.*W)
      C  = FFA*Z*(W2-XMN2+NR*OMEG*XNUS/A)/(6.*W*QVS)

C     Write all the fun Rein-Sehgal parameters
      PAR( 1) = RV
      PAR( 2) = RA
      PAR( 3) = R
      PAR( 4) = RP
      PAR( 5) = RM
      PAR( 6) = TV
      PAR( 7) = TA
      PAR( 8) = T
      PAR( 9) = TP
      PAR(10) = TM
      PAR(11) = S
      PAR(12) = B
      PAR(13) = C
      PAR(14) = XLAM

C     CWRET: now calculate the helicity amplitude
      CALL RSCLFM(IFLAG,IBLOCK,PAR,RETF)

C     Naming follows Rein-Sehgal
C     FM3C = helicity amplitude, -3 isospin
      IF (IFLAG.LE.3) THEN
         FM3C = RETF(1)
         FM1C = RETF(2)
         FP3C = RETF(3)
         FP1C = RETF(4)
         F0PC = RETF(5)
         F0MC = RETF(6)
      ELSE
         FM3 = RETF(1)
         FM1 = RETF(2)
         FP3 = RETF(3)
         FP1 = RETF(4)
         F0P = RETF(5)
         F0M = RETF(6)
      ENDIF

C     **********************************

C     STORE RESONANCE AMPLITUDES

C     This is the square root of the elasticity xE
C     Mentioned as (ii) on p88 Rein-Sehgal
      XXE = SQRT(XE)

C     CWRET: modified from RAD and RADR = 0.0001 to 1E-7
C     Eq 2.32 Rein-Sehgal squared
C     The q_pi (W) squared (p87 bottom)
      RAD =((W2-XMN2-XMPI2)**2-4.*XMN2*XMPI2)/(4.*W2)
C     The q_pi (M_res) squared (p87 bottom)
      RADR = ((XMR**2-XMN2-XMPI2)**2-4.*XMN2*XMPI2)/(4.*XMR**2)
C     CWRET: maybe rethrow here instead of approximating?

      IF (RAD .LT. 0-RS_SMNO) then
        YFN = -1
        GOTO 72
      else if (RAD .GT. 0-RS_SMNO .AND. RAD .LT. 0) then
        RAD = 0
      endif

      IF (RADR .LT. 0-RS_SMNO) then
        YFN = -1
        GOTO 72
      else if (RADR .GT. 0-RS_SMNO .AND. RADR .LT. 0) then
        RADR = 0
      endif

C     Now square root them
      PQ = SQRT(RAD)
      PQR = SQRT(RADR)

C     BR is the branching ratio for the resonance, IP is the 2l+1
C     So BRA is Rein-Sehgal 2.32
      BRA = BR*(PQ/PQR)**IP
      FBW = (BRA/(2.*PI))/((W-XMR)**2+(BRA**2)/4.)

C     The Breit-Wigner terms should be complex
      ETABW(IB) = SQRT(BRA/(2.*PI)) 
     &          * (1./CMPLX(W-XMR, BRA/2.))
     &          * (1./SQRT(FBWNO))

      IF(IFLAG.LE.3)  GOTO 991

      AM3(IB) = FM3*XXE
      AM1(IB) = FM1*XXE
      AP1(IB) = FP1*XXE
      AP3(IB) = FP3*XXE
      A0P(IB) = F0P*XXE
      A0M(IB) = F0M*XXE

      AM1C(IB) = AM1(IB)
      AM3C(IB) = AM3(IB)
      AP1C(IB) = AP1(IB)
      AP3C(IB) = AP3(IB)
      A0PC(IB) = A0P(IB)
      A0MC(IB) = A0M(IB)
      GOTO 992

C     If we're doing a CC mode
 991  CONTINUE
      AM3C(IB) = FM3C*XXE
      AM1C(IB) = FM1C*XXE
      AP1C(IB) = FP1C*XXE
      AP3C(IB) = FP3C*XXE
      A0PC(IB) = F0PC*XXE
      A0MC(IB) = F0MC*XXE
 992  CONTINUE

C     End the IBLOCK loop here; we should now have the resonances we
C     need to generate the outgoing pion
 300  CONTINUE


C     Follow the convention of Rein-Sehgal for alpha, beta, gamma on
C     last page of appendix, equation A.23
C     ALPSH = Alpha
C     BETSH = Beta
C     DELSH = Delta (not Kronecker and not resonance symbol: refers to
C     what follows alpha and beta in the paper)
C     Also found in Ravndal (1973, equation A18)
      ALPSH = W*E*SQRT(-Q2 / QVS) * U
      BETSH = W*E*SQRT(-Q2 / QVS) * V
      DELSH = W*E*(XMN / W) * SQRT(2.*U*V)
C     For anti-neutrino the only difference is swapping and sign
C     flipping the alpha and beta terms
      IF (IMODE.GT.10) THEN
        tempu = ALPSH
        ALPSH = -1*BETSH
        BETSH = -1*tempu
C       Ravndal says this last sign flip should happen, but doesn't
C       agree with what NEUT had before. Needs checking
C       DELSH = -1*delsh
      endif

C     Write all the ampltiudes from all the resonances to a 2D array
C     which will be passed to the summation code
      DO 61 NCTR=1,31
         AM(1,NCTR) = AM3C(NCTR)
         AM(2,NCTR) = AM1C(NCTR)
         AM(3,NCTR) = AP3C(NCTR)
         AM(4,NCTR) = AP1C(NCTR)
         AM(5,NCTR) = A0PC(NCTR)
         AM(6,NCTR) = A0MC(NCTR)
 61   CONTINUE

C     Throw the pion according to user input
C     MDLSPIEJ = 0 means isotropic
C     MDLSPIEJ = 1 means Delta only
C     MDLSPIEJ = 2 means full Rein-Sehgal (Appendix A)
C     MDLSPIEJ = 3 means Delta + isotropic, as in old NEUT

      IF (MDLSPIEJ .EQ. 0) THEN
        GOTO 101
C     The usual NEUT method is very similar to Delta only
      ELSE IF (MDLSPIEJ .EQ. 1 .OR. MDLSPIEJ .EQ. 3) THEN
        GOTO 102
      ELSE IF (MDLSPIEJ .EQ. 2) THEN
        GOTO 103
      ELSE
        WRITE(*,*) "WRONG MDLSPIEJ GIVEN"
        WRITE(*,*) "YOU GAVE MDLSPIEJ = ", MDLSPIEJ
        WRITE(*,*) "I NEED:" 
        WRITE(*,*) "0: ISOTROPIC"
        WRITE(*,*) "1: DELTA ONLY"
        WRITE(*,*) "2: FULL REIN-SEHGAL"
        WRITE(*,*) "3: DELTA+FLAT (old NEUT)"
        STOP
      END IF

C     ===========================================================
C     ===========================================================
C     PION EJECTION USING ISOTROPIC DISTRIBUTION
C     ===========================================================
C     Start the pion ejection calculation for ISOTROPIC distribution
 101  CONTINUE

C     Not sure how to do this for a flat distribution
C     Probably need to normalise somehow?
      YFN = 1/(4*PI)

      GOTO 72
C     Now done with this method; have valid Enu, Q2, W, theta, phi
C     Now get the momentum at GOTO 72


C     ===========================================================
C     ===========================================================
C     PION EJECTION USING DELTA ONLY, similar to old NEUT
C     ===========================================================
C     Start the pion ejection calculation for P33(1232) ONLY
C     See Rein-Sehgal 2.46; recommended only for CC1pi+1p channel
 102  CONTINUE

C     For a flat distribution, just make sure integral over
C     of #int(YFN d(costh)dphi)=1
      YFN = 1/(4*PI)

C     If we aren't dealing with a Delta resonance just make the
C     distribution flat in costheta and phi if user has requested it
      IF (MDLSPIEJ .EQ. 3 .AND. SPIDELTA .EQ. 0) GOTO 72

C     rho 3/2 3/2 = delta delta + p delta + delta p + p p
      ROP3P3 = RHOCALC(1, 3, 1, 3, ALPSH, BETSH, DELSH, AM, ETABW)
C     rho -3/2 -3/2
      ROM3M3 = RHOCALC(1,-3, 1,-3, ALPSH, BETSH, DELSH, AM, ETABW)
C     rho +1/2 +1/2
      ROP1P1 = RHOCALC(1, 1, 1, 1, ALPSH, BETSH, DELSH, AM, ETABW)
C     rho -1/2 -1/2
      ROM1M1 = RHOCALC(1,-1, 1,-1, ALPSH, BETSH, DELSH, AM, ETABW)
C     rho +3/2 +1/2
      ROP3P1 = RHOCALC(1, 3, 1, 1, ALPSH, BETSH, DELSH, AM, ETABW)
C     rho -1/2 -3/2
      ROM1M3 = RHOCALC(1,-1, 1,-3, ALPSH, BETSH, DELSH, AM, ETABW)
C     rho +3/2 -1/2
      ROP3M1 = RHOCALC(1, 3, 1,-1, ALPSH, BETSH, DELSH, AM, ETABW)
C     rho +1/2 -3/2
      ROP1M3 = RHOCALC(1, 1, 1,-3, ALPSH, BETSH, DELSH, AM, ETABW)

C     Include the multiplicative factors
      ROP3P3 = ROP3P3*ROP3P3_MULT
      ROM3M3 = ROM3M3*ROM3M3_MULT
      ROP1P1 = ROP1P1*ROP1P1_MULT
      ROM1M1 = ROM1M1*ROM1M1_MULT
      ROP3P1 = ROP3P1*ROP3P1_MULT
      ROM1M3 = ROM1M3*ROM1M3_MULT
      ROP3M1 = ROP3M1*ROP3M1_MULT
      ROP1M3 = ROP1M3*ROP1M3_MULT

C     These are the 33 resonances that will contribute to pion direction (eq 2.47 R-S)
C     RHOTIL is rho tilde, or the total equal isospin sum
      RHOTIL = ROP3P3 + ROP1P1 + ROM1M1 + ROM3M3
C     Rho tilde 3/2, 3/2
      ROTI33 = ROP3P3 + ROM3M3
C     rho tilde 3/2, 1/2
      ROTI31 = ROP3P1 - ROM1M3
C     rho tilde 3/2, -1/2
      ROT3M1 = ROP3M1 + ROP1M3


C     If not, calculate equation 2.46 (remember all rho have been
C     normalised out by rho tilde)
C     Don't use the rescaled amplitudes
      YFN = RSY00(THETA,PHI)*REALPART(RHOTIL)
     $    - (2./ SQRT(5.)) * RSY20(THETA,PHI) *
     $      (REALPART(ROTI33) - 0.5*REALPART(RHOTIL))
     $    + (4./ SQRT(10.)) * REALPART(RSY21(THETA,PHI)) 
     $      * REALPART(ROTI31)
     $    - (4./ SQRT(10.)) * REALPART(RSY22(THETA,PHI)) 
     $      * REALPART(ROT3M1)

C     Make sure the YFN integral over costheta and phi is 1
C     Include the geometric factor
      YFN = YFN/(2*SQRT(PI))
C     Then also the term multiplying RSY00
      YFN = YFN/REALPART(RHOTIL)

      GOTO 72

C     ===========================================================
C     =======================================================
C     PION EJECTION USING MULTIPLE RESONACNES, DEPENDENT ON I (isospin)
C     =======================================================
C     The most complete case, throwing differently depending on
C     interaction mode, taking P33(1232) + nearby resonances into
C     account when generating the outgoing pion
 103  CONTINUE

      IF (IFLAG .EQ. 1) THEN
        GOTO 104 ! I=3/2 (CC1pi+1p) case, P33(1232) + P33(1640)
      ELSE
        GOTO 105 ! I=1/2 case, P33(1232) + others
      END IF

C     ===========================================================
C     ===================================================
C     PION EJECTION USING P33(1232) + P33(1640) FOR I=3/2
C     ===================================================
C     Start the pion ejection calculation for P33(1232) and P33(1640) ONLY
C     Should only be done with I=3/2 channel, i.e. IFLAG = 1
C     Unforunately have to derive this by hand :( Contact Clarence if you need a hand!
C     THE I=3/2 case is covered in Rein-Sehgal eq 2.48
C
C     Resonance 1 is Delta(1232), resonance 16 is P33(1640)
C     You can find this in rsclfm or looking at the masses of resonances
C     is rscons.h
 104  CONTINUE

C     rho 3 3 = delta delta + p delta + delta p + p p
      ROP3P3 = RHOCALC(1, 3, 1, 3, ALPSH, BETSH, DELSH, AM, ETABW)
     &       + RHOCALC(1, 3, 16, 3, ALPSH, BETSH, DELSH, AM, ETABW)
     &       + RHOCALC(16, 3, 1, 3, ALPSH, BETSH, DELSH, AM, ETABW)
     &       + RHOCALC(16, 3, 16, 3, ALPSH, BETSH, DELSH, AM, ETABW)

C     rho -3 -3
      ROM3M3 = RHOCALC(1, -3, 1, -3, ALPSH, BETSH, DELSH, AM, ETABW)
     &       + RHOCALC(1, -3, 16, -3, ALPSH, BETSH, DELSH, AM, ETABW)
     &       + RHOCALC(16, -3, 1, -3, ALPSH, BETSH, DELSH, AM, ETABW)
     &       + RHOCALC(16, -3, 16, -3, ALPSH, BETSH, DELSH, AM, ETABW)

C     rho +1 +1
      ROP1P1 = RHOCALC(1, 1, 1, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &       + RHOCALC(1, 1, 16, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &       + RHOCALC(16, 1, 1, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &       + RHOCALC(16, 1, 16, 1, ALPSH, BETSH, DELSH, AM, ETABW)

C     rho -1 -1
      ROM1M1 = RHOCALC(1, -1, 1, -1, ALPSH, BETSH, DELSH, AM, ETABW)
     &       + RHOCALC(1, -1, 16, -1, ALPSH, BETSH, DELSH, AM, ETABW)
     &       + RHOCALC(16, -1, 1, -1, ALPSH, BETSH, DELSH, AM, ETABW)
     &       + RHOCALC(16, -1, 16, -1, ALPSH, BETSH, DELSH, AM, ETABW)

C     rho +3 +1
      ROP3P1 = RHOCALC(1, 3, 1, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &       + RHOCALC(1, 3, 16, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &       + RHOCALC(16, 3, 1, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &       + RHOCALC(16, 3, 16, 1, ALPSH, BETSH, DELSH, AM, ETABW)

C     rho -1 -3
      ROM1M3 = RHOCALC(1, -1, 1, -3, ALPSH, BETSH, DELSH, AM, ETABW)
     &       + RHOCALC(1, -1, 16, -3, ALPSH, BETSH, DELSH, AM, ETABW)
     &       + RHOCALC(16, -1, 1, -3, ALPSH, BETSH, DELSH, AM, ETABW)
     &       + RHOCALC(16, -1, 16, -3, ALPSH, BETSH, DELSH, AM, ETABW)

C     rho +3 -1
      ROP3M1 = RHOCALC(1, 3, 1, -1, ALPSH, BETSH, DELSH, AM, ETABW)
     &       + RHOCALC(1, 3, 16, -1, ALPSH, BETSH, DELSH, AM, ETABW)
     &       + RHOCALC(16, 3, 1, -1, ALPSH, BETSH, DELSH, AM, ETABW)
     &       + RHOCALC(16, 3, 16, -1, ALPSH, BETSH, DELSH, AM, ETABW)

C     rho +1 -3
      ROP1M3 = RHOCALC(1, 1, 1, -3, ALPSH, BETSH, DELSH, AM, ETABW)
     &       + RHOCALC(1, 1, 16, -3, ALPSH, BETSH, DELSH, AM, ETABW)
     &       + RHOCALC(16, 1, 1, -3, ALPSH, BETSH, DELSH, AM, ETABW)
     &       + RHOCALC(16, 1, 16, -3, ALPSH, BETSH, DELSH, AM, ETABW)

C     Now we construct the same variables as used in the P33(1232) case
C     Except now they have the additional P33(1640) resonance built-in

C     RHOTIL is rho tilde, or the total equal isospin sum (A.17)
      RHOTIL = ROP3P3 + ROP1P1 + ROM1M1 + ROM3M3
C     Rho tilde 3/2, 3/2
      ROTI33 = ROP3P3 + ROM3M3
C     rho tilde 3/2, 1/2
      ROTI31 = ROP3P1 - ROM1M3
C     rho tilde 3/2, -1/2
      ROT3M1 = ROP3M1 + ROP1M3

C     This is eq 2.46 of Rein-Sehgal, RSY are the spherical harmonics

      YFN = RSY00(THETA,PHI)*REALPART(RHOTIL)
     &    - (2./SQRT(5.))*RSY20(THETA,PHI)
     &      * (REALPART(ROTI33)-0.5*REALPART(RHOTIL))
     &    + (4./SQRT(10.))*REALPART(RSY21(THETA,PHI))
     &      * REALPART(ROTI31)
     &    - (4./SQRT(10.))*REALPART(RSY22(THETA,PHI))
     &      * REALPART(ROT3M1)
     &    - (4./SQRT(10.))*IMAGPART(RSY21(THETA,PHI))
     &      * IMAGPART(ROTI31)
     &    + (4./SQRT(10.))*IMAGPART(RSY22(THETA,PHI))
     &      * IMAGPART(ROT3M1)

C     Make sure the YFN integral over costheta and phi is 1
C     Include the geometric factor
      YFN = YFN/(2*SQRT(PI))
C     Then also the term multiplying RSY00
      YFN = YFN/REALPART(RHOTIL)

      GOTO 72

C     ===================================================
C     PION EJECTION USING P33(1232) + others FOR I=1/2
C     ===================================================
C     Start the pion ejection calculation for P33(1232) and others ONLY
C     Should only be done with I=1/2 channel, i.e. IFLAG != 1
C     See Appendix A Rein-Sehgal, last part
C
C     Resonance 1 is Delta(1232), resonance 3 is D13(1525)
C     Resonance 5 is S11(1540), resonance 15 is P11(1450)
C     All on neutron!!
C     You can find this in rsclfm or looking at the masses of resonances
C     is rscons.h
 105  CONTINUE
C     A00 as defined in A.19a in Rein-Sehgal
C     Here's the delta part
      A00 = RHOCALC(1, 3, 1, 3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1, 1, 1, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1,-1, 1,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1,-3, 1,-3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + 2. * (
C     Here's the S11(1540)
     &        RHOCALC(5, 3, 5, 3, ALPSH, BETSH, DELSH, AM, ETABW)
     &      + RHOCALC(5, 1, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &      + RHOCALC(5,-1, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &      + RHOCALC(5,-3, 5,-3, ALPSH, BETSH, DELSH, AM, ETABW)
C     Here's the D13(1525)
     &      + RHOCALC(3, 3, 3, 3, ALPSH, BETSH, DELSH, AM, ETABW)
     &      + RHOCALC(3, 1, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &      + RHOCALC(3,-1, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &      + RHOCALC(3,-3, 3,-3, ALPSH, BETSH, DELSH, AM, ETABW)
C     Here's the P11(1450)
     &      + RHOCALC(15, 3, 15, 3, ALPSH, BETSH, DELSH, AM, ETABW)
     &      + RHOCALC(15, 1, 15, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &      + RHOCALC(15,-1, 15,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &      + RHOCALC(15,-3, 15,-3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the bracket of the 2.*

C     A10 as defined in A.19b in Rein-Sehgal
      A10 = REALPART((2./5.) * SQRT(6.) * (
     &      RHOCALC(1, 3, 3, 3, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    - RHOCALC(1,-3, 3,-3, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    ) ! Close the second bracket
     &    + (2./15.) * SQRT(6.) * (
     &      RHOCALC(1, 1, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW) ! DelD
     &    - RHOCALC(1,-1, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW) ! DelD
     &    ) ! Close the second bracket
     &    + (4./3.) * SQRT(3.) * (
     &      RHOCALC(15, 1, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW) ! PS
     &    - RHOCALC(15,-1, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW) ! PS
     &    ) ! Close the second bracket
     &    + (4./3.) * SQRT(3.) * (
     &      RHOCALC(1, 1, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW) ! DelS
     &    + RHOCALC(1,-1, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW) ! DelS
     &    ) ! Close the second bracket
     &    + (4./3.) * SQRT(6.) * (
     &      RHOCALC(3, 1, 15, 1, ALPSH, BETSH, DELSH, AM, ETABW) ! DP
     &    + RHOCALC(3,-1, 15,-1, ALPSH, BETSH, DELSH, AM, ETABW) ! DP
     &    ) ! Close second bracket
     &    ) ! Close first bracket

C     A20 as defined in A.19c in Rein-Sehgal
C     Note difficult to tell sign on real part in paper!
      A20 = (-2./5.) * SQRT(5.) * (
     &      RHOCALC(1, 3, 1, 3, ALPSH, BETSH, DELSH, AM, ETABW) !DelDel
     &    + RHOCALC(1,-3, 1,-3, ALPSH, BETSH, DELSH, AM, ETABW) !DelDel
     &    - 0.5 * (
     &      RHOCALC(1, 3, 1, 3, ALPSH, BETSH, DELSH, AM, ETABW) !DelDel
     &    + RHOCALC(1, 1, 1, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DelDel
     &    + RHOCALC(1,-1, 1,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DelDel
     &    + RHOCALC(1,-3, 1,-3, ALPSH, BETSH, DELSH, AM, ETABW) !DelDel
     &    ) ! Close the bracket
     &    ) ! Close the bracket
     &    - (4./5.)*SQRT(5.) * (
     &      RHOCALC(3, 3, 3, 3, ALPSH, BETSH, DELSH, AM, ETABW) !DD
     &    + RHOCALC(3,-3, 3,-3, ALPSH, BETSH, DELSH, AM, ETABW) !DD
     &    - 0.5 * (
     &      RHOCALC(3, 3, 3, 3, ALPSH, BETSH, DELSH, AM, ETABW) !DD
     &    + RHOCALC(3, 1, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DD
     &    + RHOCALC(3,-1, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DD
     &    + RHOCALC(3,-3, 3,-3, ALPSH, BETSH, DELSH, AM, ETABW) !DD
     &    ) ! Close the bracket
     &    ) ! Close the bracket
     &    - REALPART(
     &      (4./5.)*SQRT(5.) * (
     &      RHOCALC(1, 1, 15, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DelP
     &    - RHOCALC(1,-1, 15,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DelP
     &    ) ! Close first bracket
     &    + (4./5.)*SQRT(10.) * (
     &      RHOCALC(3, 1, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DS
     &    - RHOCALC(3,-1, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DS
     &    ) ! Close the first bracket
     &    ) ! Close the real part bracket

C     A30 as defined in Rein-Sehgal A.19d
      A30 = REALPART(
     &      (-6./35.) * SQRT(14.) * (
     &      RHOCALC(1, 3, 3, 3, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    - RHOCALC(1,-3, 3,-3, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    ) ! Close second bracket
     &    - (18./35.) * SQRT(14.) * (
     &      RHOCALC(1, 1, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    - RHOCALC(1,-1, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    ) ! Close second bracket
     &    ) ! Close first bracket

C     A11 as defined in Rein-Sehgal A.19e
      A11 = REALPART(
     &      (-4./5.) * (
     &      RHOCALC(1, 3, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    + RHOCALC(1,-1, 3,-3, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    + RHOCALC(1, 1, 3, 3, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    + RHOCALC(1,-3, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    ) ! Close the second bracket
     &    - (4./3.) * SQRT(6.) * (
     &      RHOCALC(15, 1, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW) !PS
     &    + RHOCALC(15,-1, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW) !PS
     &    ) ! Close the second bracket 
     &    - (4./15.) * SQRT(3.) * (
     &      RHOCALC(1, 1, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    + RHOCALC(1,-1, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    ) ! Close the second bracket 
     &    + (2.) * SQRT(3.) * (
     &      RHOCALC(1, 3, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DelS
     &    - RHOCALC(1,-3, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DelS
     &    ) ! Close the second bracket 
     &    + (4.) * (
     &      RHOCALC(3, 3, 15, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DP
     &    - RHOCALC(3,-3, 15,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DP
     &    ) ! Close the second bracket 
     &    - (2.) * SQRT(6.) * ( ! Negative sign? (last line of A.19e?)
     &      RHOCALC(1, 1, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DelS
     &    - RHOCALC(1,-1, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DelS
     &    ) ! Close the second bracket 
     &    + (4.) * SQRT(3.) * (
     &      RHOCALC(3, 1, 15,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DP
     &    - RHOCALC(3,-1, 15, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DP
     &    ) ! Close the second bracket 
     &    ) ! Close the first bracket

C     A21 as defined in A.19f Rein-Sehgal
      A21 = (2./5.) * SQRT(10.) * (
     &      RHOCALC(1, 3, 1, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DelDel
     &    - RHOCALC(1,-1, 1,-3, ALPSH, BETSH, DELSH, AM, ETABW) !DelDel
     &    + 2. * (
     &      RHOCALC(3, 3, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DD
     &    - RHOCALC(3,-1, 3,-3, ALPSH, BETSH, DELSH, AM, ETABW) !DD
     &    ) ! Close the second bracket
     &    ) ! Close the first bracket
     &    + REALPART(
     &      (2./5.) * SQRT(10.) * ( ! Hard to tell if this is positive
     &      RHOCALC(1, 3, 15, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DelP
     &    + RHOCALC(1,-3, 15,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DelP
     &    ) ! Close the second bracket
     &    + (4./5.) * SQRT(5.) * (
     &      RHOCALC(3, 3, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DS
     &    + RHOCALC(3,-3, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DS
     &    ) ! Close the second bracket
     &    - (2./5.) * SQRT(30.) * (
     &      RHOCALC(1, 1, 15,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DelP
     &    + RHOCALC(1,-1, 15, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DelP
     &    ) ! Close the second bracket
     &    - (4./5.) * SQRT(15.) * (
     &      RHOCALC(3, 1, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DS
     &    + RHOCALC(3,-1, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DS
     &    ) ! Close the second bracket
     &    ) ! Close the first bracket

C     A22 as defined in A19g Rein-Sehgal
      A22 = (-2./5.) * SQRT(10.) * (
     &      RHOCALC(1, 3, 1, -1, ALPSH, BETSH, DELSH, AM, ETABW) !DelDel
     &    + RHOCALC(1, 1, 1, -3, ALPSH, BETSH, DELSH, AM, ETABW) !DelDel
     &    + (2.) * (
     &      RHOCALC(3, 3, 3, -1, ALPSH, BETSH, DELSH, AM, ETABW) !DD
     &    + RHOCALC(3, 1, 3, -3, ALPSH, BETSH, DELSH, AM, ETABW) !DD
     &    ) ! Close second bracket 
     &    ) ! Close first bracket
     &    + REALPART( (-4./5.) * SQRT(10.) * ( ! Hard to see sign
     &      RHOCALC(1, 3, 15, -1, ALPSH, BETSH, DELSH, AM, ETABW) !DelP
     &    + RHOCALC(1, -3, 15, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DelP
     &    ) ! Close second bracket
     &    - (8./5.) * SQRT(5.) * (
     &      RHOCALC(3, 3, 5, -1, ALPSH, BETSH, DELSH, AM, ETABW) !DS
     &    - RHOCALC(3, -3, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DS
     &    ) ! Close second bracket
     &    ) ! Close first bracket (REALPART)

C     A31 as defined in A.19h Rein-Sehgal
C     Delta = 1
C     P = 15
C     D = 3
C     S = 5
      A31 = REALPART(
     &      (12./35.) * SQRT(14.) * (
     &      RHOCALC(1, 3, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    + RHOCALC(1,-1, 3,-3, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    + RHOCALC(1, 1, 3, 3, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    + RHOCALC(1,-3, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    ) ! Close second bracket
     &    - (6./35.) * SQRT(42.) * (
     &      RHOCALC(1, 1, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    + RHOCALC(1,-1, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    ) ! Close second bracket
     &    ) ! Close first bracket

C     A32 as defined in A.19i
      A32 = REALPART(
     &      (-12./35.) * SQRT(35.) * (
     &      RHOCALC(1, 3, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    - RHOCALC(1, 1, 3,-3, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    + RHOCALC(1,-1, 3, 3, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    - RHOCALC(1,-3, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    ) ! Close second bracket
     &    ) ! Close first bracket

C     A33 as defined in A.19j (superfluous?)
      A33 = 0.

C     B11 as defined in A.19k
C     Delta = 1
C     P = 15
C     D = 3
C     S = 5
      B11 = IMAGPART(
     &      (4./5.) * (
     &      RHOCALC(1, 3, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    + RHOCALC(1,-1, 3,-3, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    - RHOCALC(1, 1, 3, 3, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    - RHOCALC(1,-3, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    ) ! Close the second bracket
     &    + (4./3.) * SQRT(6.) * (
     &      RHOCALC(15, 1, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW) !PS
     &    - RHOCALC(15,-1, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW) !PS
     &    ) ! Close the second bracket
     &    + (4./15.) * SQRT(3.) * (
     &      RHOCALC(1, 1, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    - RHOCALC(1,-1, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    ) ! Close the second bracket
     &    - 2. * SQRT(2.) * (
     &      RHOCALC(1, 3, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DelS
     &    + RHOCALC(1,-3, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DelS
     &    ) ! Close the second bracket
     &    - 4. * (
     &      RHOCALC(3, 3, 15, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DP
     &    + RHOCALC(3,-3, 15,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DP
     &    ) ! Close the second bracket
     &    - 2. * SQRT(6.) * (
     &      RHOCALC(1, 1, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DelS
     &    + RHOCALC(1,-1, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DelS
     &    ) ! Close the second bracket
     &    - 4. * SQRT(3.) * (
     &      RHOCALC(3, 1, 15,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DP
     &    + RHOCALC(3,-1, 15, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DP
     &    ) ! Close the second bracket
     &    ) ! Close the first bracket

C     B21 as defined in A.19l
C     Delta = 1
C     P = 15
C     D = 3
C     S = 5
      B21 = IMAGPART(
     &    - (2./5.) * SQRT(10.) * (
     &      RHOCALC(1, 3, 15, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DelP
     &    - RHOCALC(1,-3, 15,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DelP
     &    ) ! Close the second bracket
     &    - (4./5.) * SQRT(5.) * (
     &      RHOCALC(3, 3, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DS
     &    - RHOCALC(3,-3, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DS
     &    ) ! Close the second bracket
     &    + (2./5.) * SQRT(30.) * (
     &      RHOCALC(1, 1, 15,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DelP
     &    - RHOCALC(1,-1, 15, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DelP
     &    ) ! Close the second bracket
     &    + (4./5.) * SQRT(15.) * (
     &      RHOCALC(3, 1, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DS
     &    - RHOCALC(3,-1, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DS
     &    ) ! Close the second bracket
     &    ) ! Close the first bracket

C     B22 as defined in A.19m
C     Delta = 1
C     P = 15
C     D = 3
C     S = 5
      B22 = IMAGPART(
     &      (4./5.) * SQRT(10.) * (
     &      RHOCALC(1, 3, 15,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DelP
     &    + RHOCALC(1,-3, 15, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DelP
     &    ) ! Close second bracket
     &    + (8./5.) * SQRT(5.) * (
     &      RHOCALC(3, 3, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DS
     &    + RHOCALC(5,-3, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DS
     &    ) ! Close second bracket
     &    ) ! Close first bracket

C     B31 as defined in A.19n
C     Delta = 1
C     P = 15
C     D = 3
C     S = 5
      B31 = IMAGPART(
     &      (-12./35.) * SQRT(14.) * (
     &      RHOCALC(1, 3, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    + RHOCALC(1,-1, 3,-3, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    - RHOCALC(1, 1, 3, 3, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    - RHOCALC(1,-3, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    ) ! Close second bracket
     &    + (6./35.) * SQRT(42.) * (
     &      RHOCALC(1, 1, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    - RHOCALC(1,-1, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    ) ! Close second bracket
     &    ) ! Close first bracket

C     B32 as defined in A.19o
      B32 = IMAGPART(
     &      (12./35.) * SQRT(35.) * (
     &      RHOCALC(1, 3, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    - RHOCALC(1, 1, 3,-3, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    - RHOCALC(1,-1, 3, 3, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    + RHOCALC(1,-3, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW) !DelD
     &    ) ! Close second bracket
     &    ) ! Close first bracket
      
C     B33 as defined in A.19p (superfluous?)
      B33 = 0.


C     The theoretical description can be found in Rein-Sehgal, eq 2.46
C     Rein-Sehgal also provide theory for calculating outgoing pion
C     kinematics using more resonances, which is a future project!

C     W(THETA,PHI) as defined in A.18 Rein-Sehgal
      YFN = (
     &      A00*RSY00(THETA,PHI) + A10*RSY10(THETA,PHI)
     &    + A20*RSY20(THETA,PHI) + A30*RSY30(THETA,PHI)

     &    + A11*REALPART(RSY11(THETA,PHI)) 
     &    + A21*REALPART(RSY21(THETA,PHI))
     &    + A22*REALPART(RSY22(THETA,PHI)) 
     &    + A31*REALPART(RSY31(THETA,PHI))
     &    + A32*REALPART(RSY32(THETA,PHI)) 
     &    + A33*REALPART(RSY33(THETA,PHI))

     &    + B11*IMAGPART(RSY11(THETA,PHI))
     &    + B21*IMAGPART(RSY21(THETA,PHI))
     &    + B22*IMAGPART(RSY22(THETA,PHI))
     &    + B31*IMAGPART(RSY31(THETA,PHI))
     &    + B32*IMAGPART(RSY32(THETA,PHI))
     &    + B33*IMAGPART(RSY33(THETA,PHI))
     &    )

C     Make sure the YFN integral over costheta and phi is 1
C     Include the geometric factor
      YFN = YFN/(2*SQRT(PI))
C     Then also the term multiplying RSY00
      YFN = YFN/A00

      GOTO 72

 72   CONTINUE

      if (YFN.lt.0) THEN
        YFN=0
      endif

      END
