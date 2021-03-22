***********************************************************************
*     ---------------------------------------------------
      SUBROUTINE RSDT33(IMODE,INLEP,E,RETPMU,RETPPI,RETPNU,IERR)
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
*       retpmu(3)    : momentum of outgoing lepton
*       retppi(3)    : momentum of outgoing pion
*       retpnu(3)    : momentum of outgoing nucleon
*
*     ( creation date and author )
*       1993.Dec. ; for Super-Kamioka by Y.Hayato
*                   (Based on Rein and Sehgal's code)
*       1996.Nov. ; for Avoid Infinite Loop
*                   Check Cross-section after Fermi-momentum
*                   has been fixed.
*       1997.Aug. ; debug: q^2 minimum is larger than the allowed
*                          region(pneucm is wrong.)
*       2007.Dec. ; add lepton mass effects by G.Mitsuka
*       2013.10    ; P.Rodrigues & A.Bercellieadd 
*                              new form facor for Rein&Sehgal 1pi
*       2016.Jan. ; C.Wret: Bug-fix Q2 throw to not depend on W
*                           Added significant commenting for future
*                           Changed some limits to be more clear (pneucm replaced)
*                           Added lepton mass-effects in generating pion
*                           previously only accounted for in rsdcrs.F
*       2016.Oct. ; C.Wret: Add other resonance to pion ejection
*                           Eject pion relative the mode
*                           Add possibility of changing rho in ejection
*
*     ( comment )
*       Reference:
*                 Original paper:
*                   Rein and Sehgal,Ann. of Phys. 133(1981),79-153
*                 Rein           ,Z.Phys.C 35(1987),43-64
*                 FKR formalism:
*                 Feynman et al. ,Phys.Rev.D3(1971),2706-
*                 Lepton mass-effect
*                 K.S.Kuzmin et al., hep-ph/0312107
*                 Berger and Sehgal, hep-ph/0709.4378
*                 Graczyk & Sobczyk Form factors
C                 K.M. Graczyk and J.T. Sobczyk, Phys.Rev.D 77,053001 (2008)
C                    Equations are correction in erratum
C                 J.A. Nowak, arXiv:0909.3659
C                    Vector Form Factors
C                    CV1, CV2 are for simplicity, not true form factors
C                    Added the proportional factors of the vector form factors
C                    GV to CV3,CV4,CV5
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
C     The returned momentum arrays of the muon, pion and nucleon
C     after the pion ejection is completed
      REAL RETPMU(3),RETPPI(3),RETPNU(3)

      REAL SINTH,COSTH,SINPHI,COSPHI
      INTEGER*4 IDUMMY
      INTEGER*4 IERR

#include <rscons.h>
#include <neutmodel.h>
#include "necard.h"

      REAL*4 XMASPI2, XMVSPI2

      INTEGER IRSISED

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

      INTEGER IFLAG,IB,NR,IP,IBLOCK,IQ2,RES33

      REAL WMIN,E,QMAX,SMAX,W2,W,X,Q2,ES,DQ2
      INTEGER NQ

      INTEGER NTRIES

C     CWRET Multiplicative safety factor for accept-reject ceiling
      REAL SAFETY
C     CWRET added this for a W_max which is smax-xmlep or WMAX
C     WMAX parameter is set in rscons.h and is 2.0 GeV (Rein-Sehgal default)
      REAL W_max
C     CWRET added these as alternative to calculated Q2 limits
      REAL ERES,ELEP,PRES,PLEP,COS1,COS2,COSLEP

C     CWRET added these for checking Q2 is within limits
      REAL QMAXFOUND,QMINFOUND

#include <neutcrs.h>

      REAL QV,QV2,QVS,QVS2,A,V,U,V2,U2,AK,ALFA,ALFA2
      REAL BETA,BETA2,GAMA,GAMA2,FAKT
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
C     The calculator function for rho, at the bottom of this
      COMPLEX RHOCALC
      EXTERNAL RHOCALC
C     Helpers
C     REAL rhop3p3r, rhop3p1r, rhop3m1r
C     REAL rhop3p3i, rhop3p1i, rhop3m1i
C     REAL rsy21r, rsy22r, rsy21i, rsy22i

C     Accept-reject for the pion ejection
      REAL QTRPI,YRAND,MAXDST,YFN

C     For generating outgoing particles
C     Spherical harmonics (in neutcore/rsy20.F etc)
C     Some are real, others complex!
      REAL      RSY00,RSY10,RSY20,RSY30
      EXTERNAL  RSY00,RSY10,RSY20,RSY30
      COMPLEX   RSY11,RSY21,RSY22,RSY31,RSY32,RSY33
      EXTERNAL  RSY11,RSY21,RSY22,RSY31,RSY32,RSY33
C     Angles of pion in resonance rest frame
      REAL PHI,THETA,PHILEP
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

      REAL RNDSIG,DSIGMX,DSIGMXP,DSIGMXM,DSIG3,DSIG3P,DSIG3M
C     CWRET added these to agree with rsdcrs and have lepton mass-effect
C     accounted for in vector generation
      Real El
      Real Q_lab, Q2_lab, Q, nu_lab
      Real jl(4), jl_n1, jl_n2, Norm, C_L, C_R, C_S, Qs, nus, ng
      Real tempC_L
      Real Ei, Ef, nu
      Real B_BRS,C_BRS,Ctmp
      Integer LambdaIT

C     Used for accept-reject
      REAL TMPSIG,SIGMAX
      INTEGER Lambda(2), iOrg

      REAL TMPPMU(3),TMPPPI(3),TMPPNU(3)

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

C     Lepton mass effects are considered
C
C     iOrg = 0 considers lepton mass effects
C     iOrg = 1 considers original R-S model without mass effects (eg. NC, HighEnu)
C     Hardcode this to have lepton mass effects on for now
      iOrg = 0                  ! CC does need mass correction
      if (ABS(inlep).eq.12.and.E.gt.1.5 ) iOrg = 1 ! High E
      if (ABS(inlep).eq.14.and.E.gt.100.) iOrg = 1 ! High E
      if (ABS(inlep).eq.16.and.E.gt.200.) iOrg = 1 ! High E
      if (IFLAG.GE.4) iOrg = 1  ! NC

      if (IMODE.LT.10) then     ! neutrino
         Lambda(1) =  1
         Lambda(2) = -1
      else                      ! anti-neutrino
         Lambda(1) = -1
         Lambda(2) =  1
      endif

C     Set up some fundamental kinematic quantities
C     Minimum W is simply nucleon mass + pion mass
      WMIN = XMN+XMPI
      SMAX = (2 * XMN * E) + XMN2

C     CWRET:
C     We have a maximum W set by Enu (or equally SMAX) 
C     Previous NEUT had W_max = WMAX, which isn't allowed for low E_nu
C     W_max = MIN(SQRT(SMAX)+XMLEP2-2*E*XMLEP-2*XMLEP*XMN), WMAX)
      W_max = MIN(SQRT(SMAX)-XMLEP, WMAX)

C     Get the maximum cross-section for this given Enu (bottom of this file)
C     Safety = Multiplicative safety factor that multiplies the maximum
C     cross-section found
C     If found cross-section is larger than SIGMAX need SAFETY bigger
      SAFETY = 1.5
      SIGMAX = 0.
C     How many times we've tried to set the ceiling
      NTRIES = 0
      CALL GETMAXCRS_RS(E,W_max,XMLEP,SAFETY,IMODE,IORG,SIGMAX)
 111  NTRIES = NTRIES+1
      IF (NTRIES .gt. 1) THEN
        write(*,*) "Increasing ceiling to ", SIGMAX*(NTRIES+1.0)/3.0
        write(*,*) "*********************"
        SIGMAX = SIGMAX*(NTRIES+1.0)/3.0
      ENDIF

C     Set the minimum and maximum Q2 limits for throwing Q2 in
C     accept-reject
C     Important: Set ONLY so it's dependent on E_nu. If this is set to
C     be dependent on W it will significantly bias and correlate W and
C     Q2 in the kinematics.
C     This was a subtle but significant bug which caused a lot of pain!

C     QMIN can't be zero because Rein-Sehgal divides by Q2. Make very
C     small instead

      QMIN = -1E-5
      QMAX = -4*E*E

C     Try these later! (might make our space slightly smaller)
C     QMIN = m_mu * (m_mu - 2 E_nu)
C     QMAX = m_mu * m_mu - 4 * E_nu * E_nu

C     This is here to check that we don't get into an infinite loop
C     when throwing W, Q2 and RNDSIG to unphysical values or values
C     that don't pass the accept-reject method.
      INFLPCK=0

C     CWRET: 11 is where we end up when "faulty" event:
C            Either 
C            a) It isn't physically allowed (e.g. |costh| > 1)
C            b) We've been rejected in the accept-reject method.
C            Keep the INFLPCK variable to make sure we don't run into
C            infinite loops.
 11   CONTINUE

      INFLPCK=INFLPCK+1
C     CWRET: If we throw more than 1E8 times something has probably gone
C     wrong. It could be our accept-reject ceiling is ridiculously high,
C     or some kinematic problem. Stop program and unfortunately someone
C     might have to debug (sorry!)
      IF (INFLPCK .GT. 1E8) THEN
         write(*,*) "RSDT33: TOO MANY TRY/STOP..."
         IERR = 1
         write(*,*) "E = ", E, "W = ", W, "Wmin = ", WMIN,
     &   "Wmax = ", W_max, "W_max = ", W_max, "IMODE = ", IMODE,
     &   "INLEP (PDG)=", INLEP, "SIGMAX =", SIGMAX, "RNDSIG = ", RNDSIG,
     &   "DSIGMX = ", DSIGMX, "RLU(IRSISED) = ", RLU(IRSISED),
     &   "I live in src/neutcore/rsdt33.F, ca line 250"
         RETURN
      ENDIF

C     Start the throwing for accept-reject method
C     Throw W within given range (our x-value)
C     CWRET changed WMAX -> W_max
C     Since Wmax is limited by the incoming neutrino energy
      W = RLU(IRSISED)*(W_max - WMIN)+WMIN
      W2 = W*W

C     Generate a random cross-section up to SIGMAX
C     This is our random z-step in the accept-reject method
      RNDSIG = RLU(IRSISED)*SIGMAX

C     Q2 is negative by definition, and in Q2max and Q2min range
C     Idea here was to find Qmin and Qmax from parameter space search
C     Abadoned because it's faster to simply do the maths

C     Save the kinematically allowed Q2 limits in QMINFOUND and QMAXFOUND
C     This depends on Enu, W and the lepton mass. We'll reject Q2 if
C     it's outside of these limits
      CALL RES_Q2_LIMITS(SMAX, W, XMLEP, QMINFOUND, QMAXFOUND)

C     Choose a random Q2 in the range determined by Enu,
C     this is our "y-value" in the accept-reject method
      Q2 = QMIN + RLU(IRSISED)*(QMAX-QMIN)

C     CWRET:
C     Limit the allowed Q2 to between QMAX and QMIN
C     N.B.  This is significantly different to throwing Q2 in the
C           kinematically allowed Q2 range. Explanation below:
C
C           The kinematically allowed Q2 range is determined by both Enu
C           and W. Since W is included and thrown in our accept-reject 
C           as a variable this enforces a correlation between W and Q2
C           IN WHICH THE ACCEPT-REJECT METHOD IS INVALID!
C
C           This caused a bug present in NEUT since time immemorial, 
C           so please don't change this unless 100% sure!
      IF ((Q2 .LT. QMAXFOUND) .OR. (Q2 .GT. QMINFOUND)) GOTO 11

C     CWRET: Added in case of numerical accuracy problems
      IF ((Q2 .LT. QMAX) .OR. (Q2 .GT. QMIN)) THEN
        WRITE(*,*) "Invalid Q2: Q2 > Q2MAX or Q2 < Q2MIN"
        GOTO 11
      END IF

      IF (QMIN .LE. QMAX) THEN
        WRITE(*,*) "Invalid Q2LIM: Q2min >= Q2max"
        GOTO 11
      END IF

C     Rethrow W, Q2 and maximum sigma if unphysical:
C     Same checks are made in finding maximum cross-section
C     These shouldn't fire often since the Q2 check above (using the
C     kinematic Q2 limits to evaluate random Q2) implicitly include
C     kinematic checks

      ERES=(XMN2+W2-Q2)/(2*XMN)
C     Make sure resonance has at least rest-mass
      IF (ERES .LT. W) GOTO 11

      ELEP=E+XMN-ERES
C     Make sure lepton has at least rest-mass
      IF (ELEP .LT. XMLEP) GOTO 11

      PLEP=sqrt(ELEP**2-XMLEP2)
C     Make sure lepton has momentum
      IF (PLEP .LE. 0) GOTO 11

      PRES=sqrt(ERES**2-W2)
C     Make sure resonance has momentum
      IF (PRES .LE. 0) GOTO 11

      CosLep=(q2-XMLEP2+2*E*ELEP)/(2*E*PLEP)

C     Angle between lepton and neutrino
      IF (ABS(CosLep) .GT. 1.) GOTO 11

C     Angle between Resonance and Muon
      COS1 = (XMLEP2+W**2+2*ELEP*ERES - XMN2 - 2*E*XMN)
     $      /(2*PLEP*PRES)

      IF (ABS(COS1) .GT. 1.) GOTO 11

C     Angle between Delta and Neutrino
      COS2 = (XMLEP2 + XMN2 - 2*ELEP*XMN - W**2 + 2*E*ERES)
     &  /(2*E*PRES)

      IF (ABS(COS2) .GT. 1.) GOTO 11

      IF ( ((2 * XMN * E) * (2 * XMN * E + XMN2 - W2) /
     $     (2 * XMN * E + XMN2)).LT.(-Q2) ) GOTO 11

      IF ((W2 - XMN2 - Q2) / (2 * XMN) .GT. E) GOTO 11

C     Note: rsdcrs.F wants the negative sign Q2.
C     Call rsdcrs twice because lepton mass effect needs two lambdas
      CALL RSDCRS(IMODE,iOrg,Lambda(1),XMLEP,E,Q2,W,DSIGMXP,DSIG3P)
      CALL RSDCRS(IMODE,iOrg,Lambda(2),XMLEP,E,Q2,W,DSIGMXM,DSIG3M)

C     The differential cross-section for this W and Q2
C     Add the two lambdas together to create total cross-section
      DSIGMX = DSIGMXP+DSIGMXM
C     Differential cross-section from delta alone
      DSIG3  = DSIG3P +DSIG3M

C     Finally decide on keeping the Q2 and W combination based on a
C     simple accept-reject method
C     This can be expanded and significantly sped up by adaptive
C     accept-reject sampling!
      IF (RNDSIG .GT. DSIGMX) GOTO 11
C     Found a valid Enu, W, Q2 combination by now

C     Now find if this event is Delta resonance
C     Used later in the pion angle generation
C     If RES33 .EQ. 0, do isotropic flat, else do pion+nucleon ejection
C     using distribution from the Delta resonance only
      RES33 = 0
C     A random thrown xsec
      TMPSIG = RLU(IRSISED) * DSIGMX
C     Accept-reject on the Delta-only cross-section
      IF (TMPSIG.LT.DSIG3) THEN
         RES33 = 1
      ENDIF

C     Check that DSIGMX (calculated x-sec for thrown Enu, Q2, W)
C     is less than the SIGMAX used in the accept-reject method.
C     If this is the case, need to increase SIGMAX or else method fails
C     i.e. change SIGMAX = SIGMAX * (bigger number than currently)
      IF (DSIGMX .GT. SIGMAX) THEN
        write(*,*) "Calculated x-sec > accept-reject ceiling"
        write(*,*) "E = ", E, "W = ", W, "Q2 = ", Q2,
     &             "x-sec = ", DSIGMX, "SIGMAX = ", SIGMAX
        GOTO 111
      END IF

c     enu_com = e
c     q2_com = q2
c     w_com = w



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
      X = -Q2 / (2 * XMN * XNU)
      XNU2 = XNU**2
      ES = E-XNU

C     Rethrow if something weird has happened to the neutrino energy
      IF (XNU.GE.E) GOTO 11

C     This would be a disaster, stop the routine!
      IF (XNU2.LE.Q2) THEN
         WRITE(*,*) "ERROR"
         IERR=1
         RETURN
      ENDIF

C     Some Rein-Sehgal fundamental constants for a given W, Enu, Q2
C     combination
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
C     CWRET: this is kind of weird?
C     IF(V.LT.0.) V = 0.
      IF(V.LT.0.) GOTO 11
      U2 = U**2
      V2 = V**2
      AK = (W2-XMN2)/(2.*XMN)
      FAKT = W*PI/(W2-XMN2)

C     CWRET: Added to include lepton mass effects in getting the Delta
C     amplitude

C     Need to have both lambda contributions (lambda=+1,-1)
C     ---------------------------
C     --- Lepton mass effects ---

C     Lepton energy
C     El = E+XMN-(XMN2+W**2-q2)/(2*XMN)
C     Ei = (2.*XMN*E+Q2-XMLEP2)/(2.*W)
C     Ef = (2.*XMN*El-Q2-XMLEP2)/(2.*W)
C     nu = Ei - Ef
C     Q = QVS

C     Lab frame
C     If we somehow have generated a strange lepton mass, rethrow!
C     IF (El**2-XMLEP2 .LT. 0) THEN
C        GOTO 11
C     ELSE IF (El**2 - XMLEP2 .LT. 0) THEN
C        GOTO 11
C     ELSE
C        Plep = SQRT(El**2-XMLEP2)
C     ENDIF

C     Theta between incoming neutrino and outgoing lepton
C     theta = (Q2-XMLEP2+2*E*El)/(2*E*Plep)

C     If theta is forward+numerical accuracy, set to forward
C     if (theta .GT. 1.0 .AND. theta .LT. 1.0+1E-7) then
C     CWRET changed here too
C     IF (ABS(THETA) .GT. 1) THEN
C       GOTO 11
C     ENDIF

C     Some fundamental variables for the lepton mass correction
C     nu_lab = E - El
C     Q_lab = SQRT(nu_lab**2-Q2)
C     Q2_lab = Q_lab**2

C     Set up the case where we don't have a mass effect
      C_L = u
      C_R = v
      C_S = sqrt(2.*u*v)

C     Neutrino or anti neutrino: U->V and V->U

C     CWRET: this could be an error; I don't think there should be a
C     negative sign in front. Doesn't matter for calculation though,
C     since we only use these squqared
C     Upd: sign shouldn't matter since U->-V and V->-U and we only do
C     UV, UU or VV multiplications
C     CWRET: reference eq 2.18 Rein-Sehgal.

C     CWRET: replace U with C_L, V with C_R, UV with C_S
C     If anti-neturino reverse C_L and - C_R
      IF (IMODE.GT.10) THEN
         tempC_L = C_L
         C_L = -C_R
         C_R = -tempC_L
      END IF

C     Qs = 1.
C     nus = 1.

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

C     CWRET: Should start Lambda loop here for mass effect in pion accept-reject
C     Need to use both -1 and +1 Lambda
C     Only true for mass-effect case; for no mass-effect we simply get
C     same amplitude twice, so divide the amplitude by two at end of
C     loop
C     
C     This loop includes lepton mass-effects for the outgoing pion
C     direction, which was previously not done
C     It was only done for the calculation of the cross-section

C     DO 90 LambdaIT=1,2

C     Berger-Sehgal mass-effect
C     if (iOrg.ne.1) then
C     Leptonic current in RRF
C     jl_0, jl_x, i*jl_y, jl_z
C        jl(1) = (XMN-El-Lambda(LambdaIT)*Plep)/W
C        jl(2) = (Plep-Lambda(LambdaIT)*E)/Q_lab
C        jl(3) = -Lambda(LambdaIT)
C        jl(4) = ((E+Lambda(LambdaIT)*Plep)*(XMN-El)+Plep
C    $          *(Lambda(LambdaIT)*E+2.*E*theta-Plep))/(W*Q_lab)

C        jl_n1 = SQRT(1.-Lambda(LambdaIT)*theta)*
C    &          SQRT(E*(El-Lambda(LambdaIT)*Plep))
C        jl_n2 = SQRT(1.+Lambda(LambdaIT)*theta)*
C    &          SQRT(E*(El-Lambda(LambdaIT)*Plep))

C        jl(1) = jl_n1*jl(1)
C        jl(2) = jl_n2*jl(2)
C        jl(3) = jl_n2*jl(3)
C        jl(4) = jl_n1*jl(4)

C     Coefficients C_L, C_R, C_S and Qs, nus
C        Norm = SQRT(2.)*E/Q_lab*SQRT(-Q2)
C        C_L = (jl(2)+jl(3))/(SQRT(2.)*Norm)
C        C_R = (jl(3)-jl(2))/(SQRT(2.)*Norm)
C        C_S = SQRT(ABS(jl(1)**2-jl(4)**2))/Norm

C       Check that this makes sense
C        if (C_S.eq.0) then
C           iOrg = 1
C           goto 200
C        endif

C        Qs  = jl(1)/(Norm)*SQRT(-Q2)/C_S
C        nus = jl(4)/(Norm)*SQRT(-Q2)/C_S
C     endif


C     Previous code forgot lepton mass-effects here. Lepton mass-effects
C     dampen the resonance amplitude, and should be applied accordingly
C     to rsdcrs.F, otherwise we make an inconsistent calculation
C     (cross-section made with mass-effect, outgoing without
C     mass-effect!)

C     Consider only Delta resonance, so set IB = 1 (see e.g. rsclfm.F for
C     explanation of IBLOCK and IB numbering convention)
C     Adding more resonances? Need to loop IB moar!

C     Although computationally unnecessary doing the larger loop adds
C     virtually no time to generation because the W Q2 accept-reject is the
C     slowest part by far
      MAXRES = 31
      DRES = 1

 200  DO 300 IB = 1, MAXRES, DRES

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

C     CWRET: lepton mass effects weren't taken into account in 
C     generating the outgoing pion. Now uses modified U, V, UV, S, B
C     and C accordingly (C_L, C_R, C_S, S, B, C)

C     ---------------------------
C     --- Lepton mass effects ---
C     ng = (W+XMN)**2-Q2

C     New definitions of S, B, C for mass-effect
C     S = (Qs*Q-nus*nu)*(3.*W*XMN+Q2-XMN2)/(6.*Q2_lab*XMN2)*FFV
C     B = 2.*SQRT(OMEG/2.)*(Qs+2.*W*nus*Q/ng)/(6.*W*Q)*Z*FFA
C     C = ((Qs*Q-nus*nu)*(1./3.+2.*W*nu/ng)
C    &    +nus*(2.*W/3.+(NR*OMEG/3.+Q2)*2.*W/ng))/(2.*W*Q)*Z*FFA

C     Pion-pole contribution (according Berger and Sehgal ArXiV:0709.4387)
C     Correction is added to B_BRS(2M_Na should be M_Na in the reference)
C     when using mass-effect
C     B_BRS = (Qs*nu-nus*Q)/(3.*W*Q*(XMPI**2-Q2))
C    &        *sqrt(OMEG/2.)*(nu+2.*Q**2*W/ng)*Z*FFA
C     C_BRS = (Qs*nu-nus*Q)/(XMPI**2-Q2)
C    &        *(1./3.+Q2/ng+NR*OMEG/(3.*ng))*Z*FFA

C     B = B + B_BRS
C     C = C + C_BRS
C     ---------------------------
C     ---------------------------

C     Original dynamical form factors without mass-effect
C     if (iOrg .EQ. 1) then
C     CWRET: should probably be XMN2, not W2, in denominator, RS 3.11
C       S = (FFV*(3.*W*XMN-XMN2+Q2)/(6.*W2))*(-Q2/QVS2)
C       B = FFA*Z*SQRT(OMEG/2.)*(1.+XNUS/A)/(3.*W)
C     CWRET same here in denominator
C       C = FFA*Z*(W2-XMN2+NR*OMEG*XNUS/A)/(6.*W*QVS)
C     endif

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
      IF(RAD .LE. 0) RAD = RS_SMNO
      IF(RADR .LE. 0) RADR = RS_SMNO
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
C     If we're doing a NC mode
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
      ALPSH = W*E*SQRT(-Q2 / QVS) * C_L
      BETSH = W*E*SQRT(-Q2 / QVS) * C_R
      DELSH = W*E*(XMN / W) * C_S

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

C     restore C_L, C_R
C     CWRET: this could be an error; I don't think there should be a
C     negative sign in front. Doesn't matter for calculation though,
C     since we only use these squared
      IF (IMODE.GT.10) THEN
         tempC_L = C_L
         C_L = -C_R
         C_R = -tempC_L
      END IF

C     Throw the pion according to user input
C     MDLSPIEJ = 0 means isotropic
C     MDLSPIEJ = 1 means Delta only
C     MDLSPIEJ = 2 means full Rein-Sehgal (Appendix A)

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
        WRITE(*,*) "I NEED 0 (ISOTROPIC)",
     &             "1 (DELTA ONLY)", 
     &             "2 FULL REIN-SEHGAL",
     &             "3 (DELTA+FLAT old NEUT)"
        STOP
      END IF

C     ===========================================================
C     ===========================================================
C     PION EJECTION USING ISOTROPIC DISTRIBUTION
C     ===========================================================
C     Start the pion ejection calculation for ISOTROPIC distribution
 101  CONTINUE
C     WRITE(*,*) "EJECTING PION USING ISOTROPY"
      THETA = ACOS(1.-RLU(IRSISED)*2.)
      PHI   = RLU(IRSISED) * PI * 2.
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
C     WRITE(*,*) "EJECTING PION USING P33(1232) ONLY"

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

C     These are the 33 resonances that will contribute to pion direction (eq 2.47 R-S)
C     RHOTIL is rho tilde, or the total equal isospin sum
      RHOTIL = ROP3P3 + ROP1P1 + ROM1M1 + ROM3M3
C     Rho tilde 3/2, 3/2
      ROTI33 = ROP3P3 + ROM3M3
C     rho tilde 3/2, 1/2
      ROTI31 = ROP3P1 - ROM1M3
C     rho tilde 3/2, -1/2
      ROT3M1 = ROP3M1 + ROP1M3
C     This is eq 2.46 of Rein-Sehgal, RSY are the spherical harmonics
C     Normalise relative rho tilde, as Rein-Sehgal also do
      RHOP3P3 = ROTI33 / RHOTIL
      RHOP3P1 = ROTI31 / RHOTIL
      RHOP3M1 = ROT3M1 / RHOTIL

C     The maximum for the accept-reject ceiling, eq 2.46 maximised
C     (spherical harmonics maximised)
      MAXDST = SQRT(1./(4.*PI))
     &       - (2./SQRT(5.)) * (1./4.) * SQRT(5./PI) * (-1.) *
     &         (RHOP3P3 - 0.5)
     &       + (4./SQRT(10.)) * (-1./2.) * SQRT(15./(2.*PI)) *0.5 *
     &         RHOP3P1
     &       - (4./SQRT(10.)) * (1./4.) * SQRT(15./(2.*PI)) * -1.*
     &         RHOP3M1

C     If this theta, phi, yrand combination fails the accept-reject,
C     rethrow them!
C     Note: GOTO 5x will always rethrow; 50 for this simple Delta only model)
 50   CONTINUE
C     Throw CosTheta between -1 and +1
      THETA = ACOS(1.-RLU(IRSISED)*2.)
      PHI   = RLU(IRSISED) * PI * 2.

C     If we aren't dealing with a Delta resonance just make the
C     distribution flat in costheta and phi if user has requested it
      IF (MDLSPIEJ .EQ. 3 .AND. RES33 .EQ. 0) GOTO 72

C     If not, calculate equation 2.46 (remember all rho have been
C     normalised out by rho tilde)
      YFN = RSY00(THETA,PHI)
     $    - (2./ SQRT(5.)) * RSY20(THETA,PHI) *
     $      (RHOP3P3 - 0.5)
     $    + (4./ SQRT(10.)) * RSY21(THETA,PHI) * RHOP3P1
     $    - (4./ SQRT(10.)) * RSY22(THETA,PHI) * RHOP3M1

C     Do the accept-reject on the YFN distribution
      YRAND = RLU(IRSISED) * MAXDST
C     Accept-reject on YFN, rethrow
      IF (YRAND .GT. YFN) THEN
        GOTO 50
      ELSE 
        GOTO 72
      END IF
C     Now done with this method; have valid Enu, Q2, W, theta, phi
C     Now get the momentum at GOTO 72


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

C     RHOTIL is rho tilde, or the total equal isospin sum
      RHOTIL = ROP3P3 + ROP1P1 + ROM1M1 + ROM3M3
C     Rho tilde 3/2, 3/2
      ROTI33 = ROP3P3 + ROM3M3
C     rho tilde 3/2, 1/2
      ROTI31 = ROP3P1 - ROM1M3
C     rho tilde 3/2, -1/2
      ROT3M1 = ROP3M1 + ROP1M3

C     This is eq 2.46 of Rein-Sehgal, RSY are the spherical harmonics
C     Normalise relative rho tilde, as Rein-Sehgal also do
      RHOP3P3 = ROTI33 / RHOTIL
      RHOP3P1 = ROTI31 / RHOTIL
      RHOP3M1 = ROT3M1 / RHOTIL

C     Get the maximum distribution for any theta and phi
C     Needed for the accept-reject ceiling
      MAXDST = SQRT(1./(4.*PI))
     &  + (2./SQRT(5.))*SQRT(5./(16.*PI))*(3.-1.)
     &    *(REALPART(RHOP3P3) - 0.5)
     &  + (4./SQRT(10.))*SQRT(15./(8.*PI))
     &    *REALPART(RHOP3P1)
     &  + (4./SQRT(10.))*SQRT(15./(32.*PI))
     &    *REALPART(RHOP3M1)
     &  + (4./SQRT(10.))*SQRT(15./(8.*PI))
     &    *IMAGPART(RHOP3P1)
     &  + (4./SQRT(10.))*SQRT(15./(32.*PI))
     &    *IMAGPART(RHOP3M1)

C     Even though both dists (maxdst and yfn) end up with this 
C     prefactor, include for completeness
      MAXDST = MAXDST/SQRT(PI)

C     Now have the amplitudes: do the accept-reject on the pion
 51   CONTINUE
C     WRITE(*,*) "Trying again adv"

C     Throw in costheta and phi
      YRAND = RLU(IRSISED) * MAXDST
      THETA = ACOS(1.-RLU(IRSISED)*2.)
      PHI   = RLU(IRSISED) * PI * 2.

      YFN = RSY00(THETA,PHI)
     &    - (2./SQRT(5.))*RSY20(THETA,PHI)
     &      *(REALPART(RHOP3P3)-0.5)
     &    + (4./SQRT(10.))*REALPART(RSY21(THETA,PHI))
     &      *REALPART(RHOP3P1)
     &    - (4./SQRT(10.))*REALPART(RSY22(THETA,PHI))
     &      *REALPART(RHOP3M1)
     &    - (4./SQRT(10.))*IMAGPART(RSY21(THETA,PHI))
     &      *IMAGPART(RHOP3P1)
     &    + (4./SQRT(10.))*IMAGPART(RSY22(THETA,PHI))
     &      *IMAGPART(RHOP3M1)

      YFN = YFN/SQRT(PI)

      IF (YRAND .GT. YFN) THEN
        GOTO 51
      ELSE 
        GOTO 72
      END IF
C     Now done with this method; have valid Enu, Q2, W, theta, phi
C     Now get the momentum at GOTO 72


C     ===================================================
C     PION EJECTION USING P33(1232) + others FOR I=1/2
C     ===================================================
C     Start the pion ejection calculation for P33(1232) and others ONLY
C     Should only be done with I=1/2 channel, i.e. IFLAG != 1
C     See Appendix A Rein-Sehgal, last part
 105  CONTINUE
C     WRITE(*,*) "EJECTING PION USING P33(1232) + others"
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
     &      RHOCALC(1, 3, 1, 3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(1,-3, 1,-3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket
     &    + (2./15.) * SQRT(6.) * (
     &      RHOCALC(1, 1, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(1,-1, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket
     &    + (4./3.) * SQRT(3.) * (
     &      RHOCALC(15, 1, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(15,-1, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket
     &    + (4./3.) * SQRT(3.) * (
     &      RHOCALC(1, 1, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1,-1, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket
     &    + (4./3.) * SQRT(6.) * (
     &      RHOCALC(3, 1, 15, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(3,-1, 15,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close second bracket
     &    ) ! Close first bracket

C     A20 as defined in A.19c in Rein-Sehgal
      A20 = (-2./5.) * SQRT(5.) * (
     &      RHOCALC(1, 3, 1, 3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1,-3, 1,-3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - 0.5 * (
     &      RHOCALC(1, 3, 1, 3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1, 1, 1, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1,-1, 1,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1,-3, 1,-3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the bracket
     &    ) ! Close the bracket
     &    - (4./5.)*SQRT(5.) * (
     &      RHOCALC(3, 3, 3, 3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(3,-3, 3,-3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - 0.5 * (
     &      RHOCALC(3, 3, 3, 3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(3, 1, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(3,-1, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(3,-3, 3,-3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the bracket
     &    ) ! Close the bracket
     &    - REALPART(
     &      (4./5.)*SQRT(5.) * (
     &      RHOCALC(1, 1, 15, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(1,-1, 15,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close first bracket
     &    + (4./5.)*SQRT(10.) * (
     &      RHOCALC(3, 1, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(3,-1, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the first bracket
     &    ) ! Close the real part bracket

C     A30 as defined in Rein-Sehgal A.19d
      A30 = REALPART(
     &      (-6./35.) * SQRT(14.) * (
     &      RHOCALC(1, 3, 3, 3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(1,-3, 3,-3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close second bracket
     &    - (18./35.) * SQRT(14.) * (
     &      RHOCALC(1, 1, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(1,-1, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close second bracket
     &    ) ! Close first bracket

C     A11 as defined in Rein-Sehgal A.19e
      A11 = REALPART(
     &      (-4./5.) * (
     &      RHOCALC(1, 3, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1,-1, 3,-3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1, 1, 3, 3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1,-3, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket
     &    - (4./3.) * SQRT(6.) * (
     &      RHOCALC(15, 1, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(15,-1, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket 
     &    - (4./15.) * SQRT(3.) * (
     &      RHOCALC(1, 1, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1,-1, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket 
     &    + (2.) * SQRT(3.) * (
     &      RHOCALC(1, 3, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(1,-3, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket 
     &    + (4.) * (
     &      RHOCALC(3, 3, 15, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(3,-3, 15,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket 
     &    + (2.) * SQRT(6.) * (
     &      RHOCALC(1, 1, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(1,-1, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket 
     &    + (4.) * SQRT(3.) * (
     &      RHOCALC(3, 1, 15,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(3,-1, 15, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket 
     &    ) ! Close the first bracket

C     A21 as defined in A.19f Rein-Sehgal
      A21 = (2./5.) * SQRT(10.) * (
     &      RHOCALC(1, 3, 1, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(1,-1, 1,-3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + 2. * (
     &      RHOCALC(3, 3, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(3,-1, 3,-3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket
     &    ) ! Close the first bracket
     &    + REALPART(
     &      (2./5.) * SQRT(10.) * (
     &      RHOCALC(1, 3, 15, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1,-3, 15,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket
     &    + (4./5.) * SQRT(5.) * (
     &      RHOCALC(3, 3, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(3,-3, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket
     &    - (2./5.) * SQRT(30.) * (
     &      RHOCALC(1, 1, 15,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1,-1, 15, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket
     &    - (4./5.) * SQRT(15.) * (
     &      RHOCALC(3, 1, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(3,-1, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket
     &    ) ! Close the first bracket

C     A22 as defined in A19g Rein-Sehgal
      A22 = (-2./5.) * SQRT(10.) * (
     &      RHOCALC(1, 3, 1, -1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1, 1, 1, -3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + (2.) * (
     &      RHOCALC(3, 3, 3, -1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(3, 1, 3, -3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close second bracket 
     &    ) ! Close first bracket
     &    + REALPART( (-4./5.) * SQRT(10.) * (
     &      RHOCALC(1, 3, 15, -1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1, -3, 15, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close second bracket
     &    - (8./5.) * SQRT(5.) * (
     &      RHOCALC(3, 3, 5, -1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(3, -3, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close second bracket
     &    ) ! Close first bracket (REALPART)

C     A31 as defined in A.19h Rein-Sehgal
C     Delta = 1
C     P = 15
C     D = 3
C     S = 5
      A31 = REALPART(
     &      (12./35.) * SQRT(14.) * (
     &      RHOCALC(1, 3, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1,-1, 3,-3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1, 1, 3, 3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1,-3, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close second bracket
     &    - (6./35.) * SQRT(42.) * (
     &      RHOCALC(1, 3, 1,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1,-1, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close second bracket
     &    ) ! Close first bracket

C     A32 as defined in A.19i
      A32 = REALPART(
     &      (-12./35.) * SQRT(35.) * (
     &      RHOCALC(1, 3, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(1, 1, 3,-3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1,-1, 3, 3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(1,-3, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW)
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
     &      RHOCALC(1, 3, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1,-1, 3,-3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(1, 1, 3, 3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(1,-3, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket
     &    + (4./3.) * SQRT(6.) * (
     &      RHOCALC(15, 1, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(15,-1, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket
     &    + (4./15.) * SQRT(3.) * (
     &      RHOCALC(1, 1, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(1,-1, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket
     &    - 2. * SQRT(2.) * (
     &      RHOCALC(1, 3, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1,-3, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket
     &    - 4. * (
     &      RHOCALC(3, 3, 15, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(3,-3, 15,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket
     &    - 2. * SQRT(6.) * (
     &      RHOCALC(1, 1, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1,-1, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket
     &    - 4. * SQRT(3.) * (
     &      RHOCALC(3, 1, 15,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(3,-1, 15, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket
     &    ) ! Close the first bracket

C     B21 as defined in A.19l
C     Delta = 1
C     P = 15
C     D = 3
C     S = 5
      B21 = IMAGPART(
     &    - (2./5.) * SQRT(10.) * (
     &      RHOCALC(1, 3, 15, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(1,-3, 15,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket
     &    - (4./5.) * SQRT(5.) * (
     &      RHOCALC(3, 3, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(3,-3, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket
     &    + (2./5.) * SQRT(30.) * (
     &      RHOCALC(1, 1, 15,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(1,-1, 15, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket
     &    + (4./5.) * SQRT(15.) * (
     &      RHOCALC(3, 1, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(3,-1, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close the second bracket
     &    ) ! Close the first bracket

C     B22 as defined in A.19m
C     Delta = 1
C     P = 15
C     D = 3
C     S = 5
      B22 = IMAGPART(
     &      (4./5.) * SQRT(10.) * (
     &      RHOCALC(1, 3, 15,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1,-3, 15, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close second bracket
     &    + (8./5.) * SQRT(5.) * (
     &      RHOCALC(3, 3, 5,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(5,-3, 5, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close second bracket
     &    ) ! Close first bracket

C     B31 as defined in A.19n
C     Delta = 1
C     P = 15
C     D = 3
C     S = 5
      B31 = IMAGPART(
     &      (-12./35.) * SQRT(14.) * (
     &      RHOCALC(1, 3, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1,-1, 3,-3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(1, 1, 3, 3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(1,-3, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close second bracket
     &    + (6./35.) * SQRT(42.) * (
     &      RHOCALC(1, 1, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(1,-1, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close second bracket
     &    ) ! Close first bracket

C     B32 as defined in A.19o
      B32 = IMAGPART(
     &      (12./35.) * SQRT(35.) * (
     &      RHOCALC(1, 3, 3,-1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(1, 1, 3,-3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    - RHOCALC(1,-1, 3, 3, ALPSH, BETSH, DELSH, AM, ETABW)
     &    + RHOCALC(1,-3, 3, 1, ALPSH, BETSH, DELSH, AM, ETABW)
     &    ) ! Close second bracket
     &    ) ! Close first bracket
      
C     B33 as defined in A.19p (superfluous?)
      B33 = 0.


C     This is where we goto if rethrow
 52   CONTINUE
C     The theoretical description can be found in Rein-Sehgal, eq 2.46
C     Rein-Sehgal also provide theory for calculating outgoing pion
C     kinematics using more resonances, which is a future project!

C     Throw costheta between -1 and +1
      THETA = ACOS(1.-RLU(IRSISED)*2.)
      PHI   = RLU(IRSISED) * PI * 2.

C     W(THETA,PHI) as defined in A.18 Rein-Sehgal
      YFN = (1./3.) * (1./SQRT(PI)) * (
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

C     W(THETA,PHI) as defined in A.18 Rein-Sehgal, maximum
      MAXDST = (1./3.) * (1./SQRT(PI)) * (
     &      A00*SQRT(1./(4.*PI)) 
     &    + A10*(1./2.)*SQRT(3./PI)
     &    + A20*(1./4.)*SQRT(5./PI)*(3.-1.)
     &    + A30*(1./4.)*SQRT(7./PI)*(5.-3.)

     &    + A11*(1./2.)*SQRT(3./(2.*PI))
     &    + A21*(1./2.)*SQRT(15./(2.*PI))
     &    + A22*(1./4.)*SQRT(15./(2.*PI))
     &    + A31*(1./8.)*SQRT(21./PI)*(5.-1.)
     &    + A32*(1./4.)*SQRT(105./(2.*PI))
     &    + A33*(1./8.)*SQRT(35./PI)

     &    + B11*(1./2.)*SQRT(3./(2.*PI))
     &    + B21*(1./2.)*SQRT(15./(2.*PI))
     &    + B22*(1./4.)*SQRT(15./(2.*PI))
     &    + B31*(1./8.)*SQRT(21./PI)*(5.-1.)
     &    + B32*(1./4.)*SQRT(105./(2.*PI))
     &    + B33*(1./8.)*SQRT(35./PI)
     &    )

C     Throw a random number up to MAXDST for the accept-reject
      YRAND = RLU(IRSISED) * MAXDST

C     Accept-reject on YFN, rethrow
      IF (YRAND.GT.YFN) THEN
        GOTO 52
      ELSE 
        GOTO 72
      END IF

C     Now done with this method; have valid Enu, Q2, W, theta, phi
C     Now get the momentum at GOTO 72


C     ======================================
C     END DIFFERENT THROWS FOR THETA AND PHI
C     NOW BOOST TO LAB FRAME AND GET MOMENTA
C     ======================================

C     This is where we want to end up once 
C     we've got a good theta and phi for any model
 72   CONTINUE

c     th_com = theta
c     phi_com = phi

C     Boost to get the muon, pion and nucleon 3-momenta
      CALL RSLZBT(E,Q2,W,XMLEP,THETA,PHI,TMPPMU,TMPPPI,TMPPNU)

      SINTH = 0.
      COSTH = 1.
C     Note, PI = 3.1415... defined in rscons.h
C     Means the number, not the particle
      PHILEP = PI*2.*RLU(IDUMMY)
      SINPHI = sin(PHILEP)
      COSPHI = cos(PHILEP)

C     Rotate the muon
      CALL RS3DRT(TMPPMU,SINTH,COSTH,SINPHI,COSPHI)
C     Rotate the pion
      CALL RS3DRT(TMPPPI,SINTH,COSTH,SINPHI,COSPHI)
C     Rotate the nucleon
      CALL RS3DRT(TMPPNU,SINTH,COSTH,SINPHI,COSPHI)

C     Write the outgoing particle kinematics, and we're done!
      DO 80 I=1,3
         RETPMU(I) = TMPPMU(I)
         RETPPI(I) = TMPPPI(I)
         RETPNU(I) = TMPPNU(I)
 80   CONTINUE

      END


***********************************************************************
*     ---------------------------------------------------
      SUBROUTINE GETMAXCRS_RS(E, W_max, XMLEP, SAFETY, 
     &                        IMODE, IORG, SIGMAX)
*     ---------------------------------------------------
*
*     ( purpose )
*       To calculate and return the maximum cross-section in W and Q2
*       for given Enu to establish the accept-reject ceiling
*
*       Previously this was done by assuming a Delta(1232) and scanning
*       in Q2 space with a huge safety factor. Now a quick W and Q2 scan
*       is done with a small safety factor
*
*     ( input )
*       SMAX      : Maximum S for interaction
*       WMAX      : Maximum W for interaction
*       SAFETY    : Multiplicative 'safety factor' for accept-reject ceiling
*
*     ( output )
*       MAXCRS    : Maximum cross-section * ceiling
*
*     ( creation date and author )
*       2016.Jan  : C. Wret, simplifying Rein-Sehgal subroutine by cutting into
*                   separate subroutines. Also speeding up
*                   accept-reject method
*
***********************************************************************

      IMPLICIT NONE

      REAL E, W_MAX, SAFETY
      INTEGER IMODE, IORG

#include <rscons.h>
      REAL XMLEP, SMAX

      REAL WMINit, WMAXit, dW, W, Wpeak
      REAL Q2MIN, Q2MAX, dQ2, Q2, Q2peak
      INTEGER INW, NW, IQ2, NQ2

C     x-sec related variables
      REAL SIGMAX
      REAL TMPSIG, TMPSIGP, TMPSIGM, JUNK

C     Physical variable checks
      REAL ERES,ELEP,PRES,PLEP,COS1,COS2,COSLEP

C     Q2min and Q2max for a given W
      REAL Q2MINkin, Q2MAXkin

      INTEGER Lambda(2)

      if (IMODE.LT.10) then     ! neutrino
         Lambda(1) =  1
         Lambda(2) = -1
      else                      ! anti-neutrino
         Lambda(1) = -1
         Lambda(2) =  1
      endif

      SIGMAX = 0.

      SMAX = (2 * XMN * E) + XMN**2
      W_max = MIN(SQRT(SMAX)-XMLEP, WMAX)

C     Step around in W space a bit to find maximum (1.10 - 1.25)
C     Peak in W is very likely to be around Delta(1232)
C     Can be around secondary resonances too, but should be comparable
C     in cross-section to Delta(1232), and our safety factor will save
C     us in that case
      WMINit = 1.17
      WMAXit = 1.23

C     Might run into phase space where 1.15 is not allowed
      IF (W_max .LT. WMINit) THEN
        WMINit = XMN+XMPI + 1E-4
        WMAXit = W_max - 1E-4
      END IF
      IF (W_max .LT. WMAXit) THEN
        WMAXit = W_max - 1E-4
      ENDIF
        

C     Number of W steps (only need a fairly coarse scan)
      NW = 5

C     Width of step
      dW = (WMAXit - WMINit)/NW

      DO 500 INW = 0, NW
        W = WMINit + dW*INW

C       Leave these hard-coded for now
C       We're only interested in the peak, so look in low Q2 region
        Q2MIN = -0.02
        Q2MAX = -0.10

C       Get the physically allowed Q2 limits and check that the above
C       hard-coded limits are allowed. If not, set QMIN and QMAX
        CALL RES_Q2_LIMITS(SMAX, W, XMLEP, Q2MINkin, Q2MAXkin)

C       q2min>q2minkin
        IF (q2minkin .ne. q2maxkin) then
          IF (Q2MIN .GT. Q2MINkin) THEN
            Q2MIN = Q2MINkin-1E-4
          END IF
C       q2maxkin>q2max
          IF (Q2MAX .LT. Q2MAXkin) THEN
            Q2MAX = Q2MAXkin+1E-4
          END IF
        end if

C       number of Q2 slices to consider
        NQ2   = 8

C       Define DQ2 to be negative
        DQ2   = (Q2MAX-Q2MIN)/NQ2

        DO 5 IQ2 = 0,NQ2
        Q2 = Q2MIN + DQ2*IQ2

C         Do some physical event checks
C         These should all be OK since the Q2 limits take these into
C         account. However, numerical accuracy problems might still
C         occur so it's safer to keep this included

C         Now call the cross-section calculation and save x-sec
C         Note: JUNK is the cross-section from a Delta alone
C               and is not used here
        CALL RSDCRS(IMODE,iOrg,Lambda(1),XMLEP,E,Q2,W,TMPSIGP,JUNK)
        CALL RSDCRS(IMODE,iOrg,Lambda(2),XMLEP,E,Q2,W,TMPSIGM,JUNK)
        TMPSIG = TMPSIGP+TMPSIGM

C         If cross-section at given Enu, Q2 and W
C         is less than previous iteration we've found the maximum
        IF (TMPSIG.LT.SIGMAX) GOTO 5

C         Only update SIGMAX when TMPSIG > previous TMPSIG; i.e. we're
C         moving up a slope
        SIGMAX = TMPSIG
        Wpeak = W
        Q2peak = Q2

C     Finished looping Q2
 5       CONTINUE

C     Finished looping W
 500  CONTINUE

C     Finally multiply found maximum by safety factor
      SIGMAX = SIGMAX * SAFETY

C     Unfortunately this hack is needed for the case when the Q2 limits
C     are so close to each other that a cross-section isn't calculated
      IF (SIGMAX .LT. RS_SMNO) SIGMAX = RS_SMNO

      END ! end of GETMAXCRS_RS subroutine

************************************************************************
*     ---------------------------------
      COMPLEX FUNCTION RHOCALC( RES1, J1,
     &                          RES2, J2,
     &                          ALPSH, BETSH, DELSH,
     &                          AM, ETABW)
*     ---------------------------------
C
C     (Purpose)
C       Calculate the rho^{r,s}_{m,m'} for two interfering resonances
C       r,s  = two resonance (r = s is allowed), e.g. Delta and S11(1450)
C       m,m' = two isospins to be summed over
C
C       See Rein-Sehgal 1981 p.150 and p.151
C
C     (Input)
C      RES1  :  Number of the first resonance 
C               (see numbering in rsclfm.F and add 1!)
C      AMRES1:  Amplitudes of first resonance
C          J1:  Isospin of first resonance
C      ETABW1:  Eta Breit-Wigner of first resonance
C
C      RES2  :  Number of the second resonance 
C               (see numbering in rsclfm.F and add 1!)
C      AMRES2:  Amplitudes of second resonance
C          J2:  Isospin of second resonance
C      ETABW2:  Eta Breit-Wigner of second resonance
C
C      ALPSH: Alpha as defined in equation A.23 Rein-Sehgal
C      BETSH: Beta as defined in equation A.23 Rein-Sehgal
C      DELSH: Delta as defined in equation A.23 Rein-Sehgal
C
C     (Output)    
C      The rho^{RES1, RES2}_{J1,J2} of two resonances
C
C     (Creation Date and Author)
C       2016.Oct C. Wret: Add in to allow using multiple interfering 
C                         resonances in ejecting the pion in the
C                         resonance rest-frame
C
************************************************************************
      IMPLICIT NONE
C     Resonance number (found in rsclfm.F)
      INTEGER RES1,RES2
C     Isospin (J = 3 means 3/2)
      INTEGER J1,J2

      INTEGER NCTR

C     Array that holds the Breit-Wigner etas
      COMPLEX ETABW(31)

      REAL ALPSH,BETSH, DELSH

C     The alpha, beta and delta as defined in equation A.23 Rein-Sehgal

C     Actually only need to pass 2*6 helicity amplitudes
C     Holds the amplitudes
      REAL AM(6,31)
C     AMRES1(1) = AM3C for resonance 1
C     AMRES1(2) = AM1C for resonance 1
C     AMRES1(3) = AP3C for resonance 1
C     AMRES1(4) = AP1C for resonance 1
C     AMRES1(5) = A0PC for resonance 1
C     AMRES1(6) = A0MC for resonance 1
      REAL AM3C(31), AM1C(31), AP3C(31), AP1C(31), A0PC(31), A0MC(31)

C     The idea behind this function is that the way we perform 
C     the summation only depends on the isospin of the resonances, not
C     on their amplitudes or Breit-Wigners

C     Some error checking
      IF (J1.NE.3 .AND. J1.NE.1 .AND. J1.NE.-1 .AND. J1.NE.-3) THEN
        WRITE(*,*) "ERROR, J1 IS NOT 3, 1, -1, -3; J1 = ", J1
        WRITE(*,*) "EXITING"
        STOP
      END IF

      IF (J2.NE.3 .AND. J2.NE.1 .AND. J2.NE.-1 .AND. J2.NE.-3) THEN
        WRITE(*,*) "ERROR, J2 IS NOT 3, 1, -1, -3; J2 = ", J2
        WRITE(*,*) "EXITING"
        STOP
      END IF

      IF (RES1.GT.31 .OR. RES1.LT.1) THEN
        WRITE(*,*) "INVALID RES1 GIVEN; RES1 = ", RES1
        WRITE(*,*) "EXITING"
        STOP
      END IF

      IF (RES2.GT.31 .OR. RES2.LT.1) THEN
        WRITE(*,*) "INVALID RES2 GIVEN; RES2 = ", RES2
        WRITE(*,*) "EXITING"
        STOP
      END IF

C     However, rho^{r,s}_{m,m'} does depend on the amplitudes and the 
C     Breit-Wigner's, so we need to pass those to the function
C     But we don't know which ones to pass until we've got the isospins,
C     So pass all of the resonance amplitudes

      DO 99 NCTR=1,31
      AM3C(NCTR) = AM(1, NCTR)
      AM1C(NCTR) = AM(2, NCTR)
      AP3C(NCTR) = AM(3, NCTR)
      AP1C(NCTR) = AM(4, NCTR)
      A0PC(NCTR) = AM(5, NCTR)
      A0MC(NCTR) = AM(6, NCTR)
 99   CONTINUE     

C     Set to zero just to be sure!
      RHOCALC = 0.

C     Will always have the Breit-Wigner factor
      RHOCALC = 0.5*CONJG(ETABW(RES1))*ETABW(RES2)


C     3/2 xx
      IF (J1 .EQ. 3) THEN

C     3/2 3/2 
        IF (J2 .EQ. 3) THEN
          RHOCALC = RHOCALC*BETSH*BETSH*AM3C(RES1)*AM3C(RES2)
          RETURN
C     3/2 1/2
        ELSE IF (J2 .EQ. 1) THEN
          RHOCALC = RHOCALC*BETSH*DELSH*AM3C(RES1)*A0PC(RES2)
          RETURN
C     3/2 -1/2
        ELSE IF (J2 .EQ. -1) THEN
          RHOCALC = RHOCALC*(-1.)*BETSH*ALPSH*AM3C(RES1)*AP1C(RES2)
          RETURN
C     3/2 -3/2 -- THIS DOESN'T SEEM TO HAPPEN?
C       ELSE IF (J2 .EQ. -3) THEN

        END IF ! end the if on J2

C     1/2 xx
      ELSE IF (J1 .EQ. 1) THEN

C     1/2 3/2
        IF (J2 .EQ. 3) THEN
          RHOCALC = RHOCALC*BETSH*DELSH*A0PC(RES1)*AM3C(RES2)
          RETURN
C     1/2 1/2
        ELSE IF (J2 .EQ. 1) THEN
          RHOCALC = RHOCALC*(BETSH*BETSH*AM1C(RES1)*AM1C(RES2)
     &            + DELSH*DELSH*A0PC(RES1)*A0PC(RES2))
          RETURN
C     1/2 -1/2
        ELSE IF (J2 .EQ. -1) THEN
          RHOCALC = RHOCALC*(BETSH*DELSH*AM1C(RES1)*A0MC(RES2)
     &            - DELSH*ALPSH*A0PC(RES1)*AP1C(RES2))
          RETURN
C     1/2 -3/2
        ELSE IF (J2 .EQ. -3) THEN
          RHOCALC = RHOCALC*(-1.)*BETSH*ALPSH*AM1C(RES1)*AP3C(RES2)
          RETURN

        END IF ! end the if on J2


C     -1/2 xx
      ELSE IF (J1 .EQ. -1) THEN

C     -1/2 3/2
        IF (J2 .EQ. 3) THEN
          RHOCALC = RHOCALC*BETSH*BETSH*AM1C(RES1)*AM3C(RES2)
          RETURN
C     -1/2 1/2
        ELSE IF (J2 .EQ. 1) THEN
          RHOCALC = RHOCALC*(DELSH*BETSH*A0MC(RES1)*AM1C(RES2)
     &            - ALPSH*DELSH*AP1C(RES1)*A0PC(RES2))
          RETURN
C     -1/2 -1/2
        ELSE IF (J2 .EQ. -1) THEN
          RHOCALC = RHOCALC*(ALPSH*ALPSH*AP1C(RES1)*AP1C(RES2)
     &            + DELSH*DELSH*A0MC(RES1)*A0MC(RES2))
          RETURN
C     -1/2 -3/2
        ELSE IF (J2 .EQ. -3) THEN
          RHOCALC = RHOCALC*(-1.)*ALPSH*DELSH*A0MC(RES1)*AP3C(RES2)
          RETURN

        END IF ! end the if on J2


C     -3/2 xx
      ELSE IF (J1 .EQ. -3) THEN

C     -3/2 3/2 THIS DOESN'T SEEM TO HAPPEN
C       IF (J2 .EQ. 3) THEN
C     -3/2 1/2
        IF (J2 .EQ. 1) THEN
          RHOCALC = RHOCALC*(-1.)*ALPSH*BETSH*AP3C(RES1)*AM1C(RES2)
          RETURN
C     -3/2 -1/2
        ELSE IF (J2 .EQ. -1) THEN
          RHOCALC = RHOCALC*(-1.)*ALPSH*DELSH*AP3C(RES1)*A0MC(RES2)
          RETURN
C     -3/2 -3/2
        ELSE IF (J2 .EQ. -3) THEN
          RHOCALC = RHOCALC*ALPSH*ALPSH*AP3C(RES1)*AP3C(RES2)
          RETURN

        END IF ! end the if on J2

      END IF ! end the if on J1

      RETURN

      END
