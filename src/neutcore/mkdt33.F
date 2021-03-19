***********************************************************************
*     ---------------------------------------------------
      SUBROUTINE MKDT33(IMODE,INLEP,E,RETPMU,RETPPI,RETPNU)
*     ---------------------------------------------------
*
*     ( purpose )
*       calculate kinematics of final state particles for Minoo
*       Kabirnezhad's single pion model

*       Largely inspired by rsdt33.F, the Rein-Sehgal equivalent
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
*       Oct 2016: Clarence Wret and Minoo Kabirnezhad
*
*       ----------- IMODE definition ----------
*  Neutrino
*      IMODE = 1:   CC 'P - PI+'
*      IMODE = 2:   CC 'P - PI0'
*      IMODE = 3:   CC 'N - PI+'
*      IMODE = 4:   NC 'P - PI0'
*      IMODE = 5:   NC 'N - PI+'
*      IMODE = 6:   NC 'N - PI0'
*      IMODE = 7:   NC 'P - PI-'

*  anti-Neutrino
*      IMODE = 11:    CC 'N - PI-'
*      IMODE = 12:    CC 'N - PI0'
*      IMODE = 13:    CC 'P - PI-' 
*      IMODE = 14:    NC 'P - PI0'
*      IMODE = 15:    NC 'N - PI+' 
*      IMODE = 16:    NC 'N - PI0'
*      IMODE = 17:    NC 'P - PI-' 
*
***********************************************************************

      IMPLICIT NONE

      INTEGER  INFLPCK

      INTEGER IMODE
      INTEGER INLEP
      REAL RETPMU(3),RETPPI(3),RETPNU(3)

      REAL SINTH,COSTH,SINPHI,COSPHI
      INTEGER*4 IDUMMY

C     How many times we try increasing the ceiling in the accept-reject
      INTEGER NTRIES

#include <rscons.h>
#include <neutmodel.h>
#include <necard.h>

      INTEGER IRSISED

      REAL WMIN,E,QMAX,SMAX,W2,W,Q2
C     CWRET Multiplicative safety factor for accept-reject ceiling
      REAL SAFETY
C     CWRET added this for a W_max which is smax-xmlep or WMAX
C     WMAX parameter is set in rscons.h and is 2.0 GeV (Rein-Sehgal default)
      REAL W_max
C     CWRET added these as alternative to calculated Q2 limits
      REAL ERES,ELEP,PRES,PLEP,COS1,COS2,COSLEP

C     CWRET added these for checking Q2 is within limits
      REAL QMAXFOUND,QMINFOUND

      REAL XMLEP, XMLEP2
      REAL QMIN

      REAL RNDSIG,DSIGMX
      REAL SIGMAX

      REAL TMPPMU(3),TMPPPI(3),TMPPNU(3)

      INTEGER I

      REAL RLU
      EXTERNAL RLU

      REAL PHI,COSTHETA,PHILEP

      REAL K_0, Q_0

C     Set the lepton mass
      XMLEP = 0 ! Default case
C     CC is IMODE 1,2,3 for neutrino; 11,12,13 for anti-neutrino
      IF ((IMODE.LT.4).OR.(IMODE.GT.10).AND.(IMODE.LT.14)) THEN
        IF (ABS(INLEP).EQ.12) XMLEP = XME   ! Electron
        IF (ABS(INLEP).EQ.14) XMLEP = XMMU  ! Muon
        IF (ABS(INLEP).EQ.16) XMLEP = XMTAU ! Tau
      ENDIF

C     Set up some fundamental kinematic quantities
C - shift by 1 MeV 
      WMIN = XMN+XMPI
      SMAX = (2 * XMN * E) + XMN2
      XMLEP2 = XMLEP*XMLEP

C     CWRET:
C     We have a maximum W set by Enu (or equally SMAX) 
C     crs_newff had this W_max = MIN(SQRT(SMAX)-XMLEP, WMAX),
C     which I think this is wrong!
C     Previous NEUT had W_max = WMAX, which isn't allowed for low E_nu
      W_max = MIN(SQRT(SMAX)-XMLEP, WMAX)
C     Add a little to W limits for f77 precision
C     WMIN = WMIN
C     W_max = W_max

C     Get the maximum cross-section for this given Enu (bottom of this file)
C     Safety = Multiplicative safety factor that multiplies the maximum
C     cross-section found
C     If found cross-section is larger than SIGMAX we should make SAFETY
C     bigger
C     How many times we've tried to set the ceiling
      NTRIES = 0
C     Make the TRIES a factor of SAFETY
C     Come here if the maxcrs is exceeded
      SIGMAX = 0.
      SAFETY = 1.5
      CALL GETMAXCRS_MK(E, W_max, XMLEP, SAFETY, IMODE, SIGMAX)
c111  NTRIES = NTRIES+1
c     IF (NTRIES .gt. 1) THEN
c       write(*,*) "Increasing ceiling to ", SIGMAX*(NTRIES+1.0)/3.0
c       write(*,*) "*********************"
c       SIGMAX = SIGMAX*(NTRIES+1.0)/3.0
c     ENDIF

C     Set the minimum and maximum Q2 limits for throwing Q2 in
C     accept-reject
C     Important: Set ONLY so it's dependent on E_nu. If this is set to
C     be dependent on W it will significantly bias and correlate W and
C     Q2 in the kinematics.
C     This was a subtle but significant bug which caused a lot of pain!

      QMIN = 0
      QMAX = -4*E*E

C     Try these later! (might make our space slightly smaller :D)
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
         write(*,*) "MKDT33: TOO MANY TRY/STOP..."
         write(*,*) "E=", E 
         write(*,*) "W=", W, "Wmin=", WMIN, "Wmax=", W_max
         write(*,*) "Q2=", Q2, "QMIN=", QMIN, "QMAX=", QMAX
         write(*,*) "QMINFOUND=", QMINFOUND, "QMAXFOUND", QMAXFOUND
         write(*,*) "INLEP(PDG)=", INLEP
         write(*,*) "SIGMAX=", SIGMAX, "RNDSIG=", RNDSIG,
     & "DSIGMX = ", DSIGMX, "RLU(IRSISED) = ", RLU(IRSISED)
         write(*,*) "I live in src/neutcore/mkdt33.F, ca line 170"
         STOP
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
C     Add a small numerical precision to QMINFOUND and QMAXFOUND
C     Might be needed in edges of phase space (N.B. both are negative!)
C     QMIN = QMINFOUND
C     QMAX = QMAXFOUND

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
      ENDIF

      IF (QMIN .LE. QMAX) THEN
        WRITE(*,*) "Invalid Q2LIM: Q2min >= Q2max"
        GOTO 11
      ENDIF

C     Rethrow W, Q2 and maximum sigma if unphysical:
C     Same checks are made in finding maximum cross-section

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

      CosLep=(q2-XMLEP**2+2*E*ELEP)/(2*E*PLEP)
C     Angle between lepton and neutrino
      IF (ABS(CosLep) .GT. 1) GOTO 11

C     Angle between Resonance and Muon
      COS1 = (XMLEP2+W2+2*ELEP*ERES - XMN2 - 2*E*XMN)
     $      /(2*PLEP*PRES)
      IF (ABS(COS1) .GT. 1) GOTO 11

C     Angle between Delta and Neutrino
      COS2 = (XMLEP2 + XMN2 - 2*ELEP*XMN - W2 + 2*E*ERES)
     &  /(2*E*PRES)
      IF (ABS(COS2) .GT. 1) GOTO 11

C     Check the Q2
      IF ( ((2 * XMN * E) * (2 * XMN * E + XMN2 - W2) /
     $     (2 * XMN * E + XMN2) + Q2) .LT. -1.0E-5) GOTO 11

C     Check neutrino energy
      IF ((W2 - XMN2 - Q2) / (2 * XMN) .GT. E) GOTO 11

      Q_0 = (W2 - XMN2 + XMPI2)/(2.*W)
      if (Q_0*Q_0 - XMPI2 .LT. 0) THEN
        write(*,*) "Q_02 < mpi2"
        GOTO 11
      endif

      K_0 = (W2 - XMN2 + Q2)/(2.*W)
      if (K_0*K_0 - Q2 .LT. 0) then
        write(*,*) "K_02 < Q2"
        GOTO 11
      endif


C     Generate angle of Pion in ADLER FRAME randomly
      COSTHETA = 1.-RLU(IRSISED)*2.
      PHI = RLU(IRSISED) * 3.14159265 * 2.

C     If E, W, Q2, COSTHETA, PHI combination is kinemaitcally invalid, MKDCRS
C     returns 0, so DSIGMX = 0 below.
C     Then in the accept-reject, such an invalid configuration is 
C     always rejected because RNDSIG will always be greater than 0.
C     This then causes a re-throw so should be OK
      CALL MKDCRS(IMODE, E, Q2, W, COSTHETA, PHI, XMLEP, DSIGMX)

C     CWRET:  This is the final accept-reject step.
C             It takes the longest time to run.
      IF (RNDSIG .GT. DSIGMX) GOTO 11


C     Check that DSIGMX (calculated x-sec for thrown Enu, Q2, W)
C     is less than the SIGMAX used in the accept-reject method.
C     If this is the case, need to increase SIGMAX or else method fails
C     i.e. change SIGMAX = SIGMAX * (bigger number than currently)
      IF (DSIGMX .GT. SIGMAX) THEN
        write(*,*) "Calculated x-sec > accept-reject ceiling"
        write(*,*) "E = ", E, "W = ", W, "Q2 = ", Q2,
     &             "cos(theta_Ad) = ", COSTHETA, " phi_Ad = ",phi,
     &             "x-sec = ", DSIGMX, "SIGMAX = ", SIGMAX
        stop
      ENDIF

C     enu_com = e
C     q2_com = q2
C     w_com = w
C     costh_com = costheta
C     phi_com = phi

C     write(*,*) "Calculated cross-section: ", DSIGMX
      

C     Calculate kinematics (Lorentz boost from resonant rest frame to lab frame)
      CALL RSLZBT(E,Q2,W,XMLEP,ACOS(COSTHETA),PHI,TMPPMU,TMPPPI,TMPPNU)

      SINTH=0.
      COSTH=1.
      PHILEP = 3.14159265*2.*RLU(IDUMMY)
      SINPHI = sin(PHILEP)
      COSPHI = cos(PHILEP)

      CALL RS3DRT(TMPPMU,SINTH,COSTH,SINPHI,COSPHI)
      CALL RS3DRT(TMPPPI,SINTH,COSTH,SINPHI,COSPHI)
      CALL RS3DRT(TMPPNU,SINTH,COSTH,SINPHI,COSPHI)

C     Write the outgoing particle kinematics
      DO 80 I=1,3
         RETPMU(I) = TMPPMU(I)
         RETPPI(I) = TMPPPI(I)
         RETPNU(I) = TMPPNU(I)
 80   CONTINUE

      END


***********************************************************************
*     ---------------------------------------------------
      SUBROUTINE GETMAXCRS_MK(E, W_max, XMLEP, SAFETY, IMODE, SIGMAX)
*     ---------------------------------------------------
*
*     ( purpose )
*       To calculate and return the maximum cross-section in W, Q2 theta
*       and phi for given Enu to establish the accept-reject ceiling
*       For Minoo Kabirnezhad model
*
*
*     ( input )
*       SMAX        : Maximum S for interaction
*       WMAX        : Maximum W for interaction
*       SAFETY     : Multiplicative 'safety factor' for accept-reject ceiling
*
*     ( output )
*       MAXCRS      : Maximum cross-section * ceiling
*
*     ( creation date and author )
*       2016.Oct  ; Clarence Wret and Minoo Kabirnezhad
*
***********************************************************************

      IMPLICIT NONE

      REAL E, W_MAX, SAFETY
      INTEGER IMODE

#include <rscons.h>
      REAL XMLEP
      REAL SMAX

C     Related to scanning W
      REAL WMINit, WMAXit, dW, W, Wpeak
C     Related to scanning Q2
      REAL Q2MIN, Q2MAX, dQ2, Q2, Q2peak
C     Related to scanning costheta
      REAL COSTHMIN, COSTHMAX, dCOSTH, COSTH, COSTHpeak
C     Related to scanning phi
      REAL PHImin, PHImax, dPHI, PHI, PHIpeak

      INTEGER INW, NW, IQ2, NQ2, NPHI, ICOSTH, IPHI, NCOSTH

C     x-sec related variables
      REAL SIGMAX
      real diff_xsec

      real Q2minKin, Q2maxKin

      SMAX = (2 * XMN * E) + XMN**2
      W_max = MIN(SQRT(SMAX)-XMLEP, WMAX)
      SIGMAX = 0.

C     Step around in W space a bit to find maximum (1.18 - 1.25 here)
C     Peak in W is very likely to be around Delta(1232)
C     Can be around secondary resonances too, but should be comparable
C     in cross-section to Delta(1232), and our safety factor will save
C     us in that case
      WMINit = 1.17
      WMAXit = 1.23

C     Might run into phase space where 1.18 is not allowed
C     And don't place right on the boundary
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

      DO 5 INW = 0, NW
        W = WMINit + dW*INW

C       Leave these hard-coded for now
C       We're only interested in the peak, so look in low Q2 region
        Q2MIN = -0.02
        Q2MAX = -0.10

C       Get the physically allowed Q2 limits and check that the above
C       hard-coded limits are allowed. If not, set QMIN and QMAX
        CALL RES_Q2_LIMITS(SMAX, W, XMLEP, Q2MINkin, Q2MAXkin)

        IF (Q2MINkin .EQ. -9999 .AND. Q2MAXkin .EQ. 9999) THEN
          GOTO 5
        END IF

        if (q2minkin .ne. q2maxkin) then 
C       Don't start to scan at the minimum Q2
          IF (Q2MIN .GT. Q2MINkin) THEN
            Q2MIN = Q2MINkin-1E-4
          END IF
C       And don't place it right at the boundary, push it a few percent
C       to the side
          IF (Q2MAX .LT. Q2MAXkin) THEN
            Q2MAX = Q2MAXkin+1E-4
          END IF
        end if

C       number of Q2 slices to consider
        NQ2   = 8

C       Define DQ2 to be negative
        DQ2   = (Q2MAX-Q2MIN)/NQ2

        DO 50 IQ2 = 0,NQ2
        Q2 = Q2MIN + DQ2*IQ2

C         Scan cos theta
        COSTHmin = -0.99
        COSTHmax = 0.99
        NCOSTH = 5
        DCOSTH = (COSTHmax - COSTHmin)/NCOSTH
        DO 500 ICOSTH = 0, NCOSTH
        COSTH = COSTHmin + ICOSTH*dCOSTH

C           And finally scan phi
        PHImin = 0.01
        PHImax = 2*PI-.01
        NPHI = 5
        dPHI = (PHImax - PHImin)/NPHI
        DO 5000 IPHI = 0, NPHI
        PHI = PHImin + IPHI*dPHI

C             Now call the cross-section calculation and save x-sec
        call mkdcrs(imode, e, q2, w, costh, phi, xmlep, diff_xsec)

C             If cross-section at given Enu, Q2, W, COSTH, PHI
C             is less than previous iteration we've found the maximum
        IF (diff_xsec .GT. SIGMAX) THEN
C             Only update SIGMAX when TMPSIG > previous TMPSIG (=SIGMAX)
          SIGMAX = diff_xsec
          Wpeak = W
          Q2peak = Q2
          COSTHpeak = COSTH
          PHIpeak = PHI
        ENDIF

C         Finished looping PHI
 5000       CONTINUE

C       Finished looping PHI
 500      CONTINUE

C     Finished looping Q2
 50    CONTINUE

C     Finished looping W, now move on in life
 5    CONTINUE

C     Finally multiply found maximum by safety factor
        SIGMAX = SIGMAX * SAFETY

C     Unfortunately this hack is needed for the case when the Q2 limits
C     are so close to each other that a cross-section isn't calculated
        IF (SIGMAX .LT. 1E-7) SIGMAX = 1E-7

        END
