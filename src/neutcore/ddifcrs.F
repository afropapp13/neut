C     ------------------------------------------------------------------

      REAL*8 FUNCTION DDIFCRS(ENU,CO,ELEP,ITYPE)
C
C     CALCULATE DIFRENTIAL CROS SECTION
C     DSIGMA/DCO/DELEP
C     nu + Nuclear -> lepton + X
C
C     INPUT
C     ENU  : NEUTRINO ENERGY
C     CO   : COSINE OF SCATTERING ANGLE
C     ELEP : LEPTON ENERGY
C     ITYPE: e-(anti)neutrino  -> +(-)12, mu-(anti)neutrino -> +(-)14
C
C     OUTPUT
C     DDIFCRS : DIFFERENTIAL CROSS SECTION [ fm^2/MeV]
C
C

      IMPLICIT  NONE

C
#include "necard.h"
#include "neutparams.h"
#include "neutmodel.h"

C
C     ARGS
C
      REAL*4 ENU,CO,ELEP
      INTEGER ITYPE
C
C     PARAMETERS
C
      REAL*8 ALPHA, PI, HC,GF,XMN,DMN
      REAL*8 B,COSB2,PF,EBI,EBF
      INTEGER NUC
C
C     PARAMETERS
C
      PARAMETER  (
     &     ALPHA=1.D0/137.036D0, ! fine structure const
     &     PI=3.1415926D0,      ! Pi
     &     HC=197.326D0,        ! hbar * c [MEV FM]
     &     GF=1.1664D-11,       ! weak coupling IN MEV**-2
     &     B=1.76D0,            ! ??? FM
     &     COSB2=.9494D0,       ! cosine of cabibo angle
     &     XMN=939.6D0,         ! nucleon mass [MeV]
     &     DMN = 0.D0           ! dfference between init. & final
C     &     NUC = 8              ! number of proton or neutron
     &     )

C     &   PF = 225.D0,            ! Fermi momentum
C     &   EBI = 27.D0,            ! binding energy of initial state
C     &   EBF = 0.D0,             ! binding energy of final state

C
C
C     LOCAL VARS
C
      REAL*8 W,Q,Q2,PLEP,WEF,QVEC2,AP,BP,Q2EF,XM
      REAL*8 FGE,FGM,F1,F2,FA,FP
      REAL*8 T1,T2,TA,TB,T8
      REAL*8 W1,W2,WA,WB,W8
      REAL*8 A1,A2,A3,A4,A5,A6,A7
      REAL*8 B0,B1,B2
      REAL*8 C,D,X0,S1,S2,EU,EL
      REAL*8 DBLQ2
      REAL*8 COX, DSIG
      INTEGER NUSIG

      REAL*8 DAVMASS

#ifdef GEGMCORR
      real*8 DQ1,DQ2,DQ3,DQ4,DQ5,DQ6
      real*8 GEP,GEN,GMP,GMN
      real*8 XTAU
#endif

      PF=dble(PFMAX)*1.d3
      EBI=dble(VNUINI)*(-1.d3)
      EBF=dble(VNUFIN)*(-1.d3)
      DAVMASS = dble(XMAQE)

C     ! ITYPE:14=MUON, 12=ELECTRON; NUSIG: +1=NU, -1=NUBAR

      NUSIG = ISIGN(1,ITYPE)
      IF (IABS(ITYPE) .EQ. 12) XM=0.511D0
      IF (IABS(ITYPE) .EQ. 14) XM=105.7D0
      IF (IABS(ITYPE) .EQ. 16) XM=1784.1

C     ! NUC=Z OR N OF THE TARGET

      NUC = NUMATOM/2

      DDIFCRS = 0.D0
      IF(ELEP.LE.XM) RETURN
C
C
C     KINEMATIC VARS
C

      W=ENU-ELEP
      if (w.le.0.D0) return
      PLEP=DSQRT(ELEP**2-XM**2)
      EU=DSQRT(PF**2+XMN**2)
      WEF=W+EBF-EBI
      AP=EBI*(1.D0+EBF/W)
      BP=EBF*(1.D0-EBI/W)
C      Q  = DSQRT(ENU**2+PLEP**2-2*ENU*PLEP*CO)
      DBLQ2=(ENU**2+PLEP**2-2*ENU*PLEP*CO)
      if ((DBLQ2.lt.0.).and.(DBLQ2.ge.-1.e-7)) THEN
         DBLQ2=0.
      endif
      Q  = DSQRT(DBLQ2)
      QVEC2=Q**2
      Q2 = Q**2 - W**2
      Q2EF=Q**2-WEF**2  - DMN
C
C
C     FORM FACTORS
C


      FGE = (1.D0+Q2/7.1D5)**(-2)
      FGM = (1.D0+3.71D0)*FGE

#ifdef GEGMCORR

      if (q2.le.0.) q2=0.
      DQ2=Q2/(1.D6)
      DQ1=DSQRT(DQ2)
      DQ3=DQ2*DQ1
      DQ4=DQ2**2
      DQ5=DQ1*DQ4
      DQ6=DQ2**3

      GEP=1./(1-0.21867*DQ1 + 5.89885*DQ2
     $         -9.96209*DQ3 +16.23405*DQ4
     $         -9.63712*DQ5 + 2.90093*DQ6)

      XTAU=DQ2/(4.D0*(XMN/1.D3)**2)

      GEN=-1.*(-1.91304270D0)*(0.942D0*XTAU)/(1.D0+4.61D0*XTAU)
     $     *(1./(1.+DQ2/0.71)**2)

      GMP= 2.79285/(1-0.43584*DQ1 + 6.18608*DQ2
     $               -6.25097*DQ3 + 6.52819*DQ4
     $               -1.75359*DQ5 + 0.28736*DQ6)

      GMN= -1.913 /(1-0.40468*DQ1  + 5.6569*DQ2
     $               -4.6645 *DQ3  + 3.8811*DQ4)


      FGE = GEP-GEN
      FGM = GMP-GMN
#endif

      F1 = (FGE+Q2*FGM/XMN**2/4.D0)/(1.D0+Q2/XMN**2/4.D0)
      F2 = .5D0*(FGE-FGM)/XMN /(1.D0+Q2/XMN**2/4.D0)

C     FA=-1.262*(1+Q2/1032**2)**(-2)
C     FA=-1.232D0*(1+Q2/1010.D0**2)**(-2)
      FA=-1.232D0*(1+Q2/((DAVMASS*1.D3)**2))**(-2)
C     PV***
C     FP=2*XMN/XRED*FA/(Q2+139.57**2)
      FP=2*XMN*FA/(Q2+139.57D0**2)

C
C
C     T'S
C
      T1=.5D0*Q2*(F1- 2.D0*XMN*F2)**2+(2.D0*XMN**2+.5D0*Q2)*FA**2
      T2=2.D0*XMN**2*(F1**2+Q2*F2**2+FA**2)
      TA=XMN**2*(2.D0*XMN*F1*F2+ (.5D0*Q2-2.D0*XMN**2)*F2**2
     >     -2.D0*XMN*FA*FP + .5D0*Q2*FP**2)
      TB=-.5D0*T2
      T8=2.D0*XMN**2*FA*(F1- 2.D0*XMN*F2)

C
C
C     B'S etc
C

      if (q.eq.0.d0) return

      C=-WEF/Q
      D=Q2EF/(2.D0*Q*XMN)
      X0=3.D0*NUC/(4.D0*Q*PF**3)
      IF(WEF.GT.0) THEN
C         DBLS=(1.D0-C**2+D**2)
C         if ((DBLS).le.0.D0) DBLS=0.D0
         S1=XMN*(C*D+DSQRT(1.D0-C**2+D**2))/(1.D0-C**2)
C         S1=XMN*(C*D+DSQRT(DBLS))/(1.D0-C**2)
         S2=EU-WEF
         EL=MAX(S1,S2)
      ENDIF

      IF(WEF.LT.0.D0 .OR.  EU.LT.EL) RETURN

      B0=X0*(EU- EL+AP*DLOG((EU-EBI)/(EL-EBI))
     >     +BP*DLOG((EU-EBI+W)/(EL-EBI+W)))
      B1=X0/XMN*(.5D0*(EU**2-EL**2)+
     >     AP*(EU-EL+EBI*DLOG((EU- EBI)/(EL-EBI)))
     >     +BP*((EU-EL)+(EBI-W)*DLOG((EU-EBI+W)/(
     >     EL-EBI+W))))
      B2=X0/(XMN**2)*(1.D0/3.D0*(EU**3-EL**3)
     >     +AP*(.5D0*(EU**2-EL**2)+EBI*(EU-EL)
     >     +EBI**2*DLOG((EU-EBI)/(EL-EBI)))
     >     +BP*(.5D0*(EU**2-EL**2)+(EBI-W)*
     >     (EU-EL)+(EBI-W)**2*DLOG((
     >     EU-EBI+W)/(EL-EBI+W))))

C
C
C     A'S
C
      A1=B0
      A2=B2-B0
      A3=C**2*B2+2.D0*C*D*B1+D**2*B0
      A4=B2-2.D0*EBI/XMN*B1+EBI**2/(XMN**2)*B0
      A5=C*B2+(D-EBI*C/XMN)*B1-EBI*D/XMN*B0
      A6=C*B1+D*B0
      A7=B1-EBI/XMN*B0

C
C
C     W'S
C
      W1=A1*T1+.5D0*(A2- A3)*T2
      W2=(A4+2*W/Q*A5+W**2/(Q**2)*A3+.5D0*Q2/(Q**2)
     >     *(A2-A3))*T2
      WA=1.D0/(Q**2)*(1.5D0*A3-.5D0*A2)*T2
     >     +1.D0/XMN**2*A1*TA
     >     +2.D0/(XMN*Q)*A6*TB
      WB=1.D0/XMN* (A7+W/Q*A6)*TB
     >     + W/(Q**2)*(1.5D0*A3-.5D0*A2+Q/W*A5)*T2
      W8=NUSIG/XMN*(A7+W/Q*A6)*T8

C
C
C     CROSS SECTION
C
      COX=PLEP/ELEP*CO
      DSIG=COSB2*(GF*HC)**2*ELEP*PLEP/PI *.5D0 *
     >     (W2*(1.D0+COX)+(2.D0*W1+XM**2*WA)*(1.D0-COX)
     >     -2.D0*W8*(ENU+ELEP)*(1.D0-COX)
     >     +(WB+W8)*2.D0*XM**2/ELEP )
      DDIFCRS = DSIG
      RETURN
      END


