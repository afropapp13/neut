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
      REAL*8 ALPHA,PI,HC,GF,XMN,XMPI,DMN
      REAL*8 B,COSB2,TRCORA,TRCORB
      REAL*8 PF,EBI,EBF,DAMASS,DVMASS
C
      PARAMETER  (
     >     ALPHA=1.d0/137.036d0, ! fine structure const
     >     PI=3.1415926d0,      ! Pi
     >     HC=197.326d0,        ! hbar * c [MEV FM]
     >     GF=1.1664d-11,       ! weak coupling IN MEV**-2
     >     XMN=939.6d0,         ! nucleon mass [MeV]
     >     XMPI=139.57d0,       ! pion mass [MeV]
     >     DMN=0.d0,            ! dfference between init. & final
     >     B=1.76d0,            ! ??? FM
     >     COSB2=.9494d0,       ! cosine of cabibo angle
C    >     NUC=8,               ! number of proton or neutron
C    >     PF = 225.d0,         ! Fermi momentum
C    >     EBI = 27.d0,         ! binding energy of initial state
C    >     EBF = 0.d0,          ! binding energy of final state
     >     TRCORA=6.d0,         ! ??? something to do with MEC
     >     TRCORB=.34d0         ! ??? something to do with MEC
     >     )

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
      REAL*8 COX,DSIG
      INTEGER NUC,NUSIG,FFMDL


      REAL*8 GPE,GPM,GNE,GNM
      REAL*8 TMP,BBBA07,FABBBA07
      EXTERNAL BBBA07,FABBBA07

      REAL*8 FV3,FA3

#ifdef GEGMCORR
      REAL*8 DQ1,DQ2,DQ3,DQ4,DQ5,DQ6
      REAL*8 GEP,GEN,GMP,GMN
      REAL*8 XTAU
#endif

      PF=dble(PFMAX)*1.d3           ! MeV
      EBI=dble(VNUINI)*(-1.d3)      ! MeV (>0)
      EBF=dble(VNUFIN)*(-1.d3)      ! MeV (>0)
      DAMASS = dble(XMAQE)*(1.d3)   ! MeV
      DVMASS = dble(XMVQE)*(1.d3)   ! MeV

C     ! ITYPE:14=MUON, 12=ELECTRON; NUSIG: +1=NU, -1=NUBAR

      NUSIG = ISIGN(1,ITYPE)
      IF (IABS(ITYPE) .EQ. 12) XM=0.511d0
      IF (IABS(ITYPE) .EQ. 14) XM=105.7d0
      IF (IABS(ITYPE) .EQ. 16) XM=1784.1d0

C     ! NUC=Z OR N OF THE TARGET

      NUC = NUMATOM/2

      DDIFCRS = 0.d0
      IF(ELEP.LE.XM) RETURN
C
C
C     KINEMATIC VARS
C
      W=ENU-ELEP
      IF (W.LE.0.d0) RETURN

      PLEP=DSQRT(ELEP**2-XM**2)
      EU=DSQRT(PF**2+XMN**2)
      WEF=W+EBF-EBI
      AP=EBI*(1.d0+EBF/W)
      BP=EBF*(1.d0-EBI/W)
C     Q  = DSQRT(ENU**2+PLEP**2-2*ENU*PLEP*CO)
      DBLQ2=(ENU**2+PLEP**2-2*ENU*PLEP*CO)
      IF ((DBLQ2.LT.0.).AND.(DBLQ2.GE.-1.d-7)) THEN
         DBLQ2=0.
      ENDIF
      Q  = DSQRT(DBLQ2)
      QVEC2=Q**2
      Q2 = Q**2 - W**2
      Q2EF=Q**2-WEF**2  - DMN
C
C
C     FORM FACTORS
C
      FFMDL = MOD(MDLQE,10)
      IF (FFMDL.EQ.1) THEN
C        Dipole
         FGE = (1.d0+Q2/DVMASS**2)**(-2)
         FGM = (1.d0+3.71d0)*FGE

      ELSE IF (FFMDL.EQ.2) THEN
C        BBBA05 expects Q^2 in [GeV^2]
         CALL BBBA05(Q2*1.d-6,GPE,GPM,GNE,GNM)
         FGE=GPE-GNE
         FGM=GPM-GNM

      ELSE IF (FFMDL.EQ.3) THEN
C        BBBA07 expects Q^2 in [GeV^2]
         TMP=BBBA07(Q2*1.d-6,GPE,GPM,GNE,GNM)
         FGE=GPE-GNE
         FGM=GPM-GNM

      ENDIF

      IF (MOD(MDLQE,1000)/100.EQ.1) THEN
         FGM = FGM * SQRT(1.d0+TRCORA*Q2*EXP(-Q2/TRCORB))
      ENDIF

#ifdef GEGMCORR

      IF (Q2.LE.0.) Q2=0.
      DQ2=Q2/(1.d6)
      DQ1=DSQRT(DQ2)
      DQ3=DQ2*DQ1
      DQ4=DQ2**2
      DQ5=DQ1*DQ4
      DQ6=DQ2**3

      GEP=1./(1-0.21867*DQ1 + 5.89885*DQ2
     >         -9.96209*DQ3 +16.23405*DQ4
     >         -9.63712*DQ5 + 2.90093*DQ6)

      XTAU=DQ2/(4.d0*(XMN/1.d3)**2)

      GEN=-1.*(-1.91304270d0)*(0.942d0*XTAU)/(1.d0+4.61d0*XTAU)
     >     *(1.d0/(1.d0+DQ2/0.71d0)**2)

      GMP= 2.79285d0/(1.d0-0.43584d0*DQ1 + 6.18608d0*DQ2
     >                    -6.25097d0*DQ3 + 6.52819d0*DQ4
     >                    -1.75359d0*DQ5 + 0.28736d0*DQ6)

      GMN= -1.913d0 /(1.d0-0.40468d0*DQ1  + 5.6569d0*DQ2
     >                    -4.6645d0 *DQ3  + 3.8811d0*DQ4)


      FGE = GEP-GEN
      FGM = GMP-GMN
#endif

      F1 = (FGE+Q2*FGM/XMN**2/4.d0)/(1.d0+Q2/XMN**2/4.d0)
      F2 = .5d0*(FGE-FGM)/XMN /(1.d0+Q2/XMN**2/4.d0)

      IF (MDLQEAF.EQ.2) THEN
C        BBBA07 expects Q2 in [GeV^2]
         FA=FABBBA07(Q2*1.d-6)
C     ELSE IF (MDLQEAF.EQ.1) THEN
C        Dipole model
      ELSE
C        Fall back to Dipole
         FA=-1.232d0*(1.d0+Q2/DAMASS**2)**(-2)
C        FA=-1.2673d0*(1.d0+Q2/DAMASS**2)**(-2)
      ENDIF

      FP=2.d0*XMN*FA/(Q2+XMPI**2)
#if (NECOREVER == 5321) || (NECOREVER == 5322)
C     special case for 5.3.2.1 and 5.3.2.2
C     do nothing.
#elif (NECOREVER >= 533)
C     There is now a normalization term passed by the CARD file.
      FP=FPQE*FP
#endif

C
C
C     Second class currents
C
      FV3 = -SCCFV*(1.d0+Q2/DVMASS**2)**(-2)
      FA3 = SCCFA*FA

C
C
C     T'S
C
      T1=.5d0*Q2*(F1- 2.d0*XMN*F2)**2+(2.d0*XMN**2+.5d0*Q2)*FA**2
      T2=2.d0*XMN**2*(F1**2+Q2*(F2**2+FA3**2)+FA**2)
      TA=XMN**2*(2.d0*XMN*F1*F2 +(.5d0*Q2-2.d0*XMN**2)*F2**2
     >          -2.d0*XMN*FA*FP +.5d0*Q2*FP**2
     >          +FV3*(2.d0*XMN*F1+Q2*FP) +FV3**2*(.5d0*Q2+2.d0*XMN**2)
     >          +FA3*(2.d0*XMN*FA-Q2*FP) +FA3**2*(.5d0*Q2)
     >          )
      TB=-.5d0*T2
     >   -(FV3*(2.d0*XMN*F1+Q2*F2) +FA3*(2.d0*XMN*FA-Q2*FP))*XMN**2
      T8=2.d0*XMN**2*FA*(F1- 2.d0*XMN*F2)

C
C
C     B'S etc
C

      IF (Q.EQ.0.d0) RETURN

      C=-WEF/Q
      D=Q2EF/(2.d0*Q*XMN)
      X0=3.d0*NUC/(4.d0*Q*PF**3)
      IF(WEF.GT.0) THEN
         S1=XMN*(C*D+DSQRT(1.d0-C**2+D**2))/(1.d0-C**2)
         S2=KAPP*(EU-WEF)
         EL=MAX(S1,S2)
      ENDIF

      IF(WEF.LT.0.d0 .OR. EU.LT.EL) RETURN

      B0=X0*((EU-EL)
     >      +AP*DLOG((EU-EBI  )/(EL-EBI  ))
     >      +BP*DLOG((EU-EBI+W)/(EL-EBI+W))
     >      )
      B1=X0/XMN*(.5d0*(EU**2-EL**2)
     >          +AP*((EU-EL)+(EBI  )*DLOG((EU-EBI  )/(EL-EBI  )))
     >          +BP*((EU-EL)+(EBI-W)*DLOG((EU-EBI+W)/(EL-EBI+W)))
     >          )
      B2=X0/(XMN**2)*(1.d0/3.d0*(EU**3-EL**3)
     >               +AP*(.5d0*(EU**2-EL**2)+EBI*(EU-EL)
     >                   +EBI**2*DLOG((EU-EBI)/(EL-EBI))
     >                   )
     >               +BP*(.5d0*(EU**2-EL**2)
     >                   +(EBI-W)*(EU-EL)
     >                   +(EBI-W)**2*DLOG((EU-EBI+W)/(EL-EBI+W))
     >                   )
     >               )

C
C
C     A'S
C
      A1=B0
      A2=B2-B0
      A3=C**2*B2+2.d0*C*D*B1+D**2*B0
      A4=B2-2.d0*EBI/XMN*B1+EBI**2/(XMN**2)*B0
      A5=C*B2+(D-EBI*C/XMN)*B1-EBI*D/XMN*B0
      A6=C*B1+D*B0
      A7=B1-EBI/XMN*B0

C
C
C     W'S
C
      W1=A1*T1+.5d0*(A2-A3)*T2
      W2=(A4
     >     +2*W/Q*A5
     >     +W**2/(Q**2)*A3
     >     +Q2/(Q**2)*.5d0*(A2-A3)
     >   )*T2
      WA=1.d0/(Q**2)*(1.5d0*A3-.5d0*A2)*T2
     >     +1.d0/XMN**2*A1*TA
     >     +2.d0/(XMN*Q)*A6*TB
      WB=1.d0/XMN*(A7+W/Q*A6)*TB
     >     +W/(Q**2)*(1.5d0*A3-.5d0*A2+Q/W*A5)*T2
      W8=NUSIG/XMN*(A7+W/Q*A6)*T8

C
C
C     CROSS SECTION
C
      COX=PLEP/ELEP*CO
      DSIG=COSB2*(GF*HC)**2*ELEP*PLEP/PI *.5d0 *
     >     (W2*(1.d0+COX)+(2.d0*W1+XM**2*WA)*(1.d0-COX)
     >     -2.d0*W8*(ENU+ELEP)*(1.d0-COX)
     >     +(WB+W8)*2.d0*XM**2/ELEP )
      DDIFCRS = DSIG
      RETURN
      END

