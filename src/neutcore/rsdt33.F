***********************************************************************
*     ---------------------------------------------------
      SUBROUTINE RSDT33(IMODE,INLEP,E,RETPMU,RETPPI,RETPNU)
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
*
*     ( comment )
*       Reference:Rein and Sehgal,Ann. of Phys. 133(1981),79-153
*                 Rein           ,Z.Phys.C 35(1987),43-64
*                 Feynman et al. ,Phys.Rev.D3(1971),2706-
*                 K.S.Kuzmin et al., hep-ph/0312107
*                 Berger and Sehgal, hep-ph/0709.4378
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

      INTEGER  INFLPCK

      INTEGER IMODE,INLEP
      REAL RETPMU(3),RETPPI(3),RETPNU(3)

      REAL SINTH,COSTH,SINPHI,COSPHI
      INTEGER*4 IDUMMY

#include<rscons.h>
#include <neutmodel.h>

      REAL*4 XMASPI2, XMVSPI2

      INTEGER IRSISED
C      COMMON /RSRNDM/IRSISED

      REAL  AM3(50),AM1(50),AP1(50),AP3(50),A0P(50),A0M(50)
      REAL AP1C(50),AP3C(50),AM1C(50),AM3C(50),A0PC(50),A0MC(50)
      REAL FFBW(50)

      INTEGER IPP(31),NRR(31)

      INTEGER IFLAG,IB,NR,IP,IBLOCK,IQ2,RES33
      REAL WMIN,E,QMAX,SMAX,W2,W,X,Q2,ES,DQ2
      REAL QV,QV2,QVS,QVS2,A,V,U,V2,U2,TWUV,AK,ALFA,ALFA2
      REAL BETA,BETA2,GAMA,GAMA2,FAKT
C      REAL BR,FBWNO,XE,FFV,FFA,ETA
      REAL BR,FBWNO,XE,FFV,FFA
      REAL R,RA,RV,RM,RP,T,TA,TV,TP,TM,S,B,C
      REAL FM3C,FM1C,FP3C,FP1C,F0PC,F0MC
      REAL FM3,FM1,FP3,FP1,F0P,F0M
      REAL XXE,RAD,RADR,PQ,PQR,FBW,BRA

      REAL XNU,XNU2,XNUS,XLAM,XLAM2,XMR
      REAL XMLEP,QMIN,ELEPL,APLEPL,XMLEP2

C
      REAL PNEUCM,XCONS,PTAUA,PTAUA2

      REAL GV3,GV1,NFF,CV1,CV2,CV3A,CV3B,CV4,CV5
      REAL CA5

      COMPLEX ETABW(31)
      REAL ROP3P3,ROM3M3,ROP3P1,ROM1M3,ROP3M1,ROP1M3,ROP1P1,ROM1M1
      REAL RHOTIL,ROTI33,ROTI31,ROT3M1

      REAL RESSGN(31)
      DATA RESSGN/1.,-1.,-1.,-1.,-1.,1.,1.,1.,-1.,-1.,1.,1.,1.,1.,1.,
     $     1.,1.,1.,1.,1.,1.,1.,-1.,-1.,1.,1.,1.,1.,1.,1.,1. /

      REAL PAR(14)
      REAL RETF(6)

      REAL QTRPI,YRAND,MAXDST,YFN

      REAL RNDSIG,DSIGMX,DSIGMXP,DSIGMXM,DSIG3,DSIG3P,DSIG3M
      REAL TMPSIG,TMPSIGP,TMPSIGM,SIGMAX,JUNK,UU
      INTEGER Lambda(2), iOrg

      REAL TMPPMU(3),TMPPPI(3),TMPPNU(3)

      INTEGER I

C     function type

      REAL RLU
      EXTERNAL RLU
      REAL RSY00,RSY20,RSY21,RSY22
      EXTERNAL RSY00,RSY20,RSY21,RSY22

      REAL PHI,THETA,PHILEP
      REAL RHOP3P3,RHOP3P1,RHOP3M1


C

      DATA IPP/3,5,5,1,1,1,1,1,5,5,5,5,5,3,3,3,7,7,3,3,3,3,3,7,7,3,7,7,
     1     1,1,3/
      DATA NRR/0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
     1     0,0,0/


C
C     set CONSTANTS
C
      XMASPI2= XMASPI * XMASPI
      XMVSPI2= XMVSPI * XMVSPI

      XMLEP = 0.00000
      if ((IMODE.lt.4).or.((IMODE.gt.10).and.(IMODE.lt.14)))
     $     then
         IF (ABS(INLEP).EQ.12) XMLEP=XME
         IF (ABS(INLEP).EQ.14) XMLEP=XMMU
         IF (ABS(INLEP).EQ.16) XMLEP=XMTAU
      endif

      WMIN = XMN+XMPI

      PNEUCM=sqrt(E*E*XMN2/(XMN2+2*E*XMN))

      SMAX = (2 * XMN * E) + XMN2


 3    CONTINUE

C
C     set iflag
C
      IFLAG = IMODE
      IF (IMODE.GT.10) IFLAG = IFLAG -10


C     Lepton mass effects are considered
C
*     iOrg = 0 is considering lepton mass effects
*     iOrg = 1 is original R-S model(eg. NC, HighEnu)
      iOrg = 0                  ! CC
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

C
C     get maxmimum of differential cross-section
C

      SIGMAX=0.

      W=1.232
      XCONS =SMAX-XMLEP**2-W**2

      PTAUA2=( (XCONS**2-4*(XMLEP**2)*(W**2))
     $        /(4*XCONS+4*W**2+4*XMLEP**2))


      IF (PTAUA2.LT.0) THEN
         SIGMAX = 1000.
C         QMAX = (SMAX - XMN2)
C         QMIN = (XMLEP**2-2*E**2+2*E*SQRT(E**2-XMLEP**2))
C         DQ2 = -1.*QMAX/100.
         GOTO 10
      ENDIF

      PTAUA =  sqrt(PTAUA2)
      PNEUCM=  sqrt(E*E*XMN2/(XMN2+2*E*XMN))
      QMIN  =  XMLEP**2
     $       - 2*PNEUCM*sqrt(PTAUA**2+XMLEP**2)
     $       + 2*PNEUCM*PTAUA
      QMAX  =  XMLEP**2
     $       - 2*PNEUCM*sqrt(PTAUA**2+XMLEP**2)
     $       - 2*PTAUA*PNEUCM

      DQ2   =  (QMAX-QMIN)/100.
      Q2    =  QMIN


      DO 5 IQ2=0,100
         Q2 = Q2 + DQ2

C     Lepton currents are no longer conserved
C      due to the helicity of final state lepton
         CALL RSDCRS(IMODE,iOrg,Lambda(1),XMLEP,E,Q2,1.232,TMPSIGP,JUNK)
         CALL RSDCRS(IMODE,iOrg,Lambda(2),XMLEP,E,Q2,1.232,TMPSIGM,JUNK)
         TMPSIG = TMPSIGP+TMPSIGM

         IF (TMPSIG.LT.SIGMAX) THEN
            GOTO 10
         ENDIF
         SIGMAX=TMPSIG
C         write(*,*) SIGMAX
 5    CONTINUE

 10   CONTINUE

      SIGMAX = SIGMAX * 10.0

C
C     fix W,q2
C
C-FOR CHECK INFINITE LOOP
      INFLPCK=0
C-

 11   CONTINUE

      INFLPCK=INFLPCK+1
      IF (INFLPCK.gt.100000000) THEN
         CALL RSTCRS(IMODE,E,SIGMAX,TMPSIG)
         write(*,*) "RSDT33: TOO MANY TRY/STOP..."
         write(*,*) "E=",E,"/IMODE=",IMODE,"/INLEP=",INLEP,
     $              "CRSSECT=",SIGMAX
         STOP
      ENDIF

      RNDSIG = RLU(IRSISED)*SIGMAX
      W  = RLU(IRSISED)*(WMAX - WMIN)+WMIN
C      if (e.le.1.) then
C         Q2 = -1.*RLU(IRSISED)*8
C      else
         Q2 = -1.*RLU(IRSISED)*e*e*8
C      endif

 12   CONTINUE
      XCONS=SMAX-XMLEP**2-W**2
      PTAUA2=( (XCONS**2-4*(XMLEP**2)*(W**2))
     $     /(4*XCONS+4*W**2+4*XMLEP**2))

      IF (PTAUA2.lt.0) GO TO 11
      PTAUA = sqrt(PTAUA2)

C      QMIN  =  XMLEP**2
C     $      - 2*PNEUCM*sqrt(PTAUA**2+XMLEP**2)
C     $      + 2*PNEUCM*PTAUA
C      QMAX  =  XMLEP**2
C     $      - 2*PNEUCM*sqrt(PTAUA**2+XMLEP**2)
C     $      - 2*PTAUA*PNEUCM

C      Q2 = QMIN + RLU(IRSISED)*(QMAX-QMIN)

      W2 = W*W
      IF ( ((2 * XMN * E) * (2 * XMN * E + XMN2 - W2) /
     $     (2 * XMN * E + XMN2)).LT.(-Q2) ) GOTO 11

      IF (((W2+XMN2-Q2)/(2*XMN)).LE.W) GOTO 11
	  IF ((Q2-W2+XMN2+2*XMN*E)/(2*XMN).lt.XMLEP) GOTO 11

      XMLEP2 = XMLEP*XMLEP
      ELEPL=(2*E*XMN+XMN2+Q2-W2)/(2*XMN)
      APLEPL=ELEPL**2-XMLEP2
      if (APLEPL.LE.0) goto 11
      APLEPL=SQRT(APLEPL)

C-- add to avoid the case APLEPL is almost 0 / Sep.17,2001 by Y.Hayato
      if (elepl.lt.1.e-7) goto 11

      IF (abs(Q2+2*E*ELEPL-XMLEP2)/(2*E*APLEPL).GT.1.)
     $     GOTO 11

C     Lepton currents are no longer conserved
C      due to the helicity of final state lepton
      CALL RSDCRS(IMODE,iOrg,Lambda(1),XMLEP,E,Q2,W,DSIGMXP,DSIG3P)
      CALL RSDCRS(IMODE,iOrg,Lambda(2),XMLEP,E,Q2,W,DSIGMXM,DSIG3M)

      DSIGMX = DSIGMXP+DSIGMXM
      DSIG3  = DSIG3P +DSIG3M

      IF (DSIGMX.LT.RNDSIG) GOTO 11

C      write(92,*) q2

C
C     select mode(33-resonance or not.)
C
      RES33 = 0
      TMPSIG = RLU(IRSISED) * DSIGMX
      IF (TMPSIG.LT.DSIG3) THEN
         RES33 = 1
      ENDIF

C 1    call HF1(30,q2,1.)

C
C     main
C

      ROP3P3 = 0.
      ROM3M3 = 0.
      ROP3P1 = 0.
      ROM1M3 = 0.
      ROP1M3 = 0.
      ROP3M1 = 0.
      ROP1P1 = 0.
      ROM1M1 = 0.


C     -------------------------------------------------------
C     NOW  W  AND  E  ARE FIXED
C     Q2MAX = 2.*MN*E*(2.*MN*E+MN2-W2)/(2.*MN*E+MN2)
C     -------------------------------------------------------

C     KINEMATICAL QUANTITIES

      XNU = (W2 - XMN2 - Q2) / (2 * XMN)
      X = -Q2 / (2 * XMN * XNU)

      IF (XNU.GE.E) GOTO 11

      XNU2 = XNU**2
      ES = E-XNU

      IF (XNU2.LE.Q2) THEN
         WRITE(*,*) "ERROR"
         STOP
      ENDIF

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
      IF(V.LT.0.) V = 0.
      U2 = U**2
      V2 = V**2
      TWUV = 2.*U*V
      AK = (W2-XMN2)/(2.*XMN)
      ALFA = SQRT(-Q2/QV2)*U
      BETA = SQRT(-Q2/QV2)*V
      GAMA = (XMN/W)*SQRT(2.*U*V)
      ALFA2 = ALFA**2
      BETA2 = BETA**2
      GAMA2 = GAMA**2
      FAKT = W*PI/(W2-XMN2)

C     SUMMATION OVER RESONANCES

      IB = 1

      IBLOCK = IB-1
      NR = NRR(IB)
      XMR = XMRR(IB)
      BR = BRR(IB)

C        SETS UP NORMALIZATION CONSTANTS
C        G&S don't use the BW normalizing factors for their model 
C        NOTE: This should only work for G&S form factors only
      if(IFFSPI.ne.0) then
        FBWNO = 1
      else
        FBWNO = FBWNOO(IB)
      endif

C      ETA = ETAA(IB)
      XE = XEE(IB)
      IP = IPP(IB)

C     **********************************

C     FORMFACTORS AND RELATED QUANTITIES

      if (IFFSPI.ne.0) then

C     K.M. Graczyk and J.T. Sobczyk, Phys.Rev.D 77,053001 (2008)
C     Equations are correction in erratum
C     CA5(0) and XMA2 are set by the user

        CA5=   RCA5ISPI/((1-Q2/XMASPI2)**2)

C       J.A. Nowak, arXiv:0909.3659
C       Vector Form Factors
C       CV1, CV2 are for simplicity, not true form factors
C       Added the proportional factors of the vector form factors GV
C       to CV3,CV4,CV5

C       (PAR) Yargh, the variable Q2 here has a minus sign relative to
C       Q2 in the Nowak/Graczyk & Sobczyk papers

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

        FFV = 0.5 * (1 - Q2/(XMN+W)**2)**0.5
     &        * (1 - Q2/(4*XMN2))**(-1*NR)
     &        * NFF

C       Axial. Divide by Z, the arbitrary Rein-Sehgal factor, because
C       we multiply by it later

        FFA = 0.5 * SQRT(3.0) * (1 - Q2/(XMN+W)**2)**0.5
     &        * (1 - Q2/(4*XMN2))**(-1*NR)
     &        * (1 - (W2 + Q2 - XMN2)/(8*XMN2))
     &        * CA5 / Z

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

      CALL RSCLFM(IFLAG,IBLOCK,PAR,RETF)

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

      XXE = SQRT(XE)

      RAD =((W2-XMN2-XMPI2)**2-4.*XMN2*XMPI2)/(4.*W2)
      IF(RAD.LE.0.) RAD = 0.0001
      PQ = SQRT(RAD)
      RADR = ((XMR**2-XMN2-XMPI2)**2-4.*XMN2*XMPI2)/(4.*XMR**2)
      IF(RADR.LE.0.) RADR = 0.0001
      PQR = SQRT(RADR)
      BRA = BR*(PQ/PQR)**IP
      FBW = (BRA/(2.*PI))/((W-XMR)**2+(BRA**2)/4.)

      ETABW(IB) = CMPLX(W-XMR,-1.*BRA/2.)
      ETABW(IB) = ETABW(IB) / ((W-XMR)**2+(BRA**2)/4.)
      ETABW(IB) = ETABW(IB)*SQRT(BRA / (2.*PI) / FBWNO)

      FFBW(IB) = FBW/FBWNO

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
 991  CONTINUE
      AM3C(IB) = FM3C*XXE
      AM1C(IB) = FM1C*XXE
      AP1C(IB) = FP1C*XXE
      AP3C(IB) = FP3C*XXE
      A0PC(IB) = F0PC*XXE
      A0MC(IB) = F0MC*XXE
 992  CONTINUE

C     Neutrino or anti neutrino ??

      IF (IMODE.GT.10) THEN
         UU= U
         U = -V
         V = -UU
ccc this was bug 98 MAR 11 by Y.Itow
ccc         TWUV = -TWUV
      ENDIF


      ROP3P3 = ROP3P3
     $     + ( W * E * SQRT(-Q2 / QVS) * V * AM3C(1))
     $     *( W * E * SQRT(-Q2 / QVS) * V * AM3C(1))
     $     *(ABS(ETABW(1))**2)
      ROM3M3 = ROM3M3
     $     + (-W * E * SQRT(-Q2 / QVS) * U * AP3C(1))
     $     *(-W * E * SQRT(-Q2 / QVS) * U * AP3C(1))
     $     *(ABS(ETABW(1))**2)
      ROP3P1 = ROP3P1
     $     + ( W * E * SQRT(-Q2 / QVS) * V * AM3C(1))
     $     *( W * E * XMN / W * SQRT(TWUV) *A0PC(1))
     $     *(ABS(ETABW(1))**2)
      ROM1M3 = ROM1M3
     $     + ( W * E * XMN / W * SQRT(TWUV) *A0MC(1))
     $     *(-W * E * SQRT(-Q2 / QVS) * U * AP3C(1))
     $     *(ABS(ETABW(1))**2)
      ROP3M1 = ROP3M1
     $     + ( W * E * SQRT(-Q2 / QVS) * V * AM3C(1))
     $     *(-W * E * SQRT(-Q2 / QVS) * U * AP1C(1))
     $     *(ABS(ETABW(1))**2)
      ROP1M3 = ROP1M3
     $     + ( W * E * SQRT(-Q2 / QVS) * V * AM1C(1))
     $     *(-W * E * SQRT(-Q2 / QVS) * U * AP3C(1))
     $     *(ABS(ETABW(1))**2)

      ROP1P1 = ROP1P1
     $     + ( W * E * XMN / W * SQRT(TWUV) *A0PC(1))
     $     *( W * E * XMN / W * SQRT(TWUV) *A0PC(1))
     $     *(ABS(ETABW(1))**2)
     $     + ( W * E * SQRT(-Q2 / QVS) * V * AM1C(1))
     $     *( W * E * SQRT(-Q2 / QVS) * V * AM1C(1))
     $     *(ABS(ETABW(1))**2)

      ROM1M1 = ROM1M1
     $     + ( W * E * XMN / W * SQRT(TWUV) *A0MC(1))
     $     *( W * E * XMN / W * SQRT(TWUV) *A0MC(1))
     $     *(ABS(ETABW(1))**2)
     $     + (-W * E * SQRT(-Q2 / QVS) * U * AP1C(1))
     $     *(-W * E * SQRT(-Q2 / QVS) * U * AP1C(1))
     $     *(ABS(ETABW(1))**2)

C     restore value

      IF (IMODE.GT.10) THEN
         UU= U
         U = -V
         V = -UU
c this was bug/ fixed by Y.Itow 98 Mar 11
c         TWUV = -TWUV
      ENDIF

      RHOTIL = ROP1P1+ROM1M1+ROP3P3+ROM3M3
      ROTI33 = ROP3P3+ROM3M3
      ROTI31 = ROP3P1-ROM1M3
      ROT3M1 = ROP3M1+ROP1M3

      RHOP3P3 = ROTI33 / RHOTIL
      RHOP3P1 = ROTI31 / RHOTIL
      RHOP3M1 = ROT3M1 / RHOTIL

      QTRPI = SQRT(1 / (4 * 3.14))
      MAXDST=QTRPI + QTRPI * (2 * ABS(RHOP3P3) + 1.)
     $     + 4 * QTRPI * ABS(RHOP3P1)
     $     + 4 * QTRPI * ABS(RHOP3M1)

 50   CONTINUE

C     generate angle of Pion in ADLER FRAME randomly

      THETA = ACOS(1.-RLU(IRSISED)*2.)
      PHI = RLU(IRSISED) * 3.14 * 2
      YRAND = RLU(IRSISED) * MAXDST

      IF (RES33.EQ.0) GOTO 72

      YFN = RSY00(THETA,PHI)
     $     - ((2 / SQRT( 5.)) * RSY20(THETA,PHI) *
     $        (RHOP3P3 - 0.5))
     $     + ((4 / SQRT(10.)) * RSY21(THETA,PHI) * RHOP3P1)
     $     - ((4 / SQRT(10.)) * RSY22(THETA,PHI) * RHOP3M1)

      IF (YRAND.GT.YFN) GOTO 50

 72   CONTINUE

C      write(210,*) theta * 180 / 3.14,cos(theta)
C      call HF1(10,cos(theta),1.)

      CALL RSLZBT(E,Q2,W,XMLEP,THETA,PHI,TMPPMU,TMPPPI,TMPPNU)

C      write(211,*) ((PmuLx*PpiLx+PmuLy*PpiLy+PmuLz*PpiLz)/
C     $     (PmuL*PpiL))

C- Add 98/03/01
C- Debug 98/12/19

      SINTH=0.
      COSTH=1.
      PHILEP = 3.1415*2*RLU(IDUMMY)
      SINPHI = sin(PHILEP)
      COSPHI = cos(PHILEP)

      CALL RS3DRT(TMPPMU,SINTH,COSTH,SINPHI,COSPHI)
      CALL RS3DRT(TMPPPI,SINTH,COSTH,SINPHI,COSPHI)
      CALL RS3DRT(TMPPNU,SINTH,COSTH,SINPHI,COSPHI)

      DO 80 I=1,3
         RETPMU(I) = TMPPMU(I)
         RETPPI(I) = TMPPPI(I)
         RETPNU(I) = TMPPNU(I)
 80   CONTINUE

      END
