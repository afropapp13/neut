***********************************************************************
*     ---------------------------------------------------
      SUBROUTINE RSDT33E(IMODE,INLEP,E,RETPMU,RETPETA,RETPNU,IERR)
*     ---------------------------------------------------
*
*     ( purpose )
*       calculate kinematics of final state particles
*       for single eta generation
*
*     ( input )
*       imode       : interaction mode
*       inlep       : kind of incoming neutrino
*       e           : energy of incoming neutrino
*
*     ( output )
*       retpmu(3)    : momentum of outgoing lepton
*       retpeta(3)    : momentum of outgoing ETA
*       retpnu(3)    : momentum of outgoing nucleon
*       IERR         : error flag
*
*     ( creation date and author )
*       1993.Dec. ; for Super-Kamioka by Y.Hayato
*                   (Based on Rein and Sehgal's code)
*       1996.Nov. ; for Avoid Infinite Loop
*                   Check Cross-section after Fermi-momentum
*                   has been fixed.
*       1997.Aug. ; debug: q^2 minimum is larger than the allowed
*                          region(pneucm is wrong.)
*       1997.Nov. ; Made to calculate eta
*       2009.Feb. ; add lepton mass effects by G.Mitsuka
*
*
*     ( comment )
*       Reference:Rein and Sehgal,Ann. of Phys. 133(1981),79-153
*                 Rein           ,Z.Phys.C 35(1987),43-64
*                 Feynman et al. ,Phys.Rev.D3(1971),2706-
*                 K.S.Kuzmin et al., hep-ph/0312107
*                 Berger and Sehgal, hep-ph/0709.4378
*
*       ----------- IFLAG means Interaction mode ----------
*
*       IFLAG < 10 -> neutrino
*       IFLAG > 10 -> anti neutrino
*
*       IFLAG = 1:'CHARGED CURRENT, P - ETA0'
*       IFLAG = 2:'NEUTRAL CURRENT, P - ETA0'
*       IFLAG = 3:'NEUTRAL CURRENT, N - ETA0'
*
*       IFLAG =11 ... -> charge conjugate
*
***********************************************************************

      IMPLICIT NONE

      INTEGER  INFLPCK

      INTEGER IMODE,INLEP
      REAL E
      REAL RETPMU(3),RETPETA(3),RETPNU(3)

      REAL SINTH,COSTH,SINPHI,COSPHI

      INTEGER IERR
#include<rscons.h>

      INTEGER IRSISED,IDUMMY
C      COMMON /RSRNDM/IRSISED

c     Duplicate definition in COMMON /RSCONS/
c      REAL PI,OMEG,XW,EPSI

      INTEGER IPP(31),NRR(31)

      INTEGER IFLAG,IQ2
      REAL WMIN,QMAX,SMAX,W2,W,Q2,DQ2
c     E is already declared as argument.
c      REAL WMIN,E,QMAX,SMAX,W2,W,Q2,DQ2


      REAL XMLEP,QMIN,ELEPL,APLEPL,XMLEP2

C
      REAL PNEUCM,XCONS,PTAUA,PTAUA2

C

      REAL RESSGN(31)
      DATA RESSGN/1.,-1.,-1.,-1.,-1.,1.,1.,1.,-1.,-1.,1.,1.,1.,1.,1.,
     $     1.,1.,1.,1.,1.,1.,1.,-1.,-1.,1.,1.,1.,1.,1.,1.,1. /

      REAL RNDSIG,DSIGMX,DSIGMXP,DSIGMXM,DSIG3,DSIG3P,DSIG3M
      REAL TMPSIG,TMPSIGP,TMPSIGM,SIGMAX,JUNK,UU
      INTEGER Lambda(2), iOrg

      REAL TMPPMU(3),TMPPETA(3),TMPPNU(3)

      INTEGER I

C     function type

      REAL RLU
      EXTERNAL RLU
      REAL RSY00,RSY20,RSY21,RSY22
      EXTERNAL RSY00,RSY20,RSY21,RSY22

      REAL PHI,THETA,PHILEP

C

      DATA IPP/3,5,5,1,1,1,1,1,5,5,5,5,5,3,3,3,7,7,3,3,3,3,3,7,7,3,7,7,
     1     1,1,3/
      DATA NRR/0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
     1     0,0,0/


      IERR = 0

C
C     set CONSTANTS
C
      XMLEP = 0.00000
      if ((IMODE.le.1).or.((IMODE.gt.10).and.(IMODE.le.11)))
     $     then
         IF (ABS(INLEP).EQ.12) XMLEP=XME
         IF (ABS(INLEP).EQ.14) XMLEP=XMMU
         IF (ABS(INLEP).EQ.16) XMLEP=XMTAU
      endif

      WMIN = XMN+XMETA

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
      if (IFLAG.GE.2) iOrg = 1  ! NC

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

      W=1.535
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
         CALL RSDCRSE(IMODE,iOrg,Lambda(1),XMLEP,E,Q2,1.5,TMPSIGP,JUNK)
         CALL RSDCRSE(IMODE,iOrg,Lambda(2),XMLEP,E,Q2,1.5,TMPSIGM,JUNK)
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
         CALL RSTCRSET(IMODE,E,SIGMAX,TMPSIG)
         write(*,*) "RSDT33E: TOO MANY TRY/STOP..."
         write(*,*) "E=",E,"/IMODE=",IMODE,"/INLEP=",INLEP,
     $              "CRSSECT=",SIGMAX
         IERR = 1
         goto 1000
      ENDIF

C     Throw the cross-section, W and Q2
      RNDSIG = RLU(IRSISED)*SIGMAX
      W = RLU(IRSISED)*(WMAX - WMIN)+WMIN
      Q2 = -1.*RLU(IRSISED)*E*E*4

 12   CONTINUE
      XCONS=SMAX-XMLEP**2-W**2
      PTAUA2=( (XCONS**2-4*(XMLEP**2)*(W**2))
     $     /(4*XCONS+4*W**2+4*XMLEP**2))

      IF (PTAUA2.lt.0) GO TO 11
      PTAUA = sqrt(PTAUA2)
      IF (Q2 .GT. QMIN .OR. Q2 .LT. QMAX) GOTO 11

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

      IF (abs(Q2+2*E*ELEPL-XMLEP2)/(2*E*APLEPL).GT.1.)
     $     GOTO 11

ccc event occurs or not
C     Lepton currents are no longer conserved
C      due to the helicity of final state lepton
      CALL RSDCRSE(IMODE,iOrg,Lambda(1),XMLEP,E,Q2,W,DSIGMXP,DSIG3P)
      CALL RSDCRSE(IMODE,iOrg,Lambda(2),XMLEP,E,Q2,W,DSIGMXM,DSIG3M)

      DSIGMX = DSIGMXP+DSIGMXM
      DSIG3  = DSIG3P +DSIG3M

      IF (DSIGMX.LT.RNDSIG) GOTO 11


C
C     main
C

C     -------------------------------------------------------
C     NOW  W  AND  E  ARE FIXED
C     Q2MAX = 2.*MN*E*(2.*MN*E+MN2-W2)/(2.*MN*E+MN2)
C     -------------------------------------------------------


C     generate angle of Eta in ADLER FRAME randomly

      THETA = ACOS(1.-RLU(IRSISED)*2.)
      PHI = RLU(IRSISED) * 3.14 * 2

C      write(210,*) theta * 180 / 3.14,cos(theta)
C      call HF1(10,cos(theta),1.)

      CALL RSLZBTE(E,Q2,W,XMLEP,THETA,PHI,TMPPMU,TMPPETA,TMPPNU)

C      write(211,*) ((PmuLx*PpiLx+PmuLy*PpiLy+PmuLz*PpiLz)/
C     $     (PmuL*PpiL))

C- Add 98/03/01
C- Bug fixed 99/04/14

      SINTH=0.
      COSTH=1.
      PHILEP = 3.1415*2*RLU(IDUMMY)
      SINPHI = sin(PHILEP)
      COSPHI = cos(PHILEP)


      CALL RS3DRT(TMPPMU,SINTH,COSTH,SINPHI,COSPHI)
      CALL RS3DRT(TMPPETA,SINTH,COSTH,SINPHI,COSPHI)
      CALL RS3DRT(TMPPNU,SINTH,COSTH,SINPHI,COSPHI)

      DO 80 I=1,3
         RETPMU(I) = TMPPMU(I)
         RETPETA(I) = TMPPETA(I)
         RETPNU(I) = TMPPNU(I)
 80   CONTINUE
 1000 continue
      return

      END


