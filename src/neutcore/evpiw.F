**********************************************************************
*     -------------------------
      DOUBLE PRECISION FUNCTION EVPIW()
*     -------------------------
*
*     ( purpose )
*        SAVE INTERACTION CROSS-SECTION
*
*     ( input )
*       COMMON NEWORK
*          modene
*          numne (=4)
*          ipne
*          pne
*       COMMON NECARD
*       COMMON NEUTMODEL
*       COMMON NEUTPARAMS
*
*     ( output )
*       DIFCRSSECT
*
*     ( creation date and author )
*       2010.07.19 ; Taken from nemodsel
*       2010.10.13 ; Integrated into NEUTCORE - P. de Perio
*
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C    MODSEL   : INTERACTION MODE
C
C
C        ######### NEUTRINO MODE #########
C
C            ***** CHARGED CURRENT *****
C
C               -- ELASTIC --
C           1 : NEU,N --> LEPTON-,P
C
C               -- SINGLE PI FROM DELTA RESONANCE --
C          11 : NEU,P --> LEPTON-,P,PI+
C          12 : NEU,N --> LEPTON-,P,PI0
C          13 : NEU,N --> LEPTON-,N,PI+
C
C          16 : NEU,O(16) --> LEPTON-,O(16),PI+
C
C               -- SINGLE GAMMA FROM DELTA RESONANCE --
C          17 : NEU,N --> LEPTON-,P,GAMMA
C
C               -- MULTI PI (1.3 < W < 2.0 GeV) --
C          21 : NEU,(N OR P) --> LEPTON-,(N OR P),MULTI PI
C
C               -- SINGLE ETA FROM DELTA RESONANCE --
C                                     (added 97/12/01 J.Kameda)
C          22 : NEU,N --> LEPTON-,P,ETA0
C
C               -- SINGLE K FROM DELTA RESONANCE --
C                                     (added 98/02/25 J.Kameda)
C          23 : NEU,N --> LEPTON-,LAMBDA,K+
C
C               -- DEEP INELASTIC (2.0 GeV < W , JET set) --
C          26 : NEU,(N OR P) --> LEPTON-,(N OR P),MESONS
C
C            ***** NEUTAL CURRENT *****
C
C               -- SINGLE PI FROM DELTA RESONANCE --
C          31 : NEU,N --> NEU,N,PI0
C          32 : NEU,P --> NEU,P,PI0
C          33 : NEU,N --> NEU,P,PI-
C          34 : NEU,P --> NEU,N,PI+
C
C          36 : NEU,O(16) --> NEU,O(16),PI0
C
C               -- SINGLE GAMMA FROM DELTA RESONANCE --
C          38 : NEU,N --> NEU,N,GAMMA
C          39 : NEU,P --> NEU,P,GAMMA
C
C               -- MULTI PI (1.3 GeV < W < 2.0 GeV) --
C          41 : NEU,(N OR P) --> NEU,(N OR P),MULTI PI
C
C               -- SINGLE ETA FROM DELTA RESONANCE --
C                                     (added 97/12/01 J.Kameda)
C          42 : NEU,N --> NEU,N,ETA0
C          43 : NEU,P --> NEU,P,ETA0
C
C               -- SINGLE  K  FROM DELTA RESONANCE --
C                                     (added 98/02/20 J.Kameda)
C          44 : NEU,N --> NEU,LAMBDA,K0
C          45 : NEU,P --> NEU,LAMBDA,K+
C
C               -- DEEP INELASTIC (2.0 GeV < W , JET set) --
C          46 : NEU,(N OR P) --> NEU,(N OR P),MESONS
C
C               -- ELASTIC --
C          51 : NEU,P --> NEU,P
C          52 : NEU,N --> NEU,N
C
C
C        ######### ANTI NEUTRINO MODE #########
C
C            ***** CHARGED CURRENT *****
C
C               -- ELASTIC --
C          -1 : NEUBAR,P --> LEPTON+,N
C
C               -- SINGLE PI FROM DELTA RESONANCE --
C         -11 : NEUBAR,N --> LEPTON+,N,PI-
C         -12 : NEUBAR,P --> LEPTON+,N,PI0
C         -13 : NEUBAR,P --> LEPTON+,P,PI-
C
C         -16 : NEUBAR,O(16) --> LEPTON+,O(16),PI-
C
C               -- SINGLE GAMMA FROM DELTA RESONANCE --
C         -17 : NEUBAR,P --> LEPTON+,N,GAMMA
C
C               -- MULTI PI (W > 1.4 GEV) --
C         -21 : NEUBAR,(N OR P) --> LEPTON+,(N OR P),MULTI PI
C
C               -- SINGLE ETA FROM DELTA RESONANCE --
C                                     (added 97/12/01 J.Kameda)
C         -22 : NEUBAR,P --> LEPTON+,N,ETA0
C
C               -- SINGLE  K FROM DELTA RESONANCE --
C                                     (added 98/02/25 J.Kameda)
C         -23 : NEUBAR,P --> LEPTON+,LAMBDA,K0
C
C               -- DEEP INELASTIC (2.0 GeV < W , JET set) --
C         -26 : NEUBAR,(N OR P) --> LEPTON+,(N OR P),MESONS
C
C               ** NEUTAL CURRENT **
C
C               -- SINGLE PI FROM DELTA RESONANCE --
C         -31 : NEUBAR,N --> NEUBAR,N,PI0
C         -32 : NEUBAR,P --> NEUBAR,P,PI0
C         -33 : NEUBAR,N --> NEUBAR,P,PI-
C         -34 : NEUBAR,P --> NEUBAR,N,PI+
C
C         -36 : NEUBAR,O(16) --> NEUBAR,O(16),PI0
C
C               -- SINGLE GAMMA FROM DELTA RESONANCE --
C         -38 : NEUBAR,N --> NEUBAR,N,GAMMA
C         -39 : NEUBAR,P --> NEUBAR,P,GAMMA
C
C               -- MULTI PI (W > 1.4 GEV) --
C         -41 : NEUBAR,(N OR P) --> NEUBAR,(N OR P),MULTI PI
C
C               -- SINGLE ETA FROM DELTA RESONANCE --
C                                     (added 97/12/01 J.Kameda)
C         -42 : NEUBAR,N --> NEUBAR,N,ETA0
C         -43 : NEUBAR,P --> NEUBAR,P,ETA0
C
C               -- SINGLE  K  FROM DELTA RESONANCE --
C                                     (added 98/02/20 J.Kameda)
C          -44 : NEUBAR,N --> NEUBAR,LAMBDA,K0
C          -45 : NEUBAR,P --> NEUBAR,LAMBDA,K+
C
C               -- DEEP INELASTIC (2.0 GeV < W , JET set) --
C         -46 : NEUBAR,(N OR P) --> NEUBAR,(N OR P),MESONS
C
C               -- ELASTIC --
C         -51 : NEUBAR,P --> NEUBAR,P
C         -52 : NEUBAR,N --> NEUBAR,N
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C    IPAR : NEUTRINO TYPE
C        = 12  NUE
C        =-12  NUEBAR
C        = 14  NUMU
C        =-14  NUMUBAR
C    E    : ENERGY OF NEUTRINO (UNIT : GEV)
C
C    MEANING OF DATA
C    RATN   : NEUTRON RATIO IN TARGET
C    RATP   : PROTON RATIO IN TARGET
C    RATFRE : FREE PROTON RATIO IN TARGET
C    IFREE  : =0 BOUND NUCLEON, =1 FREE PROTON
C
C    IBOUND : =0 Free proton  , =1 Bound nucleon
C
      IMPLICIT NONE

#include "nework.h"
#include "rscons.h"
#include "neutcrs.h"
#include "neutmodel.h"
#include "necard.h"
#include "vcwork.h"
#include "posinnuc.h"

      integer*4 ipar,i
      real*4    e
C
C -- Preparation
C
      real*4 plep,   elep,   xmlep
      real*4 xmn01,  xmn02,  enuc1, enuc2
      real*4 q2
      real*4 w

      real yfn

C     Adler angles needed for Minoo's model
C     Also needed for reweighting Rein-Seghal ejection method
      real*4 costhetaAd, phiAd

C -- For boosting
      real*4 PNEU(3),PFABS,AMIN,BETA,EV(3),GM

      integer*4 IMOD

      evpiw = 1

C     Only works for Rein-Sehgal model
      IF (MDLSPI .NE. 1) THEN
        write(*,*) "ERROR***"
        write(*,*) "Pion/nucleon ejection dial MDLSPIEJ reweight"
        write(*,*) "only works for Rein-Sehgal model. Returning 1"
        GOTO 999
      endif

C     Santiy check the input
      IF (MDLSPIEJ .LT. 0 .OR. MDLSPIEJ .GT. 3) THEN
        write(*,*) "ERROR***"
        write(*,*) "Pion/nucleon ejection dial MDLSPIEJ reweight:", MDLSPIEJ
        write(*,*) "only works for value 0,1,2,3"
        GOTO 999
      endif

      ipar = ipne(1)

C     Only applies to resonant interaction
      if (abs(modene) .lt. 11 .or. ! CC1pi+1p starts at 11
     &    (abs(modene) .gt. 14 .and. ! CC1pi+1n ends at 14
     &     abs(modene) .lt. 31) .or. ! NCpi starts at 31
     &    (abs(modene) .gt. 34)) then
        write(*,*) "ERROR***"
        write(*,*) "Pion/nucleon ejection dial MDLSPIEJ reweight"
        write(*,*) "only works for mode 11,12,13,31,32,33,34, not: ", modene

        goto 999
      endif

      IF (MODENE.eq.11)  IMOD=1
      IF (MODENE.eq.12)  IMOD=2
      IF (MODENE.eq.13)  IMOD=3
      IF (MODENE.eq.31)  IMOD=6
      IF (MODENE.eq.32)  IMOD=4
      IF (MODENE.eq.33)  IMOD=7
      IF (MODENE.eq.34)  IMOD=5
      IF (MODENE.eq.-11) IMOD=11
      IF (MODENE.eq.-12) IMOD=12
      IF (MODENE.eq.-13) IMOD=13
      IF (MODENE.eq.-31) IMOD=16
      IF (MODENE.eq.-32) IMOD=14
      IF (MODENE.eq.-33) IMOD=17
      IF (MODENE.eq.-34) IMOD=15

C     Set lepton mass
      if (abs(ipne(3)).eq.11) then
         xmlep=XME
      else if (abs(ipne(3)).eq.13) then
         xmlep=XMMU
      else if (abs(ipne(3)).eq.15) then
         xmlep=XMTAU
      else
         xmlep=0.
      endif

      xmn01 = XMN
      xmn02 = XMN

      e = sqrt(pne(1,1)**2+pne(2,1)**2+pne(3,1)**2)

      enuc1= sqrt((pne(1,2)**2+pne(2,2)**2+pne(3,2)**2)+xmn01**2)
      plep = sqrt((pne(1,3)**2+pne(2,3)**2+pne(3,3)**2))
      elep = sqrt((pne(1,3)**2+pne(2,3)**2+pne(3,3)**2)+xmlep**2)
      enuc2= sqrt((pne(1,4)**2+pne(2,4)**2+pne(3,4)**2)+xmn02**2)

C     Calculate Q2
      q2 =   (e - elep)**2
     $     - ( (pne(1,1)-pne(1,3))**2
     $        +(pne(2,1)-pne(2,3))**2
     $        +(pne(3,1)-pne(3,3))**2)

C     Calculate W
      w = sqrt ( (e + enuc1 - elep)**2
     $          -( (pne(1,1)+pne(1,2)-pne(1,3))**2
     $            +(pne(2,1)+pne(2,2)-pne(2,3))**2
     $            +(pne(3,1)+pne(3,2)-pne(3,3))**2))

C     Check the W
      if ( w.le.XMN+XMPI ) then
        WRITE(*,*) "Resetting W: ", w, " -> ", XMN+XMPI
        w = XMN+XMPI
      else if (w.ge.2) then
        write(*,*) "Resetting W: ", w, " -> ", 2.0
        w = 2.0
      endif

C     Boost neutrino to nucleon rest frame since all the calculations in
C     Rein-Sehgal and Minoo are done in nucleon rest frame
C     Q2 and W neednt be calculated since invariant
      PFABS=SQRT(pne(1,2)**2+pne(2,2)**2+pne(3,2)**2)
      IF (abs(PFABS).gt.1E-5) THEN
C       Get relativistic beta (momentum/energy)
        BETA=PFABS/enuc1
C       Get relativistic gamma factor
        GM=1./SQRT(1.-BETA**2)
        do i =1,3
C         Copy the neutrino energy
          PNEU(i) = pne(i,1)
C         Get the unit vector of the initial state nucleon momentum
          EV(i)=-pne(i,2)/PFABS
        enddo
C       Boost neutrino to nucleon frame
        CALL MCVECBST(PNEU,0.,EV,GM)
        e = SQRT(PNEU(1)**2+PNEU(2)**2+PNEU(3)**2)
      ENDIF

C     For Rein-Sehgal model choice call Rein-Sehgal code
C     IMOD is the interaction mode (e.g. CC1pi+1p), IORIG refers to
C     using lepton mass effect, xmlep is lepton mass, e = neutrino
C     energy, q2 = -Q2, W = hadronic mass, dcp is returned positive
C     helicity differential xsec, dcm is returned negative helicity
C     differential xsec, the 33s are using only the Delta(1232) 
C     Need to calculate theta and phi
C     Get the primary nucleon and pion so we can boost into the frame
C     Get the Adler angles
      costhetaAd = -999
      phiAd = -999
      call mkAdler(costhetaAd, phiAd)
      call rscalcy(IMOD, ipne(3), e, q2, w, acos(costhetaAd), phiAd,
     &             yfn)

      evpiw = yfn

 999  continue
      RETURN
      END

