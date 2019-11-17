C#define DEBUG

      PROGRAM NEUT

      IMPLICIT NONE
C
C ++ GAMMA,PI- SCATTERING SIMULATION
C
*     ( purpose )
*        Photoproduction of pions from nuclei
*
*     ( Usage )
*        ./gampi neut.card output.hb
*
*     ( Input )
C         neut.card parameters:
C           EVCT-NEVT  : Number of photons
c           EVCT-MPV   : 1: Fixed energy, 2: Random energy
C           EVCT-PV    : Energy (or energy range for #2 above)
C         NEUT-NUMBNDP : # of bound protons (6: carbon, 8: oxygen, etc.)
C         NEUT-NUMATOM : Number of atoms in nucleus (12: carbon, etc.)
C         NEUT-PFSURF  : Fermi surface momentum for Pauli Blocking (GeV/c)
C           EVCT-RAD   : abs(Separation Energy) (MeV), set EVCT-MPOS = 2
C         NEFRMFLG, NEPAUFLG, NENEFO16, NUCRESCAT : flags
*
*     (Creation Date and Author)
*         2010.01.21 ; P.de Perio
*

#include "necard.h"
#include "necardev.h"
#include "vcwork.h"
#include "efpion.h"
#include "posinnuc.h"
#include "neutparams.h"
#include "fsihist.h"

      REAL*4    POS(3),DIR(3),PNEU(3)
      REAL*4    ENEU

      REAL*4    RLU
      external  RLU

      DOUBLE PRECISION EVPIPROB
      external EVPIPROB

      INTEGER*4 I,J,IDUM,IERR,N1,N2,IVERT,IVCVERT

      character*127 out_fname,strpid

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL*4 COST,PHI

      REAL*4 RMAX,RAD,PTOT

      INTEGER*4 IDPION,IDNUCI,IDNUCF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      INTEGER*4  ICYCLE
      INTEGER*4  LRECL
      PARAMETER (LRECL=1024)

      INTEGER*4  ISIZE
      PARAMETER (ISIZE=5000000)

      REAL*4     H(ISIZE)
      COMMON/PAWC/H

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C
C     Vector Parameter COMMON
C
      INTEGER*4   NPARVEC,IORGVEC(100),IPVEC(100),ICHVEC(100),
     $            IFLVEC(100),IMODE
      REAL*4      ABSPVEC(100),PVEC(3,100),POSV(3)
      COMMON/COM2/NPARVEC,IPVEC,IORGVEC,ICHVEC,IFLVEC,
     $            ABSPVEC,PVEC,POSV,WEIGHT,IMODE

C     Particle Index :     Definition
C            1       :  Incoming Photon (lab frame)
C            2       :  Target Nucleon  (lab frame)
C            3       :  Absorbed Photon (in 2's rest frame)
C            4       :  Recoil Nucleon  (lab frame)
C            5       :  Outgoing Pion   (lab frame)

C
C     FSI history COMMON
C
      INTEGER*4   oNVERT, oNVCVERT,oIFLGVERT(MAXVERT),
     $     oIPVERT(MAXVCVERT), oIVERTI(MAXVCVERT), oIVERTF(MAXVCVERT)
      REAL*4  oPOSVERT(3,MAXVERT), oDIRVERT(3,MAXVCVERT),
     $        oABSPVERT(MAXVCVERT), oABSTPVERT(MAXVCVERT),oFSIPROB

      COMMON /COM3/ oNVERT,oPOSVERT,oIFLGVERT
      COMMON /COM4/ oNVCVERT, oDIRVERT,
     $    oABSPVERT, oABSTPVERT,
     $    oIPVERT,oIVERTI,oIVERTF,oFSIPROB


C
C     Nucleus Target COMMON
C
      INTEGER*4   oNUMBNDN,oNUMBNDP,oNUMFREP,oNUMATOM
      COMMON /COM5/ oNUMBNDN,oNUMBNDP,oNUMFREP,oNUMATOM

ccccccccccccccccccccccccccccccccccccccc

      REAL PIN(3),PINABS,FSURF
      REAL EGAM,RADIUS,S,PGAM(3),PNI(3),PNO(3),PF,NUMERATOR
      REAL TOTENUCL,BETAGAM(3),INVNUCDIR(3),ABSBETAGAM,GAMMANUCL
      REAL ABSCM,GAMMA,BETA,BSTVEC(3),ROTPHI,ROTCOS
      REAL PGAMRST(3),PPIABS,PITTA,PIPHI,PPICM(3)
      REAL COSPI,PGAMABS,WEIGHT

      DOUBLE PRECISION EFABRHO
      EXTERNAL EFABRHO

      REAL GAMFCT
      EXTERNAL GAMFCT

      REAL EFGAMPIMPT
      EXTERNAL EFGAMPIMPT
      REAL EFGAMPIPNT
      EXTERNAL EFGAMPIPNT
      REAL EFGAMPI0PT
      EXTERNAL EFGAMPI0PT
      REAL EFGAMPI0NT
      EXTERNAL EFGAMPI0NT

      REAL EFFRMGAS
      EXTERNAL EFFRMGAS

      real sumTcrs,sumNcrs,sumPcrs

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C-- Initialize ntuple

      CALL HLIMIT(ISIZE)

      out_fname='sample.hb'
      call getarg(2,out_fname)

      CALL HROPEN(85,'CWNT',out_fname,'N',LRECL,IERR)

      IF (IERR.ne.0) THEN
         WRITE(6,*) 'HROPEN:FAILED TO OPEN FILE ', out_fname
         STOP
      ENDIF

      CALL HBNT(1,'NEUTVEC',' ')
      CALL HBNAME(1,'VECINFO',NPARVEC,
     $'NPVC[0,100],IPVC(NPVC),IORGVC(NPVC),ICHVC(NPVC),IFLVC(NPVC),
     $ABSPVC(NPVC),PVC(3,NPVC),POSVC(3),WEIGHT,IMODE')

C--Read card files

      CALL NECARD
      CALL NECARDEV

C--Main loop

      RMAX = C * 2.5

      DO 100 I=1,NEVTEVCT

      IERR = 0
      IMODE = 0

 105  CONTINUE

C--Set interaction position
      POSV(1)=-RMAX+RLU(IDUM)*RMAX*2.
      POSV(2)=-RMAX+RLU(IDUM)*RMAX*2.
      POSV(3)=-RMAX+RLU(IDUM)*RMAX*2.
      RAD = SQRT(POSV(1)**2 + POSV(2)**2 + POSV(3)**2)
      IF (RAD.GT.RMAX) GOTO 105

C     Weight by nuclear density
      IF(EFABRHO(RAD).LT.RLU(IDUM)) GOTO 105

C--Set photon energy
      if (MPVEVCT.eq.1) then
         EGAM=PVEVCT(1)
      else if (MPVEVCT.eq.2) then
         EGAM=(PVEVCT(2)-PVEVCT(1))*RLU(IDUM)+PVEVCT(1)
      else
         write(*,*) 'NEUT : MPV is not set to 1 or 2.'
         stop
      endif


C--Photon Direction in z
      PGAM(1)=0.
      PGAM(2)=0.
      PGAM(3)=EGAM

C     Clear particle stack
      CALL VCCLCM

C     Fill initial photon
      NVC      = 1
      IPVC(1)  = 22
      PVC(1,1) = PGAM(1)
      PVC(2,1) = PGAM(2)
      PVC(3,1) = PGAM(3)

      POSNUC(1,1)=posv(1)
      POSNUC(2,1)=posv(2)
      POSNUC(3,1)=posv(3)

      IORGVC(1)= 0
      IFLGVC(1)= -1
      ICRNVC(1)= 0

C     SET INITIAL NUCLEON MOMENTUM

      PF = 0.
      DIR(1) = 0.
      DIR(2) = 0.
      DIR(3) = 0.

C     Consider fermi motion
 200  if (NEFRMFLG.EQ.0) then
         CALL RNFERM(PF,IDUM)

         FSURF=EFFRMGAS(16.,8.,RAD)
         IF(PF.GT.FSURF)GO TO 200

         CALL RNDIR(DIR)
      end if

      PNI(1) = DIR(1)*PF
      PNI(2) = DIR(2)*PF
      PNI(3) = DIR(3)*PF

C     Boost gamma momentum to nucleon rest frame
C     (Use average nucleon mass.)
      TOTENUCL=SQRT(938.9187865**2+PF**2)
      DO J=1,3
         BETAGAM(J)=PNI(J)/TOTENUCL
         INVNUCDIR(J)=-1.0*DIR(J)
      ENDDO

      ABSBETAGAM=SQRT(BETAGAM(1)**2+BETAGAM(2)**2+BETAGAM(3)**2)
      GAMMANUCL=GAMFCT(ABSBETAGAM)

      DO J=1,3
         PGAMRST(J)=PGAM(J)
      ENDDO

      CALL MCVECBST(PGAMRST,0.,INVNUCDIR,GAMMANUCL)

      PGAMABS=SQRT(PGAMRST(1)**2+PGAMRST(2)**2+PGAMRST(3)**2)

C     Select Nucleon Type

C     Assume same density distribution for neutron and proton,
C     But weight by total xsec. on proton or neutron

      sumPcrs = ( efgampipnt(PGAMABS)+efgampi0pt(PGAMABS))*NUMBNDP
      sumNcrs = (efgampi0nt(PGAMABS)+efgampimpt(PGAMABS))*NUMBNDN
      sumTcrs = sumNcrs+sumPcrs

      IF(RLU(IDUM).LE.(sumPcrs/sumTcrs)) then
         IDNUCI = 2212
      ELSE
         IDNUCI = 2112
      end if

C     Store initial nucleon
      NVC      = NVC+1
      IPVC(2)  = IDNUCI
      call mcmass(IPVC(2),AMASVC(2))
      PVC(1,2) = PNI(1)
      PVC(2,2) = PNI(2)
      PVC(3,2) = PNI(3)

      POSNUC(1,2)=posv(1)
      POSNUC(2,2)=posv(2)
      POSNUC(3,2)=posv(3)

      IORGVC(2)= 0
      IFLGVC(2)= -1
      ICRNVC(2)= 0

C     Store absorbed gamma in nucleon rest frame
      NVC      = NVC+1
      IPVC(3)  = 22
      PVC(1,3) = PGAMRST(1)
      PVC(2,3) = PGAMRST(2)
      PVC(3,3) = PGAMRST(3)

      POSNUC(1,3)=posv(1)
      POSNUC(2,3)=posv(2)
      POSNUC(3,3)=posv(3)

      IORGVC(3)= 1
      IFLGVC(3)= -1
      ICRNVC(3)= 0

C     Event weight is just the total cross section
      weight = sumTcrs

C     Below threshold or above valid energy range
      if (weight.LE.0.) THEN
         IMODE = 0
         GOTO 101
      ENDIF

C     Determine pion direction/type, nucleon type,
C     and event weight by gam+p/n total cross section
 201  CALL EFGAMPIPROD(PGAMABS,IDNUCI,IDPION,IDNUCF,
     &                 COSPI,IMODE,IERR)

C     Below threshold or above valid energy range
      if (IERR.NE.0) THEN
         IMODE = 0
         GOTO 101
       ENDIF

C     Mass of recoil nucleon and pion
      call mcmass(IDNUCF,AMASVC(4))
      call mcmass(IDPION,AMASVC(5))

C     Calculate pion momentum in C.M. frame
      s = (AMASVC(4)+VNUINI*1000)**2 + 2*PGAMABS*(AMASVC(4)+VNUINI*1000)

      NUMERATOR = s**2 - 2*s*(AMASVC(4)**2+AMASVC(5)**2)
     &               + AMASVC(4)**4 + AMASVC(5)**4
     &               - 2*(AMASVC(4)**2)*(AMASVC(5)**2)

C     Threshold for reaction
      if (NUMERATOR.LT.0) THEN
         IMODE = -9
         GOTO 101
      ENDIF

      PINABS = SQRT( NUMERATOR / (4*s) )

      PITTA=ACOS(COSPI)

      PPICM(3)=PINABS*COSPI
      PPICM(1)=PINABS*SIN(PITTA)

      PIPHI=RLU(IDUM)*2*3.14
      PIN(3)=PPICM(3)
      PIN(1)=PPICM(1)*COS(PIPHI)
      PIN(2)=PPICM(1)*SIN(PIPHI)

C     Nucleon recoil in C.M. frame
      PNO(1)=-PIN(1)
      PNO(2)=-PIN(2)
      PNO(3)=-PIN(3)

C     Boost from C.M. frame into lab frame
      ABSCM = SQRT(PNI(1)**2 + PNI(2)**2 + (PNI(3)+PGAM(3))**2)
      BETA = ABSCM/( TOTENUCL + EGAM )
      GAMMA = GAMFCT(BETA)
      BSTVEC(1) = PNI(1)/ABSCM
      BSTVEC(2) = PNI(2)/ABSCM
      BSTVEC(3) = (PNI(3)+PGAM(3))/ABSCM

      CALL MCVECBST(PNO,AMASVC(4),BSTVEC,GAMMA)
      CALL MCVECBST(PIN,AMASVC(5),BSTVEC,GAMMA)

C     Pauli-blocking
      if (NEPAUFLG.EQ.0) then
         IF ((PNO(1)**2+PNO(2)**2+PNO(3)**2).LT.FSURF**2) THEN
            GOTO 201
         end if
      end if

C     Outgoing nucleon
      NVC      = NVC+1
      IPVC(4)  = IDNUCF
      PVC(1,4) = PNO(1)
      PVC(2,4) = PNO(2)
      PVC(3,4) = PNO(3)

      POSNUC(1,4)=posv(1)
      POSNUC(2,4)=posv(2)
      POSNUC(3,4)=posv(3)

      IORGVC(4)= 2
      IFLGVC(4)= 0
      ICRNVC(4)= 1

C     Add generated pion to stack
      NVC      = NVC+1
      IPVC(5)  = IDPION
      PVC(1,5) = PIN(1)
      PVC(2,5) = PIN(2)
      PVC(3,5) = PIN(3)

      POSNUC(1,5)=posv(1)
      POSNUC(2,5)=posv(2)
      POSNUC(3,5)=posv(3)

      IORGVC(5)= 2
      IFLGVC(5)= 0
      ICRNVC(5)= 1


C-- Nuclear effects
      IBOUND = 1

C     --  pion rescattering
      CALL EFCLFSI
      if (NENEFO16.EQ.0) then
         CALL NEFFECT
      end if

C     --  nucleon rescattering (flag is contained)
      CALL NRINTR

C     --  nucleon ejection (flag is contained)
      call absneject

  101 NPARVEC=NVC

      DO 40 N1=1,NVC
         IPVEC(N1)    = IPVC(N1)
         ICHVEC(N1)   = ICRNVC(N1)
         IORGVEC(N1)  = IORGVC(N1)
         IFLVEC(N1)   = IFLGVC(N1)

         ABSPVEC(N1)  = sqrt( PVC(1,N1)**2
     $        +PVC(2,N1)**2
     $        +PVC(3,N1)**2)

         PVEC(1,N1)   = PVC(1,N1)
         PVEC(2,N1)   = PVC(2,N1)
         PVEC(3,N1)   = PVC(3,N1)

 40   CONTINUE

      CALL HFNT(1)

C      if ( MOD(I,100).EQ.0 ) WRITE(*,*) 'Processed event: ',I

  100 CONTINUE

      CALL HROUT(0,ICYCLE,' ')
      CALL HREND('CWNT')

      CALL NECHADIR(COST,PHI,DIR,POSV)

      STOP
      END
