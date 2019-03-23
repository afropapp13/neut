*     -----------------------------------------------------
      SUBROUTINE EFHINELS_ISO(IP,PIN,IPF,PFI,NQN,radius,ireaction)      
*     -----------------------------------------------------
*     
*     ( purpose )
*       CALCULATE INELASTIC SCATTERING at high momentum
*
*     ( input ) 
C          IP   : INITIAL KIND OF PION
C          PIN  : INITIAL PION MOMENTUM(MeV)
C          NQN  : TARGET NUCLEON CHARGE
*
*     ( output ) 
C          PFI  : FINAL PION MOMENTUM
C          PNO  : FINAL NUCLEON MOMENTUM
*
c     ( creation & updates)
C          2010.05: P. de Perio
C          2011.02: This is deprecated
C
C
#include "vcwork.h"
#include "necard.h"
C#include "efpiinth.h"

      DIMENSION PIN(3),PFI(3),PNUC(3),PPI(3),PBUF1(3),PBUF2(3)
      DIMENSION PCM(4,2),PNO(3),DIR(3),BSTVEC(3)
      REAL PF,W,radius,PNOABS
      REAL ABSCM,BETA,GAMMA,TOTENUCL,TOTEPI
      INTEGER*4 IDUM,INUCF,IP,IPF,NQN,ireaction,INUC
      COMMON/EFLAG/IERR
      DATA AMN/939./,AMPI/139./

      write (*,*) "EFHINELS_ISO is deprecated, do not use. Stopping..."
      stop


C---  Initialize
      IPF=IP
      PFI(1)=PIN(1)
      PFI(2)=PIN(2)
      PFI(3)=PIN(3)

      if (NQN.eq.1) then
         INUC = 2212
      else 
         INUC = 2112
      end if
      INUCF=INUC

C     CHARGE EXCHANGE
      IF (ireaction.eq.3) THEN
         CALL EFHINELSCX(IP,INUC,IPF,INUCF)
      END IF


C     -- SET NUCLEON MOMENTUM AND DIRECTION in lab frame
      PF = 0.
      DIR(1) = 0.
      DIR(2) = 0.
      DIR(3) = 0.

      if (NEFRMFLG.EQ.0) then
 200     CALL RNFERM(PF,IDUM)
         
         FSURF=EFFRMGAS(16.,8.,radius)
         IF(PF.GT.FSURF)GO TO 200

         CALL RNDIR(DIR)
      END IF

      PNUC(1) = DIR(1)*PF
      PNUC(2) = DIR(2)*PF
      PNUC(3) = DIR(3)*PF 

      TOTENUCL=SQRT(AMN**2+PF**2)
      TOTEPI=SQRT(AMPI**2+PIN(1)**2+PIN(2)**2+PIN(3)**2)

C     Generate W mass
      W = sqrt( AMN**2 + AMPI**2 + 2*(TOTENUCL*TOTEPI - 
     &     (PIN(1)*PNUC(1)+PIN(2)*PNUC(2)+PIN(3)*PNUC(3)) ) )

c     Isotropic decay in CM frame
 201  CALL GDECA2(W,AMN,AMPI,PCM)

      PNO(1)=PCM(1,1)
      PNO(2)=PCM(2,1)
      PNO(3)=PCM(3,1)

      PFI(1)=PCM(1,2)
      PFI(2)=PCM(2,2)
      PFI(3)=PCM(3,2)

C     Boost from C.M. frame into lab frame
      ABSCM = SQRT((PIN(1)+PNUC(1))**2+(PIN(2)+PNUC(2))**2 
     &     + (PIN(3)+PNUC(3))**2)
      BETA = ABSCM/(TOTENUCL+TOTEPI)
      GAMMA = GAMFCT(BETA)
      BSTVEC(1) = (PIN(1)+PNUC(1))/ABSCM
      BSTVEC(2) = (PIN(2)+PNUC(2))/ABSCM
      BSTVEC(3) = (PIN(3)+PNUC(3))/ABSCM  

      CALL MCVECBST(PNO,AMN,BSTVEC,GAMMA)
      CALL MCVECBST(PFI,AMPI,BSTVEC,GAMMA)

      if (NEPAUFLG.EQ.0) then
         PNOABS = SQRT(PNO(1)**2+PNO(2)**2+PNO(3)**2)
         IF(PNOABS.LT.FSURF)GO TO 201
      END IF



C     ADD NUCLEON TO VCWORK (non-functional, need better way
C                            to keep track of indexing)
C      if (neabspiemit.eq.0) then
C         return
C      endif

C     Do not store if stack is full (save one spot for pion)
C      if (nvc.GE.maxvc-1) then
C         return
C      end if

C      ipvc(nvc+1)  = INUCF
C      icrnvc(nvc+1) = 1
C      iflgvc(nvc+1) = 0
C      ivtivc(nvc+1) = 1
C
C      iorgvc(nvc+1) = LOOPPI
C
C      pvc(1,nvc+1) = PNO(1)
C      pvc(2,nvc+1) = PNO(2)
C      pvc(3,nvc+1) = PNO(3)
C
C      nvc = nvc + 1

      RETURN
      END
