*     -----------------------------------------------------
      SUBROUTINE EFINELS(IP,PIN,IPF,PFI,radius)
*     -----------------------------------------------------
*
*     ( purpose )
*       CALCULATE INELASTIC SCATTERING
*
*     ( input )
C          IP   : INITIAL KIND OF PION
C          PIN  : INITIAL PION MOMENTUM(MeV)
C          r    : radius(fm)
*
*     ( output )
C          IPF  : FINAL KIND OF PION
C          PFI  : FINAL PION MOMENTUM
*
C     ( update )
C          2009.12.09: Added ejected nucleon - P. de Perio
C
#include "necard.h"
#include "neutparams.h"
#include "vcwork.h"
C#include "efpiinth.h"

      DIMENSION PIN(3),PFI(3),PNUC(3),PPI(3),PBUF1(3),PBUF2(3)
      COMMON/EFLAG/IERR
      INTEGER*4 INUCF
      DIMENSION PNUCF(3)
C
C     -- SET PION DIRECTION TO X-AXIS
      PPI(1)=SQRT(PIN(1)**2+PIN(2)**2+PIN(3)**2)
      PPI(2)=0.
      PPI(3)=0.
      COST=PIN(3)/PPI(1)
      THETA=ACOS(COST)
      IF(PIN(1).EQ.0. .AND. PIN(2).EQ.0.)GO TO 10
      PHI=ATAN2(PIN(2),PIN(1))
      GO TO 20
   10 PHI=0.
C     -- SET NUCLEON MOMENTUM AND DIRECTION
   20 CALL EFSETPF(IP,PPI,PNUC,MODE,radius)
      IF(IERR.NE.0)GO TO 9000
C     -- SET PION MOMENTUM
      CALL EFSCATT(IP,PPI,PNUC,MODE,IPF,PFI,INUCF,PNUCF,radius)
C     -- BACK TO INITIAL COORDINATE
C     - RANDOM ROTATION AROUND X-AXIS
C      ETA=RNDM(DUM)*6.283185
      ETA=RLU(IDUM)*6.283185

C     ROTATE PION VECTOR
      PBUF1(1)=PFI(1)
      PBUF1(2)=COS(ETA)*PFI(2)-SIN(ETA)*PFI(3)
      PBUF1(3)=SIN(ETA)*PFI(2)+COS(ETA)*PFI(3)
C     - THETA ROTATION
      PBUF2(1)=SIN(THETA)*PBUF1(1)-COS(THETA)*PBUF1(3)
      PBUF2(2)=PBUF1(2)
      PBUF2(3)=COS(THETA)*PBUF1(1)+SIN(THETA)*PBUF1(3)
C     - PHI ROTATION
      PFI(1)=COS(PHI)*PBUF2(1)-SIN(PHI)*PBUF2(2)
      PFI(2)=SIN(PHI)*PBUF2(1)+COS(PHI)*PBUF2(2)
      PFI(3)=PBUF2(3)


C     ROTATE NUCLEON VECTOR
      PBUF1(1)=PNUCF(1)
      PBUF1(2)=COS(ETA)*PNUCF(2)-SIN(ETA)*PNUCF(3)
      PBUF1(3)=SIN(ETA)*PNUCF(2)+COS(ETA)*PNUCF(3)
C     - THETA ROTATION
      PBUF2(1)=SIN(THETA)*PBUF1(1)-COS(THETA)*PBUF1(3)
      PBUF2(2)=PBUF1(2)
      PBUF2(3)=COS(THETA)*PBUF1(1)+SIN(THETA)*PBUF1(3)
C     - PHI ROTATION
      PNUCF(1)=COS(PHI)*PBUF2(1)-SIN(PHI)*PBUF2(2)
      PNUCF(2)=SIN(PHI)*PBUF2(1)+COS(PHI)*PBUF2(2)
      PNUCF(3)=PBUF2(3)

C     ADD NUCLEON TO VCWORK(non-functional, need better way
C                            to keep track of indexing)
C      if (neabspiemit.eq.0) then
C         return
C      endif
C
C     Do not store if stack is full (save one spot for pion)
C      if (nvc.GE.maxvc-1) then
C         return
C      end if
C
C      ipvc(nvc+1)  = INUCF
C      icrnvc(nvc+1) = 1
C      iflgvc(nvc+1) = 0
C      ivtivc(nvc+1) = 1
C
C      iorgvc(nvc+1) = LOOPPI
C
C      pvc(1,nvc+1) = PNUCF(1)
C      pvc(2,nvc+1) = PNUCF(2)
C      pvc(3,nvc+1) = PNUCF(3)
C
C      nvc = nvc + 1

      RETURN
C
C     -- ERROR RETURN (CAN NOT SCATTER)
C
 9000 IPF=IP
      PFI(1)=PIN(1)
      PFI(2)=PIN(2)
      PFI(3)=PIN(3)
      RETURN
      END
