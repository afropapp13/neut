*     -----------------------------------------------------
      SUBROUTINE EFCOHSCT(PPI,POUT)
*     -----------------------------------------------------
*
*     *!!!!!*   OBSOLETE   *!!!!!*
*     !*****!Don't USE THIS!*****! for low energy region
*
*     ( purpose )
*       SIMULATE COHERENT SCATTERING
*
*     ( input ) 
*       PPI  : INITIAL PION MOMENTUM
*
*     ( output ) 
*       POUT : FINAL PION MOMENTUM
C     
#include "necard.h"

      DIMENSION PPI(3),POUT(3),SK(3),PF(3),GK(3),EV(3)
      DATA AMN/939./,AMPI/139./
      DATA Q0/200./
C     

C      WRITE(*,*) "EFCOHSCT: THIS SUBROUTINE IS OBSOLETE!"

      if(NUMATOM.gt.1) then
         CALL RNFERM(PFABS,IWAVE)
         CALL RNDIR(PF)
      else
         PFABS = 0
         PF(1) = 0
         PF(2) = 0
         PF(3) = 0
      end if
C     
      PF(1)=PFABS*PF(1)
      PF(2)=PFABS*PF(2)
      PF(3)=PFABS*PF(3)
      GK(1)=PF(1)+PPI(1)
      GK(2)=PF(2)+PPI(2)
      GK(3)=PF(3)+PPI(3)
      EPI=SQRT(AMPI**2+PPI(1)**2+PPI(2)**2+PPI(3)**2)
      ENU=SQRT(AMN**2+PF(1)**2+PF(2)**2+PF(3)**2)
      GKABS=SQRT(GK(1)**2+GK(2)**2+GK(3)**2)
      BETA=GKABS/(EPI+ENU)
      GM=1./SQRT(1.-BETA**2)
      IF(GKABS.EQ.0.)GO TO 140
      EV(1)=-GK(1)/GKABS
      EV(2)=-GK(2)/GKABS
      EV(3)=-GK(3)/GKABS
      GO TO 150
  140 EV(1)=1.
      EV(2)=0.
      EV(3)=0.
  150 SK(1)=PPI(1)
      SK(2)=PPI(2)
      SK(3)=PPI(3)
      CALL MCVECBST(SK,AMPI,EV,GM)
C     
      SKABS2=SK(1)**2+SK(2)**2+SK(3)**2
C  200 RAN=RNDM(DUM)
  200 RAN=RLU(IDUM)
      Q2=-Q0**2*ALOG(RAN)
      SIN2T=Q2/4./SKABS2
      SINT=SQRT(SIN2T)
      IF(SINT.GT.1.)GO TO 200
      THETA=ASIN(SINT)
      THETA=THETA*2.
      COST=COS(THETA)
C      PHI=RNDM(DUM)*6.2831852
      PHI=RLU(IDUM)*6.2831852
      CALL NECHADIR(COST,PHI,SK,POUT)
      EV(1)=-EV(1)
      EV(2)=-EV(2)
      EV(3)=-EV(3)
      CALL MCVECBST(POUT,AMPI,EV,GM)
C     
C     To do: Knock out the nucleon (maybe) -Patrick
C
      RETURN
      END
