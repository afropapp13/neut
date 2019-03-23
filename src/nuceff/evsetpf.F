*     -----------------------------------------------------
      SUBROUTINE EVSETPF(IP,PPI,radius,INTMOD,PROB)      
*     -----------------------------------------------------
*     
*     ( purpose )
*       Calculate probability of CX or non-CX interaction 
*       for evpiprob.F.
*       Note: Must be kept in sync with efsetpf.F
*
*     ( input ) 
C         IP     : IITIAL PION KIND
C         PPI    : INITIAL PION MOMENTUM IN LAB. SYSTEM(MeV)
C         radius : radius(fm)
C         INTMOD : INTEACTION MODE = 4(CX) or 8(non-CX)
C
*     ( output ) 
C         PROB   : Probability 
*
*     ( Creation Date and Author )
*       2007.11.05 ; G.Mitsuka - add weight of #of nucleon
*
#include "necard.h"
      DIMENSION PPI(3),MODE(5),AMP(5)
      DIMENSION XMAP(7,6,4),SK(3),GK(3),PF(3),EV(3),SUMMAP(4)
      COMMON/EFLAG/IERR
      DATA AMN/939./,AMPI/139./

      REAL*4 PROB

      INTEGER IPPION, INTMOD
      REAL PNRAT(4,3)
C     PNRAT is applied to each interaction as an weight
C     pi+
      PNRAT(1,1) = float(NUMBNDP)/float(NUMATOM)
      PNRAT(2,1) = float(NUMBNDN)/float(NUMATOM)
      PNRAT(3,1) = float(NUMBNDN)/float(NUMATOM)
      PNRAT(4,1) = 0.
C     pi0
      PNRAT(1,2) = float(NUMBNDP)/float(NUMATOM)
      PNRAT(2,2) = float(NUMBNDN)/float(NUMATOM)
      PNRAT(3,2) = float(NUMBNDP)/float(NUMATOM)
      PNRAT(4,2) = float(NUMBNDN)/float(NUMATOM)
C     pi-
      PNRAT(1,3) = float(NUMBNDP)/float(NUMATOM)
      PNRAT(2,3) = float(NUMBNDN)/float(NUMATOM)
      PNRAT(3,3) = float(NUMBNDP)/float(NUMATOM)
      PNRAT(4,3) = 0.
C     
      IERR=0
      PROB=1
C     
      IF(IP.EQ.211)GO TO 10
      IF(IP.EQ.111)GO TO 20
      IF(IP.EQ.-211)GO TO 30
      RETURN
C     
   10 NMODE=3
      MODE(1)=1
      MODE(2)=2
      MODE(3)=3
      IPPION=1
      GO TO 100
   20 NMODE=4
      MODE(1)=4
      MODE(2)=5
      MODE(3)=6
      MODE(4)=7
      IPPION=2
      GO TO 100
   30 NMODE=3
      MODE(1)=8
      MODE(2)=9
      MODE(3)=10
      IPPION=3
C     -- MAKE NUCLEON MOMENTUM AND NUC-PI OPENING ANGLE MAP
  100 CONTINUE
      DO 110 K=1,NMODE
         DO 110 J=1,6
            XMAP(1,J,K)=0.
            XMAP(7,J,K)=0.
  110 CONTINUE
C     
      DO 120 I=1,5
         PFABS=FLOAT(I)*50.
         DO 130 J=1,6
            COST=-1.+0.4*FLOAT(J-1)
            if ((abs(cost).ge.1.).or.(J.eq.1).or.(J.eq.6)) then
               sint = 0.
            else
               SINT=SQRT(1.-COST**2)
            endif
            PF(1)=PFABS*COST
            PF(2)=PFABS*SINT
            PF(3)=0.
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
  140       EV(1)=1.
            EV(2)=0.
            EV(3)=0.
  150       SK(1)=PF(1)
            SK(2)=PF(2)
            SK(3)=PF(3)
            CALL MCVECBST(SK,AMN,EV,GM)
            CALL EFCALAMP(SK,GK,PPI,radius)
            CALL EFPFCOST(PF,PPI,SK,GK,EV,GM,NMODE,MODE,AMP,radius)
            DO 160 K=1,NMODE
  160       XMAP(I+1,J,K)=AMP(K)*EFFRMMOM(PFABS)
  130    CONTINUE
  120 CONTINUE
C     -- SELECT MODE
      DO 200 K=1,NMODE
         SUMMAP(K)=0.
         DO 210 I=2,6
            DO 210 J=1,6
C     Weight XMAP(i,j,mode) by PNRAT(mode,ippion)
  210    SUMMAP(K)=SUMMAP(K)+XMAP(I,J,K)*PNRAT(K,IPPION)
         IF(K.EQ.1)GO TO 200
         SUMMAP(K)=SUMMAP(K)+SUMMAP(K-1)
  200 CONTINUE


C     Calculate probability for CX/non-CX modes

C     Charged pion
      IF (NMODE.eq.3) then
C        CX (see effsquar.F)
         IF (INTMOD.eq.4) then
            PROB = (SUMMAP(3)-SUMMAP(2))/SUMMAP(NMODE)
C        non-CX
         else if (INTMOD.eq.8) then
            PROB = SUMMAP(2)/SUMMAP(NMODE)
         else
            write(*,*) "EVSETPF: Warning unknown mode:", INTMOD
         end if
C     pi0
      ELSE IF (NMODE.eq.4) then
C        CX (see effsquar.F)
         IF (INTMOD.eq.4) then
            PROB = (SUMMAP(4)-SUMMAP(2))/SUMMAP(NMODE)
C        non-CX
         else if (INTMOD.eq.8) then
            PROB = SUMMAP(2)/SUMMAP(NMODE)
         else
            write(*,*) "EVSETPF: Warning unknown mode:", INTMOD
         end if
      else
         write(*,*) "EVSETPF: Warning unknown particle:", IP
      end if

      RETURN
      END
