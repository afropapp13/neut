*     -----------------------------------------------------
      SUBROUTINE EFSETPF(IP,PPI,PFSEL,MODSEL,radius)
*     -----------------------------------------------------
*
*     ( purpose )
*       SET NUCLEON MOMENTUM AND INTERACTION MODE
*
*     ( input )
C         IP     : IITIAL PION KIND
C         PPI    : INITIAL PION MOMENTUM IN LAB. SYSTEM(MeV)
C         radius : radius(fm)
*
*     ( output )
C         PFSEL  : NUCLEON MOMENTUM IN LAB. SYSTEM(MeV)
C         MODSEL : INTEACTION MODE (SEE FUNCTION FSQUAR)
*
*     ( Creation Date and Author )
*       2007.11.05 ; G.Mitsuka - add weight of #of nucleon
*
#include "necard.h"
      DIMENSION PPI(3),PFSEL(3),MODE(5),AMP(5)
      DIMENSION XMAP(7,6,4),SK(3),GK(3),PF(3),EV(3),SUMMAP(4)
      COMMON/EFLAG/IERR
      DATA AMN/939./,AMPI/139./

      INTEGER IPPION
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
               IF (NEFKINVER.eq.0) THEN
                  XMAP(I+1,J,K)=AMP(K)*EFFRMMOM(PFABS)
               else IF (NEFKINVER.eq.1) THEN
                  XMAP(I+1,J,K)=AMP(K)*EFFRMMOM_L(radius,PFABS)
               else
                  write(*,*) 'Invalid NEFKINVER = ',NEFKINVER
                  stop
               endif
 160        CONTINUE
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
C  230 HIT=SUMMAP(NMODE)*RNDM(DUM)
  230 HIT=SUMMAP(NMODE)*RLU(IDUM)
      K=0
  240 K=K+1
      IF(HIT.LT.SUMMAP(K))GO TO 250
      IF(K.EQ.NMODE)GO TO 9000
      GO TO 240
  250 MODSEL=MODE(K)
C     -- SELECT FERMI MOMENTUM AND SCATTERING ANGLE
      if (NEFRMFLG.EQ.1) then
         PFSEL(1)=0
         PFSEL(2)=0
         PFSEL(3)=0
         return
      endif

      XMAX=0.
      DO 300 I=2,6
         DO 310 J=1,6
            IF(XMAX.GT.XMAP(I,J,K))GO TO 310
            XMAX=XMAP(I,J,K)
  310    CONTINUE
  300 CONTINUE
C  320 PFABS=RNDM(DUM)*275.
  320 PFABS=RLU(IDUM)*275.
C      COST=-1.+RNDM(DUM)*2.
      COST=-1.+RLU(IDUM)*2.
C      HIT=XMAX*RNDM(DUM)
      HIT=XMAX*RLU(IDUM)
      PFA=PFABS/50.
      COA=(COST+1.)/0.4
      NP=IFIX(PFA)+1
      NC=IFIX(COA)+1
      IF(NP.GE.7 .OR. NC.GE.6) THEN 
C         WRITE(6,600)
C 600     FORMAT(' *** ERROR IN SETPF ***')
         goto 330
      endif
      DPF=PFA+1.-FLOAT(NP)
      DCO=COA+1.-FLOAT(NC)
      DAT1=(1.-DPF)*XMAP(NP,NC,K)+DPF*XMAP(NP+1,NC,K)
      DAT2=(1.-DPF)*XMAP(NP,NC+1,K)+DPF*XMAP(NP+1,NC+1,K)
      DISTRI=(1.-DCO)*DAT1+DCO*DAT2
      IF(HIT.GT.DISTRI)GO TO 320
      if (abs(cost).ge.1.) then
         sint = 0.
      else
         SINT=SQRT(1.-COST**2)
      endif
      PFSEL(1)=PFABS*COST
      PFSEL(2)=PFABS*SINT
      PFSEL(3)=0.
      RETURN
C
C     -- ERROR RETURN
C
 9000 WRITE(6,900)HIT,SUMMAP,PPI,radius
  900 FORMAT(' ??? WARNING AT SETPF ??? HIT=',G15.7/
     &                 ' SUMMAP=',4G15.7,'  PPI=',3G15.7,
     &'  radius=',G15.7)
      IERR=1
      RETURN
      END
