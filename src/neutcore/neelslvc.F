************************************************************************
*     -------------------------------------------
      SUBROUTINE NEELSLVC(ENEU,PF,IPAR,PLEP,IERR)
*     -------------------------------------------
*
*     (Purpose)
*       set final lepton momentum vector for elastic scattering
*
*     (Input)
*       ENE    : NEUTRINO ENERGY ( GEV )
*       PF(3)  : FERMI MOMENTUM ( GEV/C )      
*       IPAR   : NEUTRINO TYPE
*                  12 : NEU-E
*                 -12 : NEU-E-BAR            
*                  14 : NEU-MU
*                 -14 : NEU-MU-BAR            
*                  16 : NEU-TAU
*                 -16 : NEU-TAU-BAR
*     
*     (Output)
*       PLEP(3) : FINAL LEPTON MOMENTUM ( GEV/C )
*       IERR    : ERROR CODE      
*
*     (Creation Date and Author)
*       1984.01.08 ; M.Nakahata
*       1987.06.23 ; N.SATO     INCLUDE TAU
*       1988.11.28 ; T.KAJITA   NUCLEAR VINDING ENERGY IS CONSIDERD
*       1995.01.30 ; K. KANEYUKI  FOR S.K.
*                    ELSLVC -> NEELSLVC
*       1996.01.11 ; K.Kaneyuki  NEPAUFLG=1 -> No Pauli Effect
*     
*       1998.02.23 ; K.Ishihara set PFERMI 217 MeV
*
************************************************************************
#include "necard.h"

      REAL PF(3),PNEU(3),EV(3),DNEU(3),DLEP(3),PNP(3),PLEP(3)
      REAL XMAP(51,11)
      REAL*8 DNELSQ2, DBLE0, DBLQ2

      EXTERNAL DNELSQ2, RLU, FNQ2MAX

      DATA AM/0.938/
      DATA PFERMI/0.217/


C--- OBSOLETE FUNCTION

      write(6,*) '*** ERROR NEELSLVC IS OBSOLETE ***'
      write(6,*) '***       USE  NEELSLPV        ***'

      STOP

C=======
C     
      IERR=0
C
C -- CALCUALTE NEUTRINO ENERGY IN NUCELON STOP SYSTEM
C

      DO 100 LOOP=0,1000000

         PNEU(1)=ENEU
         PNEU(2)=0.
         PNEU(3)=0.

C-- set pf (not considering pauli blocking or potential)
         call neelsspff(eneu,ipar,pf,LOOP,ierr)
         if (ierr.ne.0) goto 9000

         PFABS=SQRT(PF(1)**2+PF(2)**2+PF(3)**2)
         E=SQRT(PFABS**2+AM**2)
         BETA=PFABS/E
         GM=1./SQRT(1.-BETA**2)
         IF (PFABS.NE.0.) THEN
            EV(1)=-PF(1)/PFABS
            EV(2)=-PF(2)/PFABS
            EV(3)=-PF(3)/PFABS
         ELSE
            EV(1)=1.
            EV(2)=0.
            EV(3)=0.
         ENDIF
         CALL MCVECBST(PNEU,0.,EV,GM)
         ENSTOP=SQRT(PNEU(1)**2+PNEU(2)**2+PNEU(3)**2)
         IF(ENSTOP.LT.0.001 .AND. IABS(IPAR).EQ.12)GO TO 9000
         IF(ENSTOP.LT.0.112 .AND. IABS(IPAR).EQ.14)GO TO 9000
         IF(ENSTOP.LT.3.48109 .AND. IABS(IPAR).EQ.16)GO TO 9000
         DNEU(1)=PNEU(1)/ENSTOP
         DNEU(2)=PNEU(2)/ENSTOP
         DNEU(3)=PNEU(3)/ENSTOP
C     
C     -- MAKE (Q**2 -- SCATTERING ANGLE) MAP
C     
         CALL VZERO(XMAP,561)
         
         EV(1)=-EV(1)
         EV(2)=-EV(2)
         EV(3)=-EV(3)
         Q2MAXX=FNQ2MAX(ENSTOP,IPAR)
         DQ2=Q2MAXX/50.
         IPLEP=IABS(IPAR)-1
         CALL MCMASS(IPLEP,AMLEP)
         AMLEP=AMLEP*1.0E-3
         BMAX=0.
         BLAST=0.E+0
C     
C     --  SEARCH PEAK ROUGHLY
C     
         DO 10 I=0,50
            Q2=DQ2*FLOAT(I)
            EL=ENSTOP-Q2/2./AM
            PL=SQRT(EL**2-AMLEP**2)
            COST=(2.*ENSTOP*EL-AMLEP**2-Q2)/2./ENSTOP/PL
            IF(ABS(COST).GT.1.)GO TO 10
            DBLE0=ENSTOP
            DBLQ2=Q2
            AMP=DNELSQ2(DBLE0,IPAR,DBLQ2)
            IF (BMAX.LT.AMP) BMAX=AMP
            IF (AMP.LT.BLAST) GOTO 20
            BLAST=AMP
 10      CONTINUE
 20      CONTINUE
         ILAST=I-1
C     
C     --  SEARCH PEAK CLOSELY
C     
         Q2=DQ2*FLOAT(ILAST)+DQ2*0.02
         EL=ENSTOP-Q2/2./AM
         PL=SQRT(EL**2-AMLEP**2)
         COST=(2.*ENSTOP*EL-AMLEP**2-Q2)/2./ENSTOP/PL
         IF(ABS(COST).GT.1.) THEN
            AMP=0.0E+0
         ELSE
            DBLE0=ENSTOP
            DBLQ2=Q2
            AMP=DNELSQ2(DBLE0,IPAR,DBLQ2)
         ENDIF
         IF (AMP.GT.BLAST) THEN
C     
C     --     FORWARD
C     
            BMAX=AMP
            BLAST=AMP
            DO 30 I=2,50
               Q2=DQ2*FLOAT(ILAST)+DQ2*0.02*FLOAT(I)
               EL=ENSTOP-Q2/2./AM
               PL=SQRT(EL**2-AMLEP**2)
               COST=(2.*ENSTOP*EL-AMLEP**2-Q2)/2./ENSTOP/PL
               IF(ABS(COST).GT.1.)GO TO 30
               DBLE0=ENSTOP
               DBLQ2=Q2
               AMP=DNELSQ2(DBLE0,IPAR,DBLQ2)
               IF (BMAX.LT.AMP) BMAX=AMP
               IF (AMP.LT.BLAST) GOTO 50
               BLAST=AMP
 30         CONTINUE
         ELSE IF (ILAST.NE.0) THEN
C     
C     --     BACKWORD
C     
            DO 40 I=1,50
               Q2=DQ2*FLOAT(ILAST)-DQ2*0.02*FLOAT(I)
               EL=ENSTOP-Q2/2./AM
               PL=SQRT(EL**2-AMLEP**2)
               COST=(2.*ENSTOP*EL-AMLEP**2-Q2)/2./ENSTOP/PL
               IF(ABS(COST).GT.1.)GO TO 40
               DBLE0=ENSTOP
               DBLQ2=Q2
               AMP=DNELSQ2(DBLE0,IPAR,DBLQ2)
               IF (BMAX.LT.AMP) BMAX=AMP
               IF (AMP.LT.BLAST) GOTO 50
               BLAST=AMP
 40         CONTINUE
         ENDIF
 50      CONTINUE
C     
C     --  THROW DICE
C      
C      DO 100 LOOP=1,1000000
         Q2=RLU(DUM)*Q2MAXX
         DBLQ2=Q2
         HIT=RLU(DUM)*BMAX
         SIGMA=DNELSQ2(DBLE0,IPAR,DBLQ2)
         IF (SIGMA.GT.HIT) THEN
            EL=ENSTOP-Q2/2./AM
            PL=SQRT(EL**2-AMLEP**2)
            COST=(2.*ENSTOP*EL-AMLEP**2-Q2)/2./ENSTOP/PL
            PHI=6.283185*RLU(DUM)
            IF(ABS(COST).GT.1.)GO TO 100
            DO 220 KP=1,3
               IF(ABS(COST).LT.1.) CALL NECHADIR(COST,PHI,DNEU,DLEP)
               DO 210 JJ=1,3
                  PNP(JJ)=PNEU(JJ)-PL*DLEP(JJ)
 210           CONTINUE
               CALL MCVECBST(PNP,AM,EV,GM)
               PABS=SQRT(PNP(1)**2+PNP(2)**2+PNP(3)**2)
C               EABS=SQRT(PABS**2+AM**2) + (FNNUCL(PFABS)-FNNUCL(PABS))
               EABS=SQRT(PABS**2+AM**2) + (FNNPOT(1,PFABS)-FNNPOT(2,PABS))
               IF(EABS.LT.AM) THEN
                  PABS=0.001
               ELSE
                  PABS=SQRT(EABS**2-AM**2)
               END IF
               IF(KP.GE.3) GO TO 220
C               EL=ENSTOP-Q2/2./AM + (FNNUCL(PFABS)-FNNUCL(PABS))
               EL=ENSTOP-Q2/2./AM + (FNNPOT(1,PFABS)-FNNPOT(2,PABS))
               PL=0.001
               IF(EL.LT.AMLEP.AND.KP.GE.3 ) GO TO 100
               IF(EL.LT.AMLEP ) GO TO 220
               PL=SQRT(EL**2-AMLEP**2)
               COST=(2.*ENSTOP*EL-AMLEP**2-Q2)/2./ENSTOP/PL
               IF(ABS(COST).GT.1.) THEN
                  IF(KP.GE.3) GO TO 100
                  IF(KP.LT.3) GO TO 220
               END IF
  220       CONTINUE
            IF(NEPAUFLG.EQ.0 .AND. PABS.LT.PFERMI) GOTO 100
            CALL NECHADIR(COST,PHI,DNEU,DLEP)
            DO 330 JJ=1,3
               PLEP(JJ)=PL*DLEP(JJ)
 330        CONTINUE
            CALL MCVECBST(PLEP,AMLEP,EV,GM)
C     
C     --     NORMAL ENDING
C         
            RETURN
         ENDIF
  100 CONTINUE
      WRITE(6,*) ' IN NEELSLVC ( TOO MANY TRY ) '
C-Add 97/08/19 _ Y.H. for errortrap
      IERR=1
      RETURN
C     
C ++ ERROR RETURN
C
 9000 WRITE(6,900)
  900 FORMAT(' ***ERROR IN NEELSLVC(NOT ENOUGH NEUTRINO ENERGY)***')
      IERR=1
      RETURN
      END
