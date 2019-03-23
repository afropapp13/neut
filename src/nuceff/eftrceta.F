*      --------------------------------------
       SUBROUTINE EFTRCETA(POSI,PIN,NSCAT,NDIS,NABS,NABS2,PIN1,
     &                     PIN2,KIN1,KIN2)
*      --------------------------------------
*     (purpose) calculate nuclear effect for ETA
*
*     input   : POSI(3)  initial position of ETA
*             : PIN(3)  initial momentum of ETA
*     output  : NSCAT    number of scattering in nuclear
*               NDISA   flag disappear (change to pi): 1
*                                  not               : 0
*               NABS    Pion absorption              : 1
*                                  not               : 0
*               PIN1(3)  final momentum of meson 1
*               PIN2(3)  final momentum of meson 2
*               KIN1,KIN2   kind of mesons 221:  eta
*                                          211:  pi+ 
*                                         -211:  pi-
*                                          111:  pi0
*                                         
*     (Creation Date and Author)
*     2007.11.05 ; G.Mitsuka - add PFCUT for pauli blocking
*     2010.05    ; P.de Perio - Modified to work with nucleon ejection
*                             - EFTRACE arguments modified
********************************************************************
C             DIRRS-- DIRECTION OF RESONANCE MOTION IN LAB.SYS.
C             PNCMS-- MOMENTUM VECTOR OF RESONANCE
C             POSM-- POSITION OF MESON IN THE NUCLEUS IN UNIT OF FM
C             DIR1-- DIRECTION OF MESON, UNIT VECT.
C              PIN0(3)-- MOM. OF P OR N,      VIN0(3)-- VEL. OF P OR N.
C              PIN1(3)-- MOM. OF MESON,           "
C              PIN2(3)-- MOM. OF NEUTRINO         "
C           * VREL -- RELATIVE VELOCITY
#include "vcwork.h"
#include "vcvrtx.h"
#include "neutparams.h"
       DIMENSION PN(3), PM1(3), PM2(3), DIRTMP(3)
       DIMENSION DIRRS(3)
       DIMENSION POSM(3), DIR1(3)
       DIMENSION VIN0(3), VIN1(3), VIN2(3), VREL(3)
       DIMENSION PIN0(3), PIN1(3), PIN2(3)
       real POSI(3),PIN(3),PFI(3)
       INTEGER L,NSCAT,NDIS,NSTEP,dummy1,KFI,nabs,nabs2
       real apn,ap
       REAL PFCUT

C     APN and AP are not used now in eftrace.F(i.e. dummy number)
       DATA APN/16./
       DATA AP/8./
C     PFCUT is set by the scaling to 16O(225MeV/c)
       PFCUT = 260.*(PFSURF*1000./225.) ! PFSURF is GeV unit
       
      STEP=0.07
C    STEP OF CALCULATION (FM)
      PRMAS= 940.
      AM1  = 550.
      AM2  =   0.
      KIN1 = 221
      KIN2 = 0
C    PRMAS - MASS OF NUCLEON, AM1 - MASS OF MESON, AM2 - MASS OF LEPTON

CC   initialization
      do 99 L = 1,3
       POSM(L)= POSI(L)
       PIN1(L)= PIN(L)
       PIN2(L)= 0.0
  99  continue
       NSCAT = 0
       NDIS  = 0
       NSTEP = 0

C + INITIAL  POSITION CHECK
      ABSPOS = SQRT ( POSM(1)**2 + POSM(2)**2 + POSM(3)**2 )
      IF ( ABSPOS.GT.3.) GO TO 1000
C
      P1ABS=SQRT(PIN1(1)**2+PIN1(2)**2+PIN1(3)**2)
      V1ABS=SQRT(P1ABS**2/(P1ABS**2+AM1**2))
      DO 10 I=1,3
      DIR1(I)=PIN1(I)/P1ABS
      VIN1(I)=DIR1(I)*V1ABS
   10 CONTINUE
      
   30 CONTINUE
CC   eta loop
      nstep = nstep + 1
      if (nstep.gt.10000) then 
      print *,'TOO MANY STEPS'
      return
      endif

C +  FERMI MOMENTUM OF TARGET NUCLEON
      call RNFERM(PFERM,dummy1)
      CALL RNDIR(DIRTMP)
      VPABS=SQRT(PFERM**2/(PFERM**2+PRMAS**2))
      DO 41 I=1,3
      PIN0(I) = DIRTMP(I)*PFERM
      VIN0(I) = DIRTMP(I)*VPABS
   41 CONTINUE

C +CHECK THE RELATIVE VELOCITY AND IF THE MESON AND THE NUCLEON INTERACT
C   OR NOT, SINCE THE INTERACTION RATE IS PROPOTIONAL TO V(REL).

      VRELMX=V1ABS + 1.0
c            300 - ABSOLUTE MAX. OF FERMI MOMENTUM.
      VREL(1)=VIN1(1)-VIN0(1)
      VREL(2)=VIN1(2)-VIN0(2)
      VREL(3)=VIN1(3)-VIN0(3)
      VRLABS=SQRT( VREL(1)**2 + VREL(2)**2 + VREL(3)**2 )
C      IF(VRLABS/VRELMX .LT. RNDM(DUM)) GO TO 30
      IF(VRLABS/VRELMX .LT. RLU(IDUM)) GO TO 30

C
C + GET MESON-NUCLEON INV. MASS AND MOM. IN C.M.S.
      CALL AINVMS(PRMAS,AM1,PIN0,PIN1,RMASS,PCMS)
C              RMASS-- INVARIANT MASS
C              PCMS -- ABS. MOMENTUM IN CENTER OF MASS SYSTEM

C
C + CROSS SECTION OF ETA-NUCLEON
      Call SIGETA(RMASS,PCMS,SIGM)
C              SIGM-- CROSS SECTION (MB)
C
C + INTERACTION IN 0.07 FM ?
      PRINT=1.-EXP(-STEP/PASSL(0.138,SIGM))
C      IF(PRINT.LT.RNDM(DUM)) GO TO 300
      IF(PRINT.LT.RLU(IDUM)) GO TO 300

C (1)INTERACTION

             NSCAT = NSCAT + 1
555     continue
         CALL TYPDCY(RMASS,AMN,AMM1,AMM2,PN,PM1,PM2,NDIS,KIN1,KIN2)

C         AMN : NUCLEON MASS, AMM1 : MESON MASS(1), AMM2 : MESON MASS(2)
C         PN  :   "     MOM., PM1  :   "   MOM.(1), PM2  :   "   MOM.(2)
C                    -- C.M.S. OF RMASS SYSTEM
C         KIN1,KIN2 :  kind of particle
C
         PTOTX=PIN1(1)+PIN0(1)
         PTOTY=PIN1(2)+PIN0(2)
         PTOTZ=PIN1(3)+PIN0(3)
         PTOTRS=SQRT( PTOTX**2 + PTOTY**2 + PTOTZ**2 )
         DIRRS(1)=PTOTX/PTOTRS
         DIRRS(2)=PTOTY/PTOTRS
         DIRRS(3)=PTOTZ/PTOTRS
         BETARS = SQRT ( PTOTRS**2 / ( RMASS**2 + PTOTRS**2 ) )
C                  BETA OF RESONANCE
         GAM=GAMFCT(BETARS)
C                  GAMMA FACTOR

C       + BOOST OF PN & PM1 TO LAB.SYS.
         CALL MCVECBST(PN,PRMAS,DIRRS,GAM)
         PNABS= SQRT( PN(1)**2 + PN(2)**2 + PN(3)**2 )
C       + MOM. OF NUCLEON IS GREATER THAN 
C         THE FERMI ENERGY (260MEV/C FOR 16O)?
         If(PNABS.LT.PFCUT) then
         print *,'PAULI BLOCKING'
         GO TO 555
         endif

         CALL MCVECBST(PM1,AMM1,DIRRS,GAM)
         PIN1(1) = PM1(1)
         PIN1(2) = PM1(2)
         PIN1(3) = PM1(3)
         if (NDIS.ge.2) then
          CALL MCVECBST(PM2,AMM2,DIRRS,GAM)
          PIN2(1) = PM2(1)
          PIN2(2) = PM2(2)
          PIN2(3) = PM2(3)
         endif
c         write(98,*) '***************************'
c         write(98,*) NDIS,nscat,KIN1,KIN2
c         write(98,*) AMN,AMM1,AMM2
c         write(98,*) pm1(1),pm1(2),pm1(3)
c         write(98,*) pm2(1),pm2(2),pm2(3)


C         CHECK WHETHER ETA DISAPPEAR OR NOT    

C           IF DISAPPEAR -> CALCULATE NUC EFF. OF PI
C          What is happening here? (-Patrick)
          if (NDIS.ge.1) then
            print *,'ETA changes to pion IN O16'
            CALL EFTRACE(APN,AP,KIN1,POSM,PM1,NABS,NEL,NINEL,KFI,PIN1) 
            if (NABS.ne.0) then
               KIN1 = KFI
            endif
            if (NDIS.ge.2) then
               call EFTRACE(APN,AP,KIN2,POSM,PM2,NABS2,NEL,NINEL,KFI,PIN2) 
               if (NABS2.ne.0) then
                  KIN2 = KFI
               endif                  
            endif
            return
          endif


C         GO TO 1000

C
C (2)NO INT. --> GO TO NEXT STEP

 300     CONTINUE

 50      POSM(1)=POSM(1)+DIR1(1)*0.07
         POSM(2)=POSM(2)+DIR1(2)*0.07
         POSM(3)=POSM(3)+DIR1(3)*0.07
         ABSPOS = SQRT ( POSM(1)**2 + POSM(2)**2 + POSM(3)**2 )
         IF ( ABSPOS.GT.3.) GO TO 1000
         GO TO 30

C      escape 
 1000 continue
C
C
      return

      END


