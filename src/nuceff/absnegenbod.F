c********************************************************************

      SUBROUTINE ABSGENBOD(ierr)

C     SUBROUTINE TO GENERATE N-BODY EVENT
C   ACCORDING TO FERMI LORENTZ-INVARIANT PHASE SPACE
C   ADAPTED FROM FOWL (CERN W505) SEPT. 1974 BY F. JAMES
C   EVENTS ARE GENERATED IN THEIR OWN CENTER-OF-MASS,
C   BUT MAY BE TRANSFORMED TO ANY FRAME USING LOREN4
C
C   INPUT TO SUBROUTINE IS THRU COMMON BLOCK GENIN
C             NP=NUMBER OF OUTGOING PARTICLES (.LT. 19)
C             TECM=TOTAL ENERGY IN CENTER-OF-MASS
C             AMASS(I)=MASS OF ITH OUTGOING PARTICLE
C             KGENEV=1 FOR CONSTANT CROSS SECTION
C                      2 FOR FERMI ENERGY-DEPENDANCE
C
C   OUTPUT FROM SUBROUTINE IS THRU COMMON BLOCK GENOUT
C             PCM(1,I)=X-MOMENTUM IF ITH PARTICLE
C             PCM(2,I)=Y-MOMENTUM IF ITH PARTICLE
C             PCM(3,I)=Z-MOMENTUM IF ITH PARTICLE
C             PCM(4,I)=ENERGY OF ITH PARTICLE
C             PCM(5,I)=MOMENTUM OF ITH PARTICLE
C             WT=WEIGHT OF EVENT
c
c     Modified to use ranlux instead of grndm
c       2009.11.04 ; R. Tacik

      implicit none
      
      logical ierr

#include "absgeninout.h"

C      COMMON /NEGENIN / NP, TECM, AMASS(18), KGENEV
C      COMMON /NEGENOUT/ PCM(5,18) , WT
C      DIMENSION EMM(18)
C      DIMENSION RNO(50)
      REAL*4  EMM(18)
      REAL*4  RNO(50)

C--   PCM1 IS LINEAR EQUIV. OF PCM TO AVOID DOUBLE INDICES

C      DIMENSION EM(18),PD(18),EMS(18),SM(18),FFQ(18),PCM1(90)
      REAL*4   EM(18),PD(18),EMS(18),SM(18),FFQ(18),PCM1(90)

      INTEGER*4 NT
      EQUIVALENCE (NT,NP),(AMASS(1),EM(1)),(PCM1(1),PCM(1,1))

C     FFQ(N)=PI * (TWOPI)**(N-2) / (N-2)FACTORIAL

      DATA FFQ/0.,3.141592, 19.73921, 62.01255, 129.8788, 204.0131,
     2                       256.3704, 268.4705, 240.9780, 189.2637,
     3                       132.1308,  83.0202,  47.4210,  24.8295,
     4                        12.0006,   5.3858,   2.2560,   0.8859/

      INTEGER*4 KNT
      REAL*4    TWOPI
      
      DATA KNT,TWOPI/0,6.2831853073/

      integer*4 i,j,jk,ir,ndx
      integer*4 ntm1,ntm2,ntp1,ntnm4
      REAL*4    tm,tecmtm
      real*4    wtmaxq,wtmax
      real*4    emmax,emmin,esys,psave
      real*4    bang,cb,sb
      real*4    c,s
      real*4    beta, gama
      real*4    aa

      real*4    absnepdk
      external  absnepdk

C        INITIALIZATION

      ierr = .false.

      KNT = KNT + 1
      if (knt.eq.1) then
         WRITE (6,1160)
         WRITE (6,1200) NP,TECM,(AMASS(JK),JK=1,NP)
      end if
      IF (NT.LT.2)  GO TO 1001
      IF (NT.GT.18) GO TO 1002
      NTM1   = NT - 1
      NTM2   = NT - 2
      NTP1   = NT + 1
      NTNM4  = 3*NT - 4
      EMM(1) = EM(1)
      TM     = 0.0
      do I = 1, NT
         EMS(I) = EM(I)**2
         TM     = TM + EM(I)
         SM(I)  = TM
      end do
      
C        CONSTANTS DEPENDING ON TECM

      TECMTM = TECM - TM
      IF (TECMTM.LE.0.0) GO TO 1000
      EMM(NT) = TECM
      if (KGENEV.GT.1) then
         
C-- FERMI ENERGY DEPENDENCE FOR CROSS SECTION

         WTMAXQ = TECMTM**NTM2*FFQ(NT) / TECM
      else
      
C   CONSTANT CROSS SECTION AS FUNCTION OF TECM

         EMMAX = TECMTM + EM(1)
         EMMIN = 0.0
         WTMAX = 1.0
         do  I = 2, NT
            EMMIN  = EMMIN + EM(I-1)
            EMMAX  = EMMAX + EM(I)
            WTMAX  = WTMAX*ABSNEPDK(EMMAX,EMMIN,EM(I))
         end do
         WTMAXQ = 1.0 / WTMAX
      end if
      
C     CALCULATION OF WT BASED ON EFFECTIVE MASSES EMM
      
C--   FILL RNO WITH 3*NT-4 RANDOM NUMBERS,
C--   OF WHICH THE FIRST NT-2 ARE ORDERED.

      CALL ranlux(RNO,NTNM4)
      if (ntm2.lt.0) then
         return
      else if (ntm2.gt.0) then
         CALL ABSFLPSOR(RNO,NTM2)
         do J=2,NTM1
            EMM(J) = RNO(J-1)*(TECMTM) + SM(J)
         end do
      end if
      WT = WTMAXQ
      IR = NTM2
      do I = 1, NTM1
         PD(I) = ABSNEPDK(EMM(I+1),EMM(I),EM(I+1))
         WT    = WT*PD(I)
      end do
      
C--       COMPLETE SPECIFICATION OF EVENT (RAUBOLD-LYNCH METHOD)

      PCM(1,1) = 0.0
      PCM(2,1) = PD(1)
      PCM(3,1) = 0.0
      DO I = 2, NT
         PCM(1,I) = 0.0
         PCM(2,I) = -PD(I-1)
         PCM(3,I) = 0.0
         IR       = IR + 1
         BANG     = TWOPI*RNO(IR)
         CB       = COS(BANG)
         SB       = SIN(BANG)
         IR       = IR + 1
         C        = 2.0*RNO(IR) - 1.0
         S        = SQRT(1.0-C*C)
         IF (I.EQ.NT) GO TO 20
         ESYS     = SQRT(PD(I)**2+EMM(I)**2)
         BETA     = PD(I) / ESYS
         GAMA     = ESYS / EMM(I)
         do J = 1, I
            NDX = 5*J - 5
            AA  = PCM1(NDX+1)**2 + PCM1(NDX+2)**2 + PCM1(NDX+3)**2
            PCM1(NDX+5) = SQRT(AA)
            PCM1(NDX+4) = SQRT(AA+EMS(J))
            CALL ABSROTES2(C,S,CB,SB,PCM,J)
            PSAVE = GAMA*(PCM(2,J) + BETA*PCM(4,J))
            PCM(2,J) = PSAVE
         end do
         GO TO 30
20       do J = 1, I
            AA = PCM(1,J)**2 + PCM(2,J)**2 + PCM(3,J)**2
            PCM(5,J) = SQRT(AA)
            PCM(4,J) = SQRT(AA+EMS(J))
            CALL ABSROTES2(C,S,CB,SB,PCM,J)
         end do
30       CONTINUE
      end do

      RETURN
      
C          ERROR RETURNS
      
 1000 WRITE (6,1100)
      GO TO 1050
 1001 WRITE (6,1101)
      GO TO 1050
 1002 WRITE (6,1102)
 1050 WRITE (6,1150) KNT
      WRITE (6,1200) NP,TECM,(AMASS(JK),JK=1,NP)
      ierr = .true.
      return
 1100 FORMAT (28H0 AVAILABLE ENERGY NEGATIVE )
 1101 FORMAT (33H0 LESS THAN 2 OUTGOING PARTICLES )
 1102 FORMAT (34H0 MORE THAN 18 OUTGOING PARTICLES )
 1150 FORMAT (50H0 ABOVE ERROR DETECTED IN ABSGENBOD AT CALL NUMBER,I7)
 1160 FORMAT (37H0 FIRST CALL TO SUBROUTINE ABSGENBOD )
 1200 FORMAT (39H  INPUT DATA TO ABSGENBOD.         NP=   ,I6/
     +  ,8H   TECM=,E16.7,18H  PARTICLE MASSES=,5E15.5/(42X,5E15.5)
     +)
      END
