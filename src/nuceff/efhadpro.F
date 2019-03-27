*     -----------------------------------------------------
      SUBROUTINE EFHADPRO(PAIMOM,XPOS,NQPAI,NQN,NUCDIR,GAMMANUCL,IERROR)
*     -----------------------------------------------------
*
*     ( purpose )
C         MULTI-PION PRODUCTION IN PAI-N INTERACTION
C                    T.KAJTIA  27-JUN-88
*
*     ( input )
C           PAIMOM(3)  : MOMENTUM OF PION IN MEV/C
C           NQPAI  : CHARGE OF PION
C           NQN    : CHARGE OF NUCLEON
C           NUCDIR(3) : DIRECTION OF TARGET NUCLEON
C           GAMMANUCL : GAMMA OF TARGET NUCLEON (USED TO BOOST
C                       TO ORIGINAL EXP. SYSTEM)
*
*     ( output )
C           COMMON/PIINTH/
C
C     ( Updating record )
C                    M.Ishitsuka 16-MAR-05
C                        Distinguish elastic and inelastic scattering
C                    M.Ishitsuka 16-MAR-05
C                        Consider Fermi motion of target nucleon
      implicit none

      REAL PAIMOM(3),XPOS(3),PECM(3)
      REAL NUCDIR(3), GAMMANUCL
      REAL BUF(10),UV(3),PBUF(3,10)
      INTEGER IPBUF(10)
      REAL AM(9)
C      COMMON /MASS/AM(9)

#include "efpiinth.h"
C      INTEGER NPARTL,IPINT,IPIORI,LOOPPI
C      REAL PINT,XINT
C      COMMON /PIINTH/ NPARTL,IPINT(50),PINT(3,50),XINT(3,50),IPIORI(50),
C     &       LOOPPI

      REAL XMNUCL,XMPAI
      DATA XMNUCL/939. /
      DATA XMPAI/139. /
C                      MEV/C**2
      REAL XMULTP
      INTEGER MULTIP
      REAL RLU,GM,BM,EPAI,S,SQRTS,PT,PT2,FPTW,RND,RNNUM
      INTEGER NQTOT,ICHAD,IPP,NECHARGE,NL,NH,IP0
      INTEGER NQPAI,NQN,IERROR,IERTOT,IDUM,I
      REAL PABS,FNPKNO,RAN
C
      IERTOT=0
      IERROR=0
      PABS=SQRT(PAIMOM(1)**2+PAIMOM(2)**2+PAIMOM(3)**2)
      CALL EFPNMLTP(PABS,NQPAI,NQN,XMULTP,IERROR)
      IF(IERROR.EQ.1) RETURN
C
C
C -- DETERMINE MULTIPLICITY
C
   99 CONTINUE
      BUF(1)=FNPKNO(1,XMULTP)
      DO 100 I=2,10
         BUF(I)=BUF(I-1)+FNPKNO(I,XMULTP)
  100 CONTINUE
C  110 RAN=RNDM(DUM)*BUF(10)
  110 RAN=RLU(IDUM)*BUF(10)
      I=1
  150 IF(BUF(I).GE.RAN)GO TO 200
      I=I+1
      IF(I.EQ.11)GO TO 110
      GO TO 150
  200 MULTIP=I
      IF(MULTIP.LT.2) MULTIP=2
      IF(MULTIP.GT.7) MULTIP=7
      IF(MULTIP.EQ.2 .AND. NQPAI.EQ.1 .AND. NQN.EQ.1) GOTO 99
      IF(MULTIP.EQ.2 .AND. NQPAI.EQ.-1 .AND. NQN.EQ.0) GOTO 99
C
C
C -- DETERMINE MOMENTUM OF EACH PARTICLE
C
  311 CONTINUE
      AM(1)=XMNUCL
      DO 300 IPP=2,9
  300 AM(IPP)=XMPAI
      EPAI= SQRT ( PABS**2 + XMPAI**2 )
C      S=2.*XMNUCL*EPAI
C invariant mass (M.Ishitsuka)
      S=(EPAI+XMNUCL)**2 - PABS**2
      SQRTS=SQRT(S)
      CALL VCPHSVCT(SQRTS,AM,MULTIP,PBUF,IERROR)
      IF(IERROR.EQ.1 ) IERTOT=IERTOT+1
      IF(IERROR.EQ.1 .AND.IERTOT.GE.3) RETURN
      IF(IERROR.EQ.1 .AND.IERTOT.LT.3) GO TO 99
C
C -- CONSTRAINT ON THE PT OF PROTON
C -- CONSTRAINT THAT THE PROTON GO BACKWARD IN C.M.S.
      UV(1)=PAIMOM(1)/PABS
      UV(2)=PAIMOM(2)/PABS
      UV(3)=PAIMOM(3)/PABS
      PT = PBUF(1,1)*UV(1) + PBUF(2,1)*UV(2) + PBUF(3,1)*UV(3)
      PT2=(0.001*PT)**2
C                    GEV**2
      FPTW=EXP(-PT2/0.1)
C      RND=RNDM(DUM)
      RND=RLU(IDUM)
      IF(RND.GT.FPTW) GO TO 311
      IF(PT.GT.0.) GO TO 311
C
C
C -- LORENTZ BOOST
C
      GM=SQRT(S+PABS**2)/SQRTS
      DO 400 IPP=1,MULTIP
         BM=AM(IPP)
         PECM(1)=PBUF(1,IPP)
         PECM(2)=PBUF(2,IPP)
         PECM(3)=PBUF(3,IPP)
         CALL MCVECBST ( PECM, BM, UV, GM)
C BOOST TO ORIGINAL EXP. SYSTEM
         CALL MCVECBST ( PECM, BM, NUCDIR, GAMMANUCL)
         PBUF(1,IPP)=PECM(1)
         PBUF(2,IPP)=PECM(2)
         PBUF(3,IPP)=PECM(3)
  400 CONTINUE
C
C
C -- DETERMINE PARTICLE (CHARGE)
C
      NQTOT= NQPAI+NQN
  430 ICHAD=0
      DO 410 IPP=1,MULTIP
         IF(AM(IPP).EQ.XMNUCL) THEN
C            IF(RNDM(DUM).GT.0.5) THEN
            IF(RLU(IDUM).GT.0.5) THEN
               IPBUF(IPP)=2212
            ELSE
               IPBUF(IPP)=2112
            END IF
         ELSE
C            RNNUM=RNDM(DUM)
            RNNUM=RLU(IDUM)
            IF(RNNUM.GT.0.667) THEN
               IPBUF(IPP)=211
            ELSE IF(RNNUM.GT.0.333) THEN
               IPBUF(IPP)=111
            ELSE
               IPBUF(IPP)=-211
            END IF
         END IF
         ICHAD=ICHAD+NECHARGE(IPBUF(IPP))
  410 CONTINUE
      IF(ICHAD.NE.NQTOT) GO TO 430
C REJECT ELASTIC SCATTERING
      IF(MULTIP.EQ.2 .AND. NQPAI.EQ.NECHARGE(IPBUF(2))) GOTO 430
C
C
C -- FILLING
C
      IF(NPARTL.GE.50)NPARTL=49
      NL=NPARTL+1
      NH=NPARTL+MULTIP
      IF(NH.GT.50) THEN
         WRITE(6,600)
  600 FORMAT('  ********** EFHADPRO WARNING: MULTIPLICITY > 50 IN HAD-INT ******')
         NH=50
      END IF
C
      DO 500 IPP=NL,NH
         IP0=IPP-NPARTL
         IPIORI(IPP)=LOOPPI
         PINT(1,IPP)=PBUF(1,IP0)
         PINT(2,IPP)=PBUF(2,IP0)
         PINT(3,IPP)=PBUF(3,IP0)
         XINT(1,IPP)=XPOS(1)
         XINT(2,IPP)=XPOS(2)
         XINT(3,IPP)=XPOS(3)
         IPINT(IPP)=IPBUF(IP0)
  500 CONTINUE
      NPARTL=NH

      IERROR=0
C
C
      RETURN
      END
