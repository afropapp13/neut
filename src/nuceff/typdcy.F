C   19/10/85 603241822  MEMBER NAME  TYPDCY   (FORT)        FORTRAN77
C
      SUBROUTINE TYPDCY(RMASS,AMN,AMM1,AMM2,PN,PM1,PM2,NDIS,KIN1,KIN2)
C
C   DETERMINE THE FOLLOWING KINEMATICS DECAYED FROM N(1440,1535,1710)
C                --> modified to consider only N(1535) 
C
C         AMN : NUCLEON MASS, AMM1 : MESON MASS(1), AMM2 : MESON MASS(2)
C         PN  :   "     MOM., PM1  :   "   MOM.(1), PM2  :   "   MOM.(2)
C                    -- C.M.S. OF RMASS SYSTEM
C              KIN1,KIM2   kind of particle
Cccccc         NDIS : DISAPPER OR NOT
C
C  INPUT : RMASS - RESONANCE MASS
C
C      T.KAJITA   19-OCT-85
C      J.KAMEDA   1997.11.
C
      DIMENSION DIR(3), DRTEMP(3)
      DIMENSION PN(3), PM1(3), PM2(3)
      DIMENSION BR1(10),BR2(10),BR3(10),BR4(10)
      INTEGER NDIS

C
       DATA BR1 / 0.35 ,0.42 , 0.1  ,0.2  , 0., 0., 0., 0., 0., 0./
       DATA BR2 / 0.333 ,0.333 , 1.0 ,0.0  , 0., 0., 0., 0., 0., 0./
       DATA BR3 / 1.   ,6.   , 8.  ,9.  , 0., 0., 0., 0., 0., 0./
       DATA BR4 / 0.33 ,0.67 , 1.0  ,0. , 0., 0., 0., 0., 0., 0./
       NDIS = 0

C BRANCHING RATIO  N-ETA, N-PAI, DELTA-PAI, N-ROU
C
C
C   N(1535)
C  110 RAN=RNDM(DUM)
  110 RAN=RLU(IDUM)
      IF(  BR1(1).GT.RAN)                GO TO 500
      IF( (BR1(1)+BR1(2)).GT.RAN)        GO TO 1000
      IF( (BR1(1)+BR1(2)+BR1(3)).GT.RAN) GO TO 2000
                                         GO TO 3000

C
C
C I) N+ETA
 500  AMN=939.
      AMM1=550.
      AMM2=140.
      KIN1 = 221
      kIN2 = 0
      NDIS = 0      
      CALL DCY2BD(RMASS, AMN, AMM1, P1)
      CALL RNDIR(DIR)
      DO 510 K=1,3
      PN(K)=P1*DIR(K)
      PM1(K)=-P1*DIR(K)
      PM2(K)=0.
 510  CONTINUE

      RETURN


C
C
C I) N+PAI
 1000 AMN=939.
      AMM1=140.
      AMM2=140.
      NDIS = 1      
      KIN2 = 0
      CALL DCY2BD(RMASS, AMN, AMM1, P1)
      CALL RNDIR(DIR)
      DO 1010 K=1,3
        PN(K)=P1*DIR(K)
        PM1(K)=-P1*DIR(K)
        PM2(K)=0.
 1010 CONTINUE

C      rnkind = rndm(dum)
      rnkind = RLU(idum)
      IF( (BR2(1)+BR2(2)).LE.RNKIND)  KIN1 = -211
      IF( (BR2(1)+BR2(2)).GT.RNKIND)  KIN1 = 211
      IF(  BR2(1).GT.RNKIND)          KIN1 = 111
      RETURN
C
C
C II) DELTA+PAI
 2000 AMN=939.
      AMM1=140.
      AMM2=140.
      NDIS = 2      
 2020 DEMASS=RSMASS(RMASS,1232.,AMM2,115.)
      IF ( DEMASS.LE.(AMN+AMM1+1.) ) GO TO 2020
      IF ( DEMASS.GE.(RMASS-AMM1-1.) ) GO TO 2020
      CALL DCY2BD(RMASS, AMM1, DEMASS, P2)
      CALL RNDIR(DIR)
      DO 2010 K=1,3
      PM2(K)=-P2*DIR(K)
      DRTEMP(K)=DIR(K)
 2010 CONTINUE
      BETA = SQRT ( P2**2/(P2**2+DEMASS**2) )
      GAM=GAMFCT(BETA)
C
C  DECAY OF DELTA

      CALL DCY2BD(DEMASS, AMN, AMM1, P22)
      CALL RNDIR(DIR)
      PN(1)=DIR(1)*P22
      PN(2)=DIR(2)*P22
      PN(3)=DIR(3)*P22
C
      PM1(1)=-DIR(1)*P22
      PM1(2)=-DIR(2)*P22
      PM1(3)=-DIR(3)*P22

      CALL MCVECBST(PM1, AMM1, DRTEMP, GAM)
      CALL MCVECBST(PN, AMN, DRTEMP, GAM)

cccc particle ID      

C      rnkind = rndm(dum)*9.0
      rnkind = RLU(idum)*9.0
      IF(  BR3(3).LE.RNKIND) then
       KIN1 = 111
       KIN2 = -211
      endif
      IF(  BR3(3).GT.RNKIND) then
       KIN1 = 111
       KIN2 = 111
      endif
      IF(  BR3(2).GT.RNKIND) then
       KIN1 = 211
       KIN2 = -211
      endif
      IF(  BR3(1).GT.RNKIND) then
       KIN1 = 211
       KIN2 = 111
      endif

      RETURN

C
C III) N+RHO
 3000 AMN=939.
      AMM1=140.
      AMM2=140.
      NDIS = 3      
 3020 RHMASS=RSMASS(RMASS,769.,AMN,154.)
      IF ( RHMASS.LT.(AMM1+AMM2+1.) ) GO TO 3020
      IF ( RHMASS.GT.(RMASS-AMN-1.) ) GO TO 3020
      CALL DCY2BD(RMASS, AMN, RHMASS, P3)
      CALL RNDIR(DIR)
      DO 3010 K=1,3
      PN(K)=-P3*DIR(K)
      DRTEMP(K)=DIR(K)
 3010 CONTINUE
      BETA = SQRT ( P3**2/(P3**2+RHMASS**2) )
      GAM=GAMFCT(BETA)
C
C DECAY OF RHO
      CALL DCY2BD(RHMASS, 140., 140., P32)
      CALL RNDIR(DIR)
      PM1(1)=DIR(1)*P32
      PM1(2)=DIR(2)*P32
      PM1(3)=DIR(3)*P32
C
      PM2(1)=-DIR(1)*P32
      PM2(2)=-DIR(2)*P32
      PM2(3)=-DIR(3)*P32
      CALL MCVECBST(PM1, AMM1, DRTEMP, GAM)
      CALL MCVECBST(PM2, AMM2, DRTEMP, GAM)
C
ccc particle ID
C      rnkind = rndm(dum)
      rnkind = rlu(idum)
      IF(  BR4(2).LE.RNKIND) then
       KIN1 = 111
       KIN2 = -211
      endif
      IF(  BR4(2).GT.RNKIND) then
       KIN1 = 211
       KIN2 = -211
      endif
      IF(  BR4(1).GT.RNKIND) then
       KIN1 = 211
       KIN2 = 111
      endif


      RETURN
      END

