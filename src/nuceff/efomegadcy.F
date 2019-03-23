      SUBROUTINE EFOMEGADCY(POMEGA,MULDCY,IPBUF,PBUF)
*********************************************************
*      DECAY OMEGA MESON  -> PI+ PI- PI0 (88.8%)
*                            PI0 GAMMA   (8.5%)
*                            PI+ PI-     (2.2%)
*      LORENTZ BOOST IS CONSIDERED IN THIS PROGRAM
*     (Input)
*          POMEGA(3) : MOMENTUM OF OMEGA MESON
*
*     (Output)
*          MULDCY    : NUMBER OF DECAY PRODUCTS
*          IBUF(I)   : PARTICLE CODE OF I-TH PRODUCTS
*          PBUF(3,I) : MOMENTUM OF I-TH PRODUCTS
*
*     (Creation Date and Author)
*       2002.08.02    M.Ishitsuka
*********************************************************
      IMPLICIT NONE
C INPUT
      REAL POMEGA(3)
C OUTPUT
      INTEGER MULDCY
      REAL PBUF(3,10)
      INTEGER IPBUF(10)

      REAL DIROMEGA(3)
      REAL MASSOMEGA
      DATA MASSOMEGA /782./
      REAL PABS, BETA, GAM

      REAL RNBRANCH, BRANCH(3), SUMBRANCH
      DATA BRANCH /0.888,0.085,0.021/
      REAL P2BDABS
      REAL DIR2BD(3), P2BD1(3), P2BD2(3)
      REAL P3BD1(3), P3BD2(3), P3BD3(3)
      REAL MASSPI0, MASSPI, MASSGAMMA
      DATA MASSPI0 /135./
      DATA MASSPI  /140./
      DATA MASSGAMMA /0./

      REAL GAMFCT
      REAL RLU

      INTEGER I,K
      INTEGER IDUM
      
C  calculate the direction and gamma of omega meson
      PABS = SQRT(POMEGA(1)**2+POMEGA(2)**2+POMEGA(3)**2)
      BETA = SQRT(PABS**2/(MASSOMEGA**2+PABS**2))
      GAM = GAMFCT(BETA)
      
      IF (PABS.LE.0.) PABS=0.0001
      DO I=1,3
         DIROMEGA(I) = POMEGA(I)/PABS
      END DO

C determine the decay mode
      RNBRANCH = RLU(IDUM)
      SUMBRANCH = BRANCH(1) + BRANCH(2) + BRANCH(3) 
      IF(RNBRANCH .LE. BRANCH(1)/SUMBRANCH) THEN
CC pi+ , pi-, pi0, 3 body decay
         MULDCY = 3
         CALL DCY3BD(MASSOMEGA,MASSPI,MASSPI,MASSPI0,P3BD1,P3BD2,P3BD3)
         IPBUF(1) = 211
         IPBUF(2) = -211
         IPBUF(3) = 111
         CALL MCVECBST(P3BD1,MASSPI,DIROMEGA,GAM)
         CALL MCVECBST(P3BD2,MASSPI,DIROMEGA,GAM)
         CALL MCVECBST(P3BD3,MASSPI0,DIROMEGA,GAM)
         DO K=1,3
            PBUF(K,1)=P3BD1(K)
            PBUF(K,2)=P3BD2(K)
            PBUF(K,3)=P3BD3(K)
         ENDDO
      ELSE IF(RNBRANCH .LE. (BRANCH(1)+BRANCH(2))/SUMBRANCH) THEN
CC pi0 , gamma,    2 body decay
         MULDCY = 2
         CALL DCY2BD(MASSOMEGA,MASSPI0,MASSGAMMA,P2BDABS)
         CALL RNDIR(DIR2BD)
         IPBUF(1) = 111
         IPBUF(2) = 22
         DO K=1,3
            P2BD1(K)=P2BDABS*DIR2BD(K)
            P2BD2(K)=-P2BDABS*DIR2BD(K)
         ENDDO
         CALL MCVECBST(P2BD1,MASSPI0,DIROMEGA,GAM)
         CALL MCVECBST(P2BD2,MASSGAMMA,DIROMEGA,GAM)
         DO K=1,3
            PBUF(K,1)=P2BD1(K)
            PBUF(K,2)=P2BD2(K)
         ENDDO
      ELSE
CC pi+ , pi-,      2 body decay
         MULDCY = 2
         CALL DCY2BD(MASSOMEGA,MASSPI,MASSPI,P2BDABS)
         CALL RNDIR(DIR2BD)
         IPBUF(1) = 211
         IPBUF(2) = -211
         DO K=1,3
            P2BD1(K)=P2BDABS*DIR2BD(K)
            P2BD2(K)=-P2BDABS*DIR2BD(K)
         ENDDO
         CALL MCVECBST(P2BD1,MASSPI,DIROMEGA,GAM)
         CALL MCVECBST(P2BD2,MASSPI,DIROMEGA,GAM)
         DO K=1,3
            PBUF(K,1)=P2BD1(K)
            PBUF(K,2)=P2BD2(K)
         ENDDO
      ENDIF

      RETURN
      END
