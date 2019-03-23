**********************************************************************
*     -------------------------------------
      SUBROUTINE VCRAMBO(N,ET,XM,P,WT,IERR)
*     -------------------------------------
*
*     ( purpose )
*      A DEMOCRATIC MULTI-PARTICLE PHASE SPACE GENERATOR
*
*     ( input )
*        N      : NUMBER OR PARTICLES (>1, IN THIS VERSION <101)
*        ET     : TOTAL CENTORE-OF-MASS ENERGY
*        XM     : PARTICLE MASSES ( DIM=100 )
*
*     ( output )
*        P      : MOMENTUM
*        WTMAX  : WEIGHT OF THE EVENT
*        IERR   : ERROR CODE
*
*     ( creation date and author )
*       ????. ?.?? ; First version by S.D.Ellis, R.Kleiss W.J.Stirling
*       1989. 7.24 ; mod for KAMIOKANDE by K.Kaneyuki
*       1992. 8. 6 ; mod for S.K.       by K.Kaneyuki
*                    RAMBO -> VCRMBO      
*       1995. 2. 7 ; change name
*                    VCRMBO -> VCRAMBO      
*       2001. 2. 8 ; RNRMBO -> RM48
*
***********************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  XM(100),P(4,100),Q(4,100),Z(100),R(4),
     &     B(3),P2(100),XM2(100),E(100),V(100),IWARN(5),IER(2)
      DATA ACC   /1.D-14/
      DATA ITMAX /6/
      DATA IBEGIN/0/
      DATA IWARN /5*0/
      DATA IER   /2*0/
C     
C     -- INITIALIZATION STEP: FACTORIALS FOR THE PHASE SPACE WEIGHT
C     
      IERR=0
      IF (IBEGIN.NE.0) GOTO 103
      IBEGIN=1
      TWOPI=8.*DATAN(1.D0)
      PO2LOG=DLOG(TWOPI/4.)
      Z(2)=PO2LOG
      DO 101 K=3,100
         Z(K)=Z(K-1)+PO2LOG-2.*DLOG(DFLOAT(K-2))
 101  CONTINUE
      DO 102 K=3,100
         Z(K)=(Z(K)-DLOG(DFLOAT(K-1)))
 102  CONTINUE
C     
C     -- CHECK ON THE NUMBER OF PARTICLES
C     
 103  IF (N.GT.1 .AND. N.LT.101) GOTO 104
      PRINT 1001,N
      IF (IER(1).LE.5) THEN
         IER(1)=IER(1)+1
         IERR=1
         RETURN
      ELSE
         STOP
      ENDIF
C     
C     -- CHECK WHETHER TOTAL ENERGY IS SUFFICIENT; COUNT NONZERO MASSES
C     
 104  XMT=0.
      NM=0
      DO 105 I=1,N
         IF (XM(I).NE.0.D0) NM=NM+1
         XMT=XMT+DABS(XM(I))
 105  CONTINUE
      IF (XMT.LE.ET) GOTO 201
      PRINT 1002,XMT,ET
      PRINT *,(XM(I),I=1,N)
      IF (IER(2).LE.5) THEN
         IER(2)=IER(2)+1
         IERR=1
         RETURN
      ELSE
         STOP
      ENDIF
C     
C     -- THE PARAMETER VALUES ARE NOW ACCEPTED
C     
C     GENERATE N MASSLESS MOMENTA IN INFINITE PHASE SPACE
C     
 201  DO 202 I=1,N
C         CALL RNRMBO(XRNB)
         CALL RM48(XRNB,1)
         C=2.*XRNB-1.
         S=DSQRT(1.-C*C)
C         CALL RNRMBO(XRNB)
         CALL RM48(XRNB,1)
         F=TWOPI*XRNB
C         CALL RNRMBO(XRNB1)
C         CALL RNRMBO(XRNB2)
         CALL RM48(XRNB1,1)
         CALL RM48(XRNB2,1)
         Q(4,I)=-LOG(XRNB1*XRNB2)
         Q(3,I)=Q(4,I)*C
         Q(2,I)=Q(4,I)*S*DCOS(F)
         Q(1,I)=Q(4,I)*S*DSIN(F)
 202  CONTINUE
C     
C     -- CALCULATE THE PARAMETERS OF THE CONFORMAL TRANSFORMATION
C     
      DO 203 I=1,4
         R(I)=0.
 203  CONTINUE
      DO 204 I=1,N
         DO 204 K=1,4
            R(K)=R(K)+Q(K,I)
 204     CONTINUE
         RMAS=DSQRT(R(4)**2-R(3)**2-R(2)**2-R(1)**2)
         DO 205 K=1,3
            B(K)=-R(K)/RMAS
 205     CONTINUE
         G=R(4)/RMAS
         A=1./(1.+G)
         X=ET/RMAS
C     
C     -- TRANCEFORM THE Q'S CONFORMALLY INTO THE P'S
C     
         DO 207 I=1,N
            BQ=B(1)*Q(1,I)+B(2)*Q(2,I)+B(3)*Q(3,I)
            DO 206 K=1,3
               P(K,I)=X*(Q(K,I)+B(K)*(Q(4,I)+A*BQ))
 206        CONTINUE
            P(4,I)=X*(G*Q(4,I)+BQ)
 207     CONTINUE
C     
C     -- CALCULATE WEIGHT AND POSSIBLE WARNINGS
C     
         WT=PO2LOG
         IF (N.NE.2) WT=(2.*N-4)*DLOG(ET)+Z(N)
         IF (WT.GE.-180.D0) GOTO 208
         IF (IWARN(1).LE.5) PRINT 1004,WT
         IWARN(1)=IWARN(1)+1
 208     IF (WT.LE.174.D0) GOTO 209
         IF (IWARN(2).LE.5) PRINT 1005,WT
         IWARN(2)=IWARN(2)+1
C     
C     -- RETURN FOR WEIGHTED MASSLESS MOMENTA
C     
 209     IF (NM.NE.0) GOTO 210
         WT=DEXP(WT)
         RETURN
C     
C     -- MASSIVE PARTICLES : RESCALE THE MOMENTA BY A FACTOR X
C     
 210     XMAX=DSQRT(1.-(XMT/ET)**2)
         DO 301 I=1,N
            XM2(I)=XM(I)**2
            P2(I)=P(4,I)**2
 301     CONTINUE
         ITER=0
         X=XMAX
         ACCU=ET*ACC
 302     F0=-ET
         G0=0.
         X2=X*X
         DO 303 I=1,N
            E(I)=DSQRT(XM2(I)+X2*P2(I))
            F0=F0+E(I)
            G0=G0+P2(I)/E(I)
 303     CONTINUE
         IF (DABS(F0).LE.ACCU) GOTO 305
         ITER=ITER+1
         IF (ITER.LE.ITMAX) GOTO 304
         PRINT 1006,ITMAX
         GOTO 305
 304     X=X-F0/(X*G0)
         GOTO 302
 305     DO 307 I=1,N
            V(I)=X*P(4,I)
            DO 306 K=1,3
               P(K,I)=X*P(K,I)
 306        CONTINUE
            P(4,I)=E(I)
 307     CONTINUE
C     
C     -- CALCULATE THE MASS-EFFECT WEIGHT FACTOR
C     
         WT2=1.
         WT3=0.
         DO 308 I=1,N
            WT2=WT2*V(I)/E(I)
            WT3=WT3+V(I)**2/E(I)
 308     CONTINUE
         WTM=(2.*N-3.)*DLOG(X)+DLOG(WT2/WT3*ET)
C     
C     -- RETURN FOR WEIGHTED MASSIVE MOMENTA
C     
         WT=WT+WTM
         IF (WT.GE.-180.D0) GOTO 309
         IF (IWARN(3).LE.5) PRINT 1004,WT
         IWARN(3)=IWARN(3)+1
 309     IF (WT.LE.174.D0) GOTO 310
         IF (IWARN(4).LE.5) PRINT 1005,WT
         IWARN(4)=IWARN(4)+1
 310     WT=DEXP(WT)
         RETURN
         
 1001    FORMAT(' RAMBO FAILS : # OF PARTICLES =',I5,' IS NOT ALLOWED')
 1002    FORMAT(' RAMBO FAILS : TOTAL MASS =',D15.6,' IS NOT',
     &        ' SMALLER THAN TOTAL ENERGY =',D15.6)
 1004    FORMAT(' RAMBO WARNS : WEIGHT = EXP(',F20.9,') MAY UNDERFLOW')
 1005    FORMAT(' RAMBO WARNS : WEIGHT = EXP(',F20.9,') MAY OVERFLOW')
 1006    FORMAT(' RAMBO WARNS : ',I3,' ITERATIONS DIS NOT GIVE THE',
     &        'DESIRED ACCURACY =',D15.6)
         END
