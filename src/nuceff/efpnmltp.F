*     -----------------------------------------------------
      SUBROUTINE EFPNMLTP(PAIMOM,NQPAI,NQN,XMULTP,IERROR)
*     -----------------------------------------------------
*
*     ( purpose )
C        MEAN XMULTPLICITY OF (PAI+NUCLEON) IN
C        PI-N INTERACTION AT HIGH ENERGIES (>800MEV).
C                    T.KAJITA  27-JUN-88
*
*     ( input )
C           PAIMOM : MOMENTUM OF PION IN MEV/C ( SHOULD BE > 800MEV)
C           NQPAI  : CHARGE OF PION     <=== NO MEANING AT PRESENT
C           NQN    : CHARGE OF NUCLEON  <=== NO MEANING AT PRESENT
C
*     ( output )
C           XMULTP : MEAN MULTIPLICITY OF (PAI+N)
C           IERROR : =0 NO ERROR   =1 ERROR IN (NQPAI+NQN) OR PAIMOM
*
C     DATA ARE FROM J.WHITMORE PHYS.REP.27C(1976)188.
C
C
      DATA XMNUCL/0.939/
      DATA XMPAI/0.139/
C
C -- STATEMENT FUNCTION (CHARGED XMULTPLICITY FOR PAI(+)P AND PAI(-)P)
      FMUL2(ECM2) = 1.94 + 0.28*ALOG(ECM2) + 0.130*(ALOG(ECM2))**2
      FMUL0(ECM2) = 0.59 + 0.81*ALOG(ECM2) + 0.074*(ALOG(ECM2))**2
C
C
      XMULTP=0.
      IERROR=0
      IF(ABS(NQPAI+NQN).GT.2) THEN
         IERROR=1
         RETURN
      END IF

*     IF(PAIMOM.LT.500.) THEN
      IF(PAIMOM.LT.340.) THEN
         IERROR=1
         RETURN
      END IF

C
C
C CALCULATE S (=ECM**2)
      EPAI= SQRT ( (0.001*PAIMOM)**2 + XMPAI**2 )
C      S=2.*XMNUCL*EPAI
C invariant mass
      S=(epai+xmnucl)**2 - (0.001*PAIMOM)**2
C
C CALCULATE XMULTPLICITY
C ++ TOTAL CHARGE = 0
C     IF(ABS(NQPAI+NQN).EQ.0) THEN
C     XMULTP = ( FMUL0(S)-0.5 )*3./2. + 0.5
C     RETURN
C     END IF
C
C ++ TOTAL CHARGE = 1
C     IF(ABS(NQPAI+NQN).EQ.1) THEN
      XMUL0  = ( FMUL0(S)-0.5 )*3./2. + 0.5
      XMUL2  = ( FMUL2(S)-0.5 )*3./2. + 0.5
      XMULTP = ( XMUL0+XMUL2 ) / 2.
      RETURN
C     END IF
C
C ++ TOTAL CHARGE = 2
C     IF(ABS(NQPAI+NQN).EQ.2) THEN
C     XMULTP = ( FMUL2(S)-0.5 )*3./2. + 0.5
C     RETURN
C     END IF
C
C
C     RETURN
      END
