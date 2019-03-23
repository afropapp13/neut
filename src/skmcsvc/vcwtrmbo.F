**********************************************************************
*     ----------------------------------
      SUBROUTINE VCWTRMBO(N,ET,XM,WTMAX)
*     ----------------------------------
*
*     ( purpose )
*       CALCULATE MAX WEIGHT FOR RAMBO
*
*     ( input )
*        N      : NUMBER OR PARTICLES (>1, IN THIS VERSION <101)
*        ET     : TOTAL CENTORE-OF-MASS ENERGY
*        XM     : PARTICLE MASSES ( DIM=100 )
* 
*     ( output )
*        WTMAX  : MAX WEIGHT OF THE EVENT
*
*     ( creation date and author )
*       1989. 7.25 ; First version by K.Kaneyuki
*       1992. 8. 6 ; mod for S.K.
*                    WTRMBO -> VCWRMB      
*       1995. 2. 7 ; change names
*                    VCWRMB -> VCWTRMBO
*
**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XM(N),Z(100)
      DATA IBEGIN/0/
C     
C --  INITIALIZATION STEP: FACTORIALS FOR THE PHASE SPACE WEIGHT
C     
      IF (IBEGIN.EQ.0) THEN
         IBEGIN=1
         TWOPI=8.*DATAN(1.D0)
         PO2LOG=DLOG(TWOPI/4.)
         Z(2)=PO2LOG
         DO 101 K=3,100
            Z(K)=Z(K-1)+PO2LOG-2.*DLOG(DFLOAT(K-2))
  101    CONTINUE
         DO 102 K=3,100
            Z(K)=(Z(K)-DLOG(DFLOAT(K-1)))
  102    CONTINUE
      ENDIF
C     
C --  CHECK WHETHER TOTAL ENERGY IS SUFFICIENT; COUNT NONZERO MASSES
C     
      XMT=0.
      NM=0
      DO 105 I=1,N
         IF (XM(I).NE.0.D0) NM=NM+1
         XMT=XMT+DABS(XM(I))
  105 CONTINUE
C     
C --  RETURN FOR WEIGHTED MASSLESS MOMENTA
C     
      WTMAX=PO2LOG
      IF (N.NE.2) WTMAX=(2.*N-4)*DLOG(ET)+Z(N)
      IF (NM.EQ.0) THEN
         WTMAX=DEXP(WTMAX)
         RETURN
      ENDIF
C     
C --  RETURN FOR WEIGHTED MASSIVE MOMENTA
C     
      XMAX=DSQRT(1.-(XMT/ET)**2)
      WTM=(2.*N-3.)*DLOG(XMAX)
      WTMAX=WTMAX+WTM
      WTMAX=DEXP(WTMAX)

      RETURN
      END
