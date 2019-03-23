************************************************************************
*     ---------------------------------
      SUBROUTINE NEELSVCT(IPAR,E,DNEUT)
*     ---------------------------------
*
*     (Purpose)
*       Vector generation for elastic event
*
*     (Input)
*       E      : NEUTRINO ENERGY ( GEV )
*       IPAR   : NEUTRINO TYPE
*                  12 : NEU-E
*                 -12 : NEU-E-BAR            
*                  14 : NEU-MU
*                 -14 : NEU-MU-BAR            
*                  16 : NEU-TAU
*                 -16 : NEU-TAU-BAR
*       DNEUT(3) : DIRECTION OF NEUTRINO
*     
*     (Output)
*       COMMON /NEWORK/
*
*     (Creation Date and Author)
*       19??.??.?? ; M.NAKAHATA
*       1987.06.23 ; N.SATO     INCLUDE TAU
*       1995.01.31 ; K. KANEYUKI  FOR S.K.
*                    ELSVCT -> NEELSVCT
*     
************************************************************************
C
#include "nework.h"
      DIMENSION DNEUT(3),DLEP(3)
      INTEGER IPMODE(4,6)/ 12, 2112, 11, 2212, 14, 2112, 13, 2212,
     &                     16, 2112, 15, 2212,-12, 2212,-11, 2112,
     &                    -14, 2212,-13, 2112,-16, 2212,-15, 2112/
      INTEGER IORG(4)/ 0, 0, 1, 2/
      INTEGER IFLG(4)/-1,-1, 0, 0/
      INTEGER ICRN(4)/ 0, 0, 1, 1/
      DATA PMASS /0.938/
C     
      IF(IPAR.EQ. 12)II=1
      IF(IPAR.EQ. 14)II=2
      IF(IPAR.EQ. 16)II=3
      IF(IPAR.EQ.-12)II=4
      IF(IPAR.EQ.-14)II=5
      IF(IPAR.EQ.-16)II=6
C     
C     -- SET NEUTRINO DIRECTION
C     
C      CALL NERNNEUT(E,DNEUT,IPAR)
C     
C     -- SET Q**2
C     
   10 Q2=RNELQ2(IPAR,E)
C     
C     -- EXECUTE OTHER VECTORS
C     
      EL=E-Q2/2./PMASS

      IF(IABS(IPAR).EQ.12) IPLEPT=11
      IF(IABS(IPAR).EQ.14) IPLEPT=13
      IF(IABS(IPAR).EQ.16) IPLEPT=15
      CALL MCMASS(IPLEPT,AMLEP)
      AMLEP=AMLEP*1.0E-3
      PL=SQRT(EL**2-AMLEP**2)
      COST=(2.*E*EL-AMLEP**2-Q2)/2./E/PL
      IF(ABS(COST).GT.1.) THEN
         WRITE(6,600)COST
 600     FORMAT(' *** WARNING ***  IN NEELSVCT  COST=',G15.7)
         write(6,*) IPAR,COST,PL,AMLEP,E,EL,Q2,PL
         GO TO 10
      ENDIF   

      CALL RNROT(COST,DNEUT,DLEP)
C     
C     -- NOW STORE VECTOR
C     
      IF( IPAR.GT.0 ) MODENE= 1
      IF( IPAR.LT.0 ) MODENE=-1
      NUMNE = 4
      DO 20 I=1,4
         IPNE(I)=IPMODE(I,II)
         IORGNE(I)=IORG(I)
         IFLGNE(I)=IFLG(I)
         ICRNNE(I)=ICRN(I)
         DO 30 J=1,3
            IF (I.EQ.1) PNE(J,I)=E*DNEUT(J)
            IF (I.EQ.2) PNE(J,I)=0.0
            IF (I.EQ.3) PNE(J,I)=PL*DLEP(J)
            IF (I.EQ.4) PNE(J,I)=PNE(J,1)+PNE(J,2)-PNE(J,3)
 30      CONTINUE
 20   CONTINUE

      RETURN
      END
