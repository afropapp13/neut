************************************************************************
*     --------------------------------------
      SUBROUTINE NENCELSVNT(IPAR,E,MODE,DNEUT)
*     --------------------------------------
*
*     (Purpose)
*       VECTOR GENERATION FOR N.C. ELASTIC EVENT ON FREE NUCLEONS
*
*     (Input)
*       IPAR   : NEUTRINO TYPE
*                  12 : NEU-E
*                 -12 : NEU-E-BAR            
*                  14 : NEU-MU
*                 -14 : NEU-MU-BAR            
*                  16 : NEU-TAU
*                 -16 : NEU-TAU-BAR
*       E      : NEUTRINO ENERGY ( GEV )
*       DNEUT(3) : DIRECTION OF NEUTRINO
*     
*     (Output)
*       COMMON /NEWORK/
*
*     (Creation Date and Author)
*       2011.03.15 ; Y. Hayato ( NC vector gen. to use correct formula )
*     
************************************************************************
      IMPLICIT NONE

#include "nework.h"
#include "neutmodel.h"

      INTEGER*4 IPAR,MODE
      REAL E,DNEUT(3)      
      
      INTEGER IORG(4)/ 0, 0, 1, 2/
      INTEGER IFLG(4)/-1,-1, 2, 0/
      INTEGER ICRN(4)/ 0, 0, 1, 1/

      INTEGER*4 ITARG,I,J

      REAL*4    PMASS
      REAL*4    Q2,EL,PL,COST,DLEP(3)

      REAL*4    RNNCELQ2
      EXTERNAL  RNNCELQ2


      if (  mod(MDLQE,100)/10.eq.0 ) then
         write(*,*) 'Invalid QE model (MDLQE) selected.',MDLQE
         STOP
      endif

      NUMNE = 4
      IPNE(1)=IPAR
      IPNE(3)=IPAR
      IF (ABS(MODE).EQ.51) THEN
         IPNE(2)=2212
         IPNE(4)=2212
      ELSE IF (ABS(MODE).EQ.52) THEN
         IPNE(2)=2112
         IPNE(4)=2112
      else
         write(*,*) 'Invalid mode for NENCELSVNT',MODE
         stop
      ENDIF
      
      DO 10 I=1,4
         IORGNE(I)=IORG(I)
         IFLGNE(I)=IFLG(I)
         ICRNNE(I)=ICRN(I)
 10   CONTINUE
      
      MODENE = MODE

C---- kinematics part ------------

      ITARG = IPNE(2)
      CALL MCMASSGV(ITARG,PMASS)

   20 Q2=RNNCELQ2(IPAR,E,ITARG)
      
      EL=E-Q2/2./PMASS

      PL=EL
      COST=(2.*E*EL-Q2)/2./E/PL
      IF(ABS(COST).GT.1.) THEN
         WRITE(6,600)COST
 600     FORMAT(' *** WARNING ***  IN NENCELSVNT COST=',G15.7)
         write(6,*) IPAR,COST,PL,E,EL,Q2,PL
         GO TO 20
      ENDIF   

      CALL RNROT(COST,DNEUT,DLEP)
C     
C     -- NOW STORE VECTOR
C     
      DO 40 I=1,4
         DO 30 J=1,3
            IF (I.EQ.1) PNE(J,I)=E*DNEUT(J)
            IF (I.EQ.2) PNE(J,I)=0.0
            IF (I.EQ.3) PNE(J,I)=PL*DLEP(J)
            IF (I.EQ.4) PNE(J,I)=PNE(J,1)+PNE(J,2)-PNE(J,3)
 30      CONTINUE
 40   CONTINUE

      RETURN
      END
