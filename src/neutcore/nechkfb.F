**********************************************************************
*     ------------------------------------
      SUBROUTINE NECHKFB(W,MULP,NFPI,ICHK)
*     ------------------------------------
*     
*     ( purpose )
*       CHECK THE FORWARD AND BACKWARD PION MULTIPLICITY
*
*     ( input ) 
*       W    ; INVAVIANT MASS OF HAD. SYSTEM
*       MULP ; TOTAL PION MULTI PLICITY
*       NFPI ; FORWARD PION MULTIPLICITY
*
*     ( output ) 
*        OUTPUT  ICHK ; =0 NO (THIS F. B. DIST.IS NOT GOOD -> TRY AGAIN)
*                     ; =1 OK
*
*     ( creation date and author )
*       1985.11.09 ; First version by T.KAJITA
*       1995.03.03 ; K.Kaneyuki for S.K.
*                    CHKFB -> NECHKFB
*
**********************************************************************
C#include "nemacr.h"
C      COMMON /NVECT/ IMOD,LVECT,IP(20),PIN(3,20)
      DIMENSION PRODST(10)
      DIMENSION EQUDST(10),WEIDST(10)
C
C     EQUDST -- PROB. OF TAKING N FORWARD MULTIP. OUT OF MULP
C     ASSUMING THE EQUAL F. AND B. DIST.
C     WEIDST -- PROB. OF TAKING N FORWARD MULTIP. OUT OF MULP
C     ASSUMING THE F. AND B. AN-ASYMMETRY
C     
      ICHK=0
      CALL VZERO(EQUDST,10)
      CALL VZERO(WEIDST,10)
      CALL VZERO(PRODST,10)
      IF(W.GT.1.4) GO TO 1
      ICHK=1
      RETURN
C     
    1 DO 10 I=0,MULP
         XN=FLOAT(MULP+1)
         XR=FLOAT(I+1)
         XNR=FLOAT(MULP-I+1)
         EQUDST(I+1)=GAMMA(XN)/GAMMA(XR)/GAMMA(XNR)*0.5**MULP
   10 CONTINUE
C     
C***  IF(IP(1).GT.0 ) THEN
C     + NEUTRINO  ( REF. D.ZIEMINSKA ET AL., P.R.D27 (1983) 47 )
C     HERE THE ASSUMPTION IS NEU+P AND NEU+NTHE SAME RESULT, SEE REF.
C***  AVPIF=0.975*ALOG10(W**2)-0.35     +0.35
C***  AVPIB=0.50*ALOG10(W**2)-0.10     +0.35
C     THIS 0.35 HAS NO STRONG MEANING
C     
C***  ELSE
C     + ANTI NEUTRINO ( REF. S.BARLAG ET AL., Z.PHYS.C11 (1982) 283 )
C     HERE THE ASSUMPTION IS NEU+P AND NEU+NTHE SAME RESULT, SEE REF.
      AVPIF=0.95*ALOG10(W**2)+0.35
      AVPIB=0.225*ALOG10(W**2)+0.50
C     
C***  END IF
C     
      IF(AVPIF.LE.0.1 .OR. AVPIB.LE.0.1) PRF=0.5
      PRF = AVPIF / (AVPIF+AVPIB)
      IF(PRF.LT.0.5) PRF=0.5
C     
      DO 20 I=0,MULP
         XN=FLOAT(MULP+1)
         XR=FLOAT(I+1)
         XNR=FLOAT(MULP-I+1)
         WEIDST(I+1)=GAMMA(XN)/GAMMA(XR)/GAMMA(XNR)
     &        *PRF**I * (1.-PRF)**(MULP-I)
   20 CONTINUE
C     
      PROMAX=0.
      DO 30 I=1,MULP+1
         PRODST(I)= WEIDST(I)/EQUDST(I)
         IF(PRODST(I).GT.PROMAX)PROMAX=PRODST(I)
   30 CONTINUE
C     
C     
      PRO=PRODST(NFPI+1)/PROMAX
      IF(PRO.GT.RLU(DUM)) ICHK=1
C     
      RETURN
      END
