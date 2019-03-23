************************************************************************
*     ------------------------------------------------
      SUBROUTINE NEGAMNUC(MODEGM,NGAMMA,EGAMMA,DIRGAM,*)
*     ------------------------------------------------
C
*     (Purpose)
C     
C       DETERMINE NUMBER OF GAMMAS(NGAMMA), ENERGY OF GAMMAS(EGAMMA(10))
C       AND DIRERCTION OF GAMMAS (DIRGAM(3,10)) FROM NUCLEAR EXITATION.
C
*     (Input)
C       INPUT ; MODEGM -- MODE OF NUCLEAR EXITATION GAMMA
C               ( SEE NEMODGAM.F)
C     
*     (Output)
*       NGAMMA         ; NUMBER of GAMMAS
*       EGAMMA(MAXVC)  ; ENERGY of each GAMMA
*       DIRGAM(3,MAXVC); DIRECTION of each GAMMA
C
*     (Creation Date and Author)
C       08-JAN-90    T.KAJITA
C       23-JUN-96    Y.Hayato (S.K. version)
C
C     
#include "vcvrtx.h"
#include "vcwork.h"

      DIMENSION EGAMMA(MAXVC), DIRGAM(3,MAXVC), LIST(MAXVC), DIR3(3)
      REAL*8  ETEMP(MAXVC)
C     
      IF (MODEGM.LT.0) RETURN 1
      NGAMMA=0
      IF (MODEGM.EQ.0) RETURN
C     
      IF(MODEGM.EQ.1) THEN
         NGAMMA = INT (2.*RLU(DUM)) +1
         ETOT=5.+6.*RLU(DUM)
         IF(NGAMMA.EQ.1) THEN
            EGAMMA(1)=ETOT
         ELSE
            EGAMMA(1)=ETOT*RLU(DUM)
            EGAMMA(2)=ETOT - EGAMMA(1)
         END IF
         GO TO 1000
      END IF
C     
      ETOT=5.+15.*RLU(DUM)
      IF(ETOT.LT.10.) THEN
         NGAMMA = INT (2.*RLU(DUM)) +1
      ELSE
 46      NGAMMA = INT (3.*RLU(DUM)) +2
         IF(ETOT.GT.16. .AND. NGAMMA.EQ.1) GO TO 46
      END IF
      IF(NGAMMA.EQ.1) THEN
         EGAMMA(1)=ETOT
      ELSE
         DO 30 I=1,(NGAMMA-1)
            ETEMP(I)=ETOT*RLU(DUM)
            LIST(I)=I
 30      CONTINUE
         ETEMP(NGAMMA)=ETOT
         LIST(NGAMMA)=NGAMMA
         CALL SORTDL(ETEMP,LIST,1,NGAMMA,1)
         DO 31 I=1,NGAMMA
            IF(I.EQ.1) EGAMMA(I)=ETEMP(LIST(I))
            IF(I.NE.1) EGAMMA(I)=ETEMP(LIST(I))-ETEMP(LIST(I-1))
 31      CONTINUE
      END IF
C     
C     
 1000 CONTINUE
      DO 1001 I=1,NGAMMA
         CALL RNDIR(DIR3)
         DIRGAM(1,I)=DIR3(1)
         DIRGAM(2,I)=DIR3(2)
         DIRGAM(3,I)=DIR3(3)
 1001 CONTINUE
C     
C     
      RETURN
      END
