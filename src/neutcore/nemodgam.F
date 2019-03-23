************************************************************************
*     ------------------------------------------------
      FUNCTION NEMODGAM(DUMMY)
*     ------------------------------------------------
C
*     (Purpose)
C       RETURN MODE OF NUCLEAR EXCITATION GAMMA 
C
C     (Input)
C       None
C      
C     (Output)
C       MODE OF NUCLEAR EXITETION GAMMA
C             NEMODGAM= -1    ERROR
C                     =  0    NO GAMMA (INTERACTION ON FREE PROTON)
C                     =  1    GAMMAS FROM 1 NUCLEOUS
C                     =  2    GAMMAS FROM 1 OR 2 NUCLEOUS
C
*     (Creation Date and Author)
C       08-JAN-90    T.KAJITA
C       23-JUN-96    Y.Hayato (S.K. version)
C
C
C
C      COMMON /GTVECT/ POSV(3),IPV(50),PINV(3,50),PABSV(50)
C      COMMON /NVECT/ IMOD,LVECT,IPORI(20),PINORI(3,20)
C
#include "nework.h"
#include "vcwork.h"
#include "vcvrtx.h"

      NEMODGAM= -1
C      CALL GETVCT
C      CALL GTNEUT(INPAR,&100)
C
C++ INTERACTION ON FREE PROTON
      NEMODGAM=0
      IF(sqrt(PVC(1,2)**2+PVC(2,2)**2+PVC(3,2)**2)
     $   .LT.0.0001) RETURN
C
C++ ELASTIC INTERACTION ON BOUND NUCLEON
      NEMODGAM=1
      IF(IABS(MODENE).EQ. 1) RETURN
      IF(IABS(MODENE).EQ.51) RETURN
      IF(IABS(MODENE).EQ.52) RETURN
C
C++ CHECK IF PIONS ARE ABSORPED WITHIN THE NUCLEUS
C--  ( Add delta absorption case / 091005 Y.Hayato )
      DO 10 I=1,MAXVC
         IF ((IFLGVC(I).eq.3).and.
     $        ((abs(IPVC(I)).eq.211).or.(IPVC(I).eq.111).or.
     $         (IPVC(I).eq.2224).or.(IPVC(I).eq.2214).or.
     $         (IPVC(I).eq.2114).or.(IPVC(I).eq.1114))) THEN
            GOTO 15
         ENDIF
C     10  IF(IPV(I).EQ.2000) GO TO 15
 10   CONTINUE
      RETURN
C
C++ ABSORPTION OF PIONS
  15  NEMODGAM=2
C
 100  RETURN
      END
