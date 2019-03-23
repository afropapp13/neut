**********************************************************************
*     ------------------------
      SUBROUTINE NECPNEWK(POS)
*     ------------------------
*     
*     ( purpose )
*       COPY NEWORK TO VCWORK, VCVRTX
*
*     ( input )
*       POS(3)  : INTERACTION POSITION      
*       COMMON /NEWORK/
*
*     ( output ) 
*       COMMON /VCWORK/
*       COMMON /VCVRTX/
*
*     ( creation date and author )
*       1995.03.15 ; First version by K.Kaneyuki
*       2006.02.23 ; Add implicit none
*
**********************************************************************
      IMPLICIT NONE

#include "vcwork.h"
#include "vcvrtx.h"
#include "nework.h"
      REAL POS(3) 

      integer*4 I,J

C
C --  COPY NEWORK TO VCWORK
C     
      NVC=NUMNE
      CALL UCOPY(POS,POSVC,3)
      
      DO 20 I=1,NUMNE
         IPVC(I)=IPNE(I)
         CALL MCMASS(IPNE(I),AMASVC(I))
         DO 30 J=1,3
            PVC(J,I)=PNE(J,I)*1.0E+3
 30      CONTINUE
         IORGVC(I)=IORGNE(I)
         IFLGVC(I)=IFLGNE(I)
         ICRNVC(I)=ICRNNE(I)
         TIMVC(I) =0.0
         CALL UCOPY(POS,POSIVC(1,I),3)
         IVTIVC(I)=1
         IF (IFLGNE(I).EQ.-1 .OR. IFLGNE(I).EQ.1 .OR.
     &       IFLGNE(I).EQ.3  .OR. IFLGNE(I).EQ.4 .OR.
     &       IFLGNE(I).EQ.5 ) THEN
            CALL UCOPY(POS,POSFVC(1,I),3)
            IVTFVC(I)=1
         ENDIF
 20   CONTINUE
C
C --  COPY NEWORK TO VCVRTX
C     
      NVTXVC=1
      CALL UCOPY(POS,PVTXVC(1,1),3)
      IFLVVC(1)=1
      IPARVC(1)=0
      TIMVVC(1)=0.0
      
      RETURN
      END
