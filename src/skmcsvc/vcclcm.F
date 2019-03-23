C     
      SUBROUTINE VCCLCM
C     
C     ( purpose )
C       Clear common VCWORK,VCVRTX
C
C     ( input ) 
C       none
C
C     ( output ) 
C       COMMON /VCWORK/,/VCVRTX/
C
C     ( creation date and author )
C       1992.6.5 ; First version by K.Kaneyuki
C
C     ( comment )
C      
C
#include "vcwork.h"
#include "vcvrtx.h"

      NVC=0
      CALL VZERO(POSVC ,3)
      CALL VZERO(IPVC  ,MAXVC)
      CALL VZERO(AMASVC,MAXVC)
      CALL VZERO(PVC   ,3*MAXVC)
      CALL VZERO(IORGVC,MAXVC)
      CALL VZERO(IFLGVC,MAXVC)
      CALL VZERO(ICRNVC,MAXVC)
      CALL VZERO(TIMVC ,MAXVC)
      CALL VZERO(POSIVC,3*MAXVC)
      CALL VZERO(IVTIVC,MAXVC)
      CALL VZERO(POSFVC,3*MAXVC)
      CALL VZERO(IVTFVC,MAXVC)

      NVTXVC=0
      CALL VZERO(PVTXVC,3*MAXVC)
      CALL VZERO(IFLVVC,MAXVC)
      CALL VZERO(IPARVC,MAXVC)
      CALL VZERO(TIMVVC,MAXVC)

      RETURN
      END
