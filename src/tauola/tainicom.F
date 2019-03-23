***********************************************************************
*     -------------------
      SUBROUTINE TAINICOM
*     -------------------
*     
*     ( purpose )
*       set initial parameter for tauola
*
*     ( input ) 
*       none
*
*     ( output ) 
*       COMMON  /DECPAR/ /PARMAS/ /TAUBRA/ /TAUNPI/ /JAKI/ /IDFC/
*               /INOUT/ /IDPART/ /TAURAD/ /QEDPRM/ /TAUKLE/
*
*     ( creation date and author )
*       1995.10.24 ; First version by K.Kaneyuki
*
*     ( comment )
*
***********************************************************************
#include "taucom.h"
C
C     READ CARD
C
      CALL TACARD
C      
C     COMMON /PARMAS/
C
      CALL TAINIMAS
C      
C     COMMON /DECPAR/,/QEDPAR/
C
      CALL TAINIPAR
C      
C     COMMON  /TAUBRA/,/DECPAR/,/TAUKLE/
C
      CALL TAINIDCY

      RETURN
      END
      
