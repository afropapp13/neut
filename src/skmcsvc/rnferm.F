***********************************************************************
*     ------------------------------
      SUBROUTINE RNFERM(PFERM,IWAVE)
*     ------------------------------
*     
*     ( purpose )
*       Fermi motion in O16
*
*     ( input ) 
*       none
*
*     ( output ) 
*       PFERM : Fermi momentum ( MeV/c )
*       IWAVE : S-WAVE - 0
*               P-WAVE - 1
*
*     ( creation date and author )
*       1983.12.12 ; First version by M.Nakahata
*                  ; Mod by T.Kajita ( add IWAVE )
*       1992. 9.18 ; Mod by K.Kaneyuki for S.K.
*       2001. 2. 9 ; RNDM -> RLU by K.K.
*
*     ( comment )
*
***********************************************************************

#include "necard.h"
      
      INTEGER IDUMMY

      PFERM=0

      if (NEFRMFLG.EQ.1) RETURN

C      IF (RNDM(DUMMY).LT.0.25) THEN
      IF (RLU(IDUMMY).LT.0.25) THEN
         CALL RNSWAV(PFERM)
         IWAVE=0
      ELSE
         CALL RNPWAV(PFERM)
         IWAVE=1
      ENDIF

      RETURN
      END
