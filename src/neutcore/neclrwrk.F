**********************************************************************
*     -------------------
      SUBROUTINE NECLRWRK
*     -------------------
*     
*     ( purpose )
*       Clear common NEWORK
*
*     ( input ) 
*       none
*
*     ( output ) 
*       COMMON /NEWORK/
*
*     ( creation date and author )
*       1995.2.26 ; First version by K.Kaneyuki
*
*     ( comment )
*      
**********************************************************************

#include "nework.h"

      MODENE = 0
      NUMNE = 0
      CALL VZERO(IPNE  ,MAXNE)
      CALL VZERO(PNE   ,3*MAXNE)
      CALL VZERO(IORGNE,MAXNE)
      CALL VZERO(IFLGNE,MAXNE)
      CALL VZERO(ICRNNE,MAXNE)

      RETURN
      END
