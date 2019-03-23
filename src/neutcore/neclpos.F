**********************************************************************
*     -------------------
      SUBROUTINE NECLPOS
*     -------------------
*     
*     ( purpose )
*       Clear common POSINNUC
*
*     ( input ) 
*       none
*
*     ( output ) 
*       COMMON /POSINNUC/
*
*     ( creation date and author )
*       2006.6.16 ; Y.Hayato
*
*     ( comment )
*      
**********************************************************************

#include "vcwork.h"
#include "posinnuc.h"

      integer*4 I,J

      IBOUND = -1
      do 100 I=1,MAXVC
         do 200 J=1,3
            POSNUC(J,I) = -9999.
 200     continue
 100  continue

      RETURN
      END
