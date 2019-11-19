************************************************************************
*     --------
*     RESOLVECARDNAME.F
*     --------
*
*     (Purpose)
*       FUNCTION for getting the right cardname
*
*     (Creation Date and Author)
*       2019.10.23 ; L. Pickering
************************************************************************
      SUBROUTINE RESOLVECARDNAME

      implicit none

#include "cardname.h"

!     If an external function has already set the card name, then use that
      IF (ISSET.EQV..FALSE.) THEN
!       Otherwise, check the gfortran runtime CLI args
        IF (IARGC() .GE. 1) THEN
          CALL GETARG(1, FCARD)
        ELSE
!         If all else fails, use a default
          FCARD = 'neut.card'
        END IF
        ISSET = .TRUE.
      END IF
      WRITE(*,*) "Reading ", FCARD(1:LEN_TRIM(FCARD)), " as card"
      END
