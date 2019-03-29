***********************************************************************
*     -------------------
      SUBROUTINE TAINIPAR
*     -------------------
*
*     ( purpose )
*       set initial parameter for tauola
*
*     ( input )
*       none
*
*     ( output )
*       COMMON  /DECPAR/,/QEDPRM/,/TAURAD/,/INOUT/,/IDPART/
*
*     ( creation date and author )
*       1995.10.24 ; First version by K.Kaneyuki
*
*     ( comment )
*
***********************************************************************
#include "taucom.h"
      REAL*8 PI8

C
C     COMMON /DECPAR/
C
      GFERMI = 1.16637E-5
      GV     = 1.0
      GA     = -1.0
      CCABIB = 0.975
      SCABIB = SQRT(1.-CCABIB**2)
      PI     = 4.0*ATAN(1.)
      GAMEL  = GFERMI**2*AMTAU**5/(192.0*PI**3)

C
C     COMMON /QEDPRM/
C
      PI8    = 4.D0*DATAN(1.D0)
      ALFINV = 137.03604D0
      ALFPI  = 1D0/(ALFINV*PI8)
      XK0    = 0.1D0

C
C     COMMON /TAURAD/
C
      XK0DEC = 1.0D-3
      ITDKRC = 1

C
C     COMMON /INOUT/
C
      INUT = 5
      IOUT = 6
C
C     COMMON /IDFC/
C
      IDFF = -15

C
C     COMMON /IDPART/
C
C      IA1 = 20213

      RETURN
      END

