***********************************************************************
*     -------------------
      SUBROUTINE TAINIMAS
*     -------------------
*
*     ( purpose )
*       set mass for tauola
*
*     ( input )
*       none
*
*     ( output )
*       COMMON  /PARMAS/
*
*     ( creation date and author )
*       1995.10.24 ; First version by K.Kaneyuki
*                  ; bug-fix by F.Kato
*
*     ( comment )
*
***********************************************************************
#include "taucom.h"

C
C     COMMON /PARMAS/
C

C     LEPTON

      CALL MCMASSGV(15,AMTAU)
      AMNUTA = 0.010
      CALL MCMASSGV(13,AMMU)
      AMNUMU = 0.0
      CALL MCMASSGV(11,AMEL)
      AMNUE  = 0.0

C     PI

      CALL MCMASSGV(111,AMPIZ)
      CALL MCMASSGV(211,AMPI)

C     RHO, A1

      CALL MCMASSGV(213,AMRO)
      GAMRO  = 0.145
      AMA1 = 1.251
      GAMA1 = 0.599

C     K, K*

      CALL MCMASSGV(321,AMK)
      CALL MCMASSGV(311,AMKZ)
      CALL MCMASSGV(323,AMKST)
      GAMKST = 0.0513

      RETURN
      END
