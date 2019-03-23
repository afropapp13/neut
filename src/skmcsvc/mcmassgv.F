***********************************************************************
*     -----------------------------
      SUBROUTINE MCMASSGV(IP,XMASS)
*     -----------------------------
*     
*     ( purpose )
*       set mass in GeV
*
*     ( input ) 
*       IP : particle code
*
*     ( output ) 
*       XMASS : mass ( GeV/c**2 )
*
*     ( creation date and author )
*       1995.10.24 ; First version by K.Kaneyuki
*
*     ( comment )
*
***********************************************************************
      IMPLICIT NONE

      INTEGER*4 IP
      REAL*4 XMASS

      CALL MCMASS(IP,XMASS)
      XMASS = XMASS * 1.0E-3
      RETURN
      END
      
