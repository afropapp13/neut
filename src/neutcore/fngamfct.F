**********************************************************************
*     -----------------------
      FUNCTION FNGAMFCT(BETA)
*     -----------------------
*     
*     ( purpose )
*       GET GAMMA FACTOR FROM BETA (V/C)
*
*     ( input ) 
*       BETA    : V/C
*
*     ( output ) 
*       FNGAMFCT : GAMMA FACTORR
*
*     ( creation date and author )
*       1985.10.19 ; First version by T.KAJITA
*       1995.03.03 ; K.Kaneyuki for S.K.
*                    GAMFCT -> FNGAMFCT
*
**********************************************************************

      IF (ABS(BETA).LT.1.) THEN
         FNGAMFCT = 1. / SQRT(1.-BETA**2)
      ELSE
         FNGAMFCT=10000.
      ENDIF

      RETURN
      END
