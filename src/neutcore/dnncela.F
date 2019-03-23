************************************************************************
*     --------------------------
      REAL*8 FUNCTION DNNCELA(Q2)
*     --------------------------
*
*     (Purpose)
*       elastic without pauli principle 
*
*     (Input)
*       Q2     : Q**2
*       E      : NEUTRINO ENERGY ( GEV )      
*       IPAR   : NEUTRINO TYPE
*                  12 : NEU-E
*                 -12 : NEU-E-BAR            
*                  14 : NEU-MU
*                 -14 : NEU-MU-BAR            
*                  16 : NEU-TAU
*                 -16 : NEU-TAU-BAR
*     
*     (Output)
*       DNELA  : Q**2 DISTRIBUTION OF ELASTIC SCATTERING
*
*     (Creation Date and Author)
*       2011.03.12 ; Y.Hayato
*                    FNNCELA ( modify FNELA )
*
************************************************************************
      IMPLICIT NONE

      REAL*8   Q2

      REAL*8   DNNCELSQ2
      EXTERNAL DNNCELSQ2

      REAL*8   E
      INTEGER  IPAR,ITARG
      COMMON /NEINT3/ E,IPAR,ITARG
      DNNCELA=DNNCELSQ2(E,IPAR,ITARG,Q2)

      RETURN
      END
