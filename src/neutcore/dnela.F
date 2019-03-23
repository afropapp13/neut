************************************************************************
*     --------------------------
      REAL*8 FUNCTION DNELA(Q2)
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
*       1983.11.20 ; M.Nakahata
*       1987.05.21 ; N.Sato 
*       1995.01.30 ; K. KANEYUKI  FOR S.K.
*                    FUNELA -> DNELA      
*
************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL DNELSQ2
      COMMON /NEINT3/ E,IPAR
      DNELA=DNELSQ2(E,IPAR,Q2)

      RETURN
      END
