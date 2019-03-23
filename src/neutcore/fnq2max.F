************************************************************************
*     ------------------------
      FUNCTION FNQ2MAX(E,IPAR)
*     ------------------------
*
*     (Purpose)
*       Calc. maximum value of Q**2
*
*     (Input)
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
*       FNQ2MAX : maximum value of Q**2
*
*     (Creation Date and Author)
*       1983.11.20 ; M.Nakahata
*       1987.05.21 ; N.Sato
*       1995.01.30 ; K. KANEYUKI  FOR S.K.
*                    Q2MAX -> FNQ2MAX
*     
************************************************************************
C
C    PM    ; PROTON MASS
C    EM    ; ELECTRON MASS
C    UM    ; MUON MASS
C    TM    ; TAU MASS
C     
C      COMMON /TEST/ENSTOP
      DATA PM/0.938/
      DATA EM/0.511E-3/
      DATA UM/0.1057/
C      DATA TM/1.7771/
CCCC  DATA TM/1.7842/     MODIFICATION  K.Kaneyuki 2-FEB-95
      DATA TM/1.7842/
C     
      IF(IABS(IPAR).EQ.12)AM=EM
      IF(IABS(IPAR).EQ.14)AM=UM
      IF(IABS(IPAR).EQ.16)AM=TM
C     
      A=1.+2.*E/PM
      B=2.*E**2-(1.+E/PM)*AM**2
      C=AM**4
      
C**********************************
C     TEST
C**********************************
      TESTV=B**2-A*C
      IF(TESTV.LT.0.) THEN
         PRINT *,' TESTV< 0 IN Q2MAX ********',ENSTOP
         FNQ2MAX=0.
         RETURN
      END IF
C**********************************
      FNQ2MAX=(B+SQRT(B**2-A*C))/A
      RETURN
      END
