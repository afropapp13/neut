***********************************************************************
*     ---------------------------------------------------
      FUNCTION RSSGCRS(E,IPAR,ITYPE)
*     ---------------------------------------------------
*
*     ( purpose )
*       return cross-section of single-pion production
*
*     ( input )
*       E           : Energy of incident neutrino
*       IPAR        : Kind of incoming neutrino
*       ITYPE       : interaction type
*
*     ( output )
*       cross-section for each mode
*
*     ( creation date and author )
*       1993.Dec. ; for Super-Kamioka by Y.Hayato
*                   (Based on Rein and Sehgal's code)
*       1995.Oct. ; tau-neutrino(temporary)
*       2007.08.23 ; G.Mitsuka add delta->gamma decay
*
*     ( comment )
*       Reference:Rein and Sehgal,Ann. of Phys. 133(1981),79-153
*                 Rein           ,Z.Phys.C 35(1987),43-64
*                 Feynman et al. ,Phys.Rev.D3(1971),2706-
*
*        ----------------
*        IPAR     : Particle CODE
*                   Currentry it doesn't have no meaning
*             = 12: nu-e
*             =-12: nu-e-bar
*             = 14: nu-my
*             =-14: nu-mu-bar
*             = 16: nu-tau
*             =-16: nu-tau-bar
*        ITYPE    : REACTION TYPE
*             *** CHARGED CURRENT ***
*             =1  : NEU,N    --> LEPTON-,P,gamma
*             =2  : NEUBAR,P --> LEPTON+,N,gamma
*             *** NEUTRAL CURRENT ***
*             =11 : NEU,N --> NEU,N,gamma
*             =12 : NEU,P --> NEU,P,gamma
*             =13 : NEUBAR,N --> NEUBAR,N,gamma
*             =14 : NEUBAR,P --> NEUBAR,P,gamma
*
*        E        : NEUTRINO ENERGY(GeV)
***********************************************************************

      IMPLICIT NONE

      REAL RSSGCRS
      REAL E
      INTEGER IPAR,ITYPE
      

      
      INTEGER ICONV(18)
      
      DATA ICONV/ 1,11,99,99,99,99,99,99,99,99,
     $            3, 2,12,13,99,99,99,99/

      REAL SIGMA,SIGMA3

      INTEGER IMODE
      
      IF ((ITYPE.GT.14).or.(ITYPE.LT.1)) GO TO 999

      IMODE = ICONV(ITYPE)

      IF (IMODE.EQ.99) GO TO 999

C  CHECK KIND OF PARTICLE (ADD July 13/Y.H)
      IF ((IPAR.gt.0).and.(IMODE.gt.10)) GO TO 9999
      IF ((IPAR.lt.0).and.(IMODE.lt.10)) GO TO 9999

C-MODEFY TO SUPPORT e/mu/tau differences
C      IF (ABS(IPAR).ne.16) THEN
C         CALL RSTCRS(IMODE,E,SIGMA,SIGMA3)

C  FOR ELECTRON NEUTRINO
      IF (ABS(IPAR).eq.12) THEN
         CALL RSTCRSG(IMODE,E,SIGMA,SIGMA3)
         GOTO 500
      ENDIF
C  FOR MU NEUTRINO
      IF (ABS(IPAR).eq.14) THEN
         CALL RSTCRSGM(IMODE,E,SIGMA,SIGMA3)
         GOTO 500
      ENDIF
C  FOR TAU NEUTRINO
      IF (ABS(IPAR).eq.16) THEN
         CALL RSTCRSGT(IMODE,E,SIGMA,SIGMA3)
         GOTO 500
      ENDIF

 500  RSSGCRS=SIGMA

      RETURN

 999  WRITE(*,*)
     $ 'RSSGCRS:*** Invalid ITYPE(Interaction Code):ITYPE=',
     $ ITYPE,' ***'
      STOP

 9999 WRITE(*,*)
     $ 'Pi Interaction mode mismatch:IPAR=',IPAR,':ITYPE=',ITYPE
      RSSGCRS=0.0E0
      RETURN
      END
      

      
