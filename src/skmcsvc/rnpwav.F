***********************************************************************
*     ------------------------
      SUBROUTINE RNPWAV(PFERM)
*     ------------------------
*     
*     ( purpose )
*       random Fermi momentum of P-wave
*
*     ( input ) 
*       none
*
*     ( output ) 
*       PFERM : Fermi momentum ( MeV/c )
*
*     ( creation date and author )
*                  ; First version by M.Nakahata
*       1992. 9.18 ; Mod by K.Kaneyuki for S.K.
*       1995. 2. 7 ; change names
*       2001. 2. 9 ; RNDM -> RLU by K.K.
*     
*     ( comment )
*
***********************************************************************
      DATA DPMAX/150./
      DATA PMAX/275./
      INTEGER IDUMMY
C
C --  SET MAX VALUE
C
      CALL VCPPWAVE(DPMAX,DISMAX)
C
C --  SET PFERM
C
 10   CONTINUE
C      P=PMAX*RNDM(DUMMY)
      P=PMAX*RLU(IDUMMY)
      CALL VCPPWAVE(P,HIT)
C      IF (RNDM(DUMMY)*DISMAX.GT.HIT) GOTO 10
      IF (RLU(IDUMMY)*DISMAX.GT.HIT) GOTO 10

      PFERM=P

      RETURN
      END
