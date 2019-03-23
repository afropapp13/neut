***********************************************************************
*     ------------------------------
      SUBROUTINE RNSWAV(PFERM)
*     ------------------------------
*     
*     ( purpose )
*       random Fermi momentum of S-wave
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
*       2001. 2. 9 ; RNDM -> RLU 
*
*     ( comment )
*
***********************************************************************
      DATA DPMAX/111./
      DATA PMAX/225./
C
C --  SET MAX VALUE
C
      CALL VCPSWAVE(DPMAX,DISMAX)
C
C --  SET PFERM
C
 10   CONTINUE
C      P=PMAX*RNDM(DUMMY)
      P=PMAX*RLU(IDUMMY)
      CALL VCPSWAVE(P,HIT)
C      IF (RNDM(DUMMY)*DISMAX.GT.HIT) GOTO 10
      IF (RLU(IDUMMY)*DISMAX.GT.HIT) GOTO 10

      PFERM=P

      RETURN
      END
