***********************************************************************
*     ------------------------------
      SUBROUTINE RNDIR(DIR)
*     ------------------------------
*     
*     ( purpose )
*       generate random direction 
*
*     ( input ) 
*       none
*
*     ( output ) 
*       DIR(3) : random direction
*
*     ( creation date and author )
*       1980. 9.14 ; First version by K.Arisaka
*       1992. 9.18 ; Mod by K.Kaneyuki for S.K.
*       2001. 2. 9 ; RNDM -> RLU K.K.
*
*     ( comment )
*
***********************************************************************
      REAL DIR(3)
      INTEGER IDUMMY

    1 SS=0.
      DO 10 I=1,3
C         DIR(I)=2.*RNDM(DUMMY)-1.
         DIR(I)=2.*RLU(IDUMMY)-1.
         SS=SS+DIR(I)**2
   10 CONTINUE
      S=SQRT(SS)
      IF(S.GT.1.) GO TO 1

      DO 20 I=1,3
         DIR(I)=DIR(I)/S
 20   CONTINUE

      RETURN
      END
