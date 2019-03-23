**********************************************************************
*     -------------------------------------------
      SUBROUTINE VCPHSVCT(W,AMAS,NBODY,PMOM,IERR)
*     -------------------------------------------
*
*     ( purpose )
*       VECTOR GENERATION OF NBODY DECAY
*       ONLY PHASE SPACE IS CONSIDERD.
*
*     ( input )
*       W       : INVARIANT MASS
*       AMAS(10): MASS OF DECAY PARTICLES
*       NBODY   : NBODY DECAY ( P->E+ PI0 : NBODY=2 )
*
*     ( output )
*       PMOM(10): MOMENTUM OF DECAY PARTCLES
*       IERR    : ERROR CODE
*
*     ( creation date and author )
*       19??. ?. ? ; original version : phsvct@monlib
*       1989.08.03 ; K.Kaneyuki   using RAMBO
*       1992. 8. 6 ; K.Kaneyuki   mod. for S.K.
*                    PHSVCT -> VCPHSP
*       1995. 2. 7 ; change name
*                    VCPHSP -> VCPHSVCT      
*       2001. 2. 9 ; RNDM -> RLU
*
*     ( comment )
*     
*********************************************************************

      DIMENSION PMOM(3,10),AMAS(10)
C     
C     -- FOR RAMBO
C     
      REAL*8    ET,XM(10),PRMB(4,10),WT,WTMAX
      DATA ITMAX/1000/
      
      IERR=0
C     
C     -- CLEAR BUFFER
C     
      DO 10 I=1,10
         DO 20 J=1,4
            PRMB(J,I)=0.0D0
   20    CONTINUE
   10 CONTINUE
C     
C     --  PARAMTER CHECK
C     
      XMT=0.
      DO 25 I=1,NBODY
         XMT=XMT+ABS(AMAS(I))
   25 CONTINUE

      IF (XMT.GT.W) THEN
         IERR=1
         GOTO 80
      ENDIF
C     
C     -- PARAMETER OF PHSVCT => RAMBO
C     
      ET=W
      DO 30 I=1,NBODY
         XM(I)=AMAS(I)
   30 CONTINUE
C     
C     -- SET MAX WEIGHT OF RAMBO TO WTMAX
C     
      CALL VCWTRMBO(NBODY,ET,XM,WTMAX)
C     
C     -- CALCULATE PHASE SPACE
C     
      DO 35 J=1,10
         DO 40 I=1,ITMAX
            CALL VCRAMBO(NBODY,ET,XM,PRMB,WT,IERR)
            IF (IERR.NE.0) THEN
               WRITE(6,*) ' RAMBO ERROR ',IERR
               GOTO 80
            ENDIF
C            IF (WT.GE.WTMAX*RNDM(DUM)) GOTO 50
            IF (WT.GE.WTMAX*RLU(IDUM)) GOTO 50
   40    CONTINUE
         WTMAX=WTMAX*0.1
   35 CONTINUE
      WRITE(6,*) ' ITERATION ERROR AT PHSVCT : ',I
      WRITE(6,*) ' LAST RAMBO WEIGHT : ',WT,' WTMAX : ',WTMAX
      WRITE(6,*) NBODY,' BODY : ',(AMAS(I),I=1,NBODY)
      IERR=1
      GOTO 80
C
C --  COPY MOMENTUM TO PMOM
C
   50 DO 60 I=1,NBODY
         DO 70 J=1,3
            PMOM(J,I)=PRMB(J,I)
   70    CONTINUE
   60 CONTINUE

   80 CONTINUE
      
      RETURN
      END
