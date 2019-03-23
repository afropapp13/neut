************************************************************************
*     -----------------------
      FUNCTION RNNCELQ2(IPAR,E,ITARG)
*     -----------------------
*
*     (Purpose)
*       Give random Q**2 for elastic scattering
*
*     (Input)
*       IPAR   : NEUTRINO TYPE
*                  12 : NEU-E
*                 -12 : NEU-E-BAR            
*                  14 : NEU-MU
*                 -14 : NEU-MU-BAR            
*                  16 : NEU-TAU
*                 -16 : NEU-TAU-BAR
*       E      : NEUTRINO ENERGY ( GEV )
*     
*     (Output)
*       RNNCELQ2 : random Q**2
*
*     (Creation Date and Author)
*       19??.??.?? ; M.NAKAHATA
*       1995.02.01 ; K. KANEYUKI  FOR S.K.
*                    RNELQ2 -> RNELQ2 ( SAME )
*       1995.02.03 ; K. KANEYUKI  BUG FIX
*     
************************************************************************
      IMPLICIT NONE

      REAL*4    RNNCELQ2
      REAL*4    E
      INTEGER*4 IPAR,ITARG

      REAL*4    FNNCQ2MAX, RLU
      REAL*8    DNNCELSQ2
      EXTERNAL  FNNCQ2MAX, DNNCELSQ2, RLU

      REAL*8    DBQ2, DBE, SIGMA, BLAST, BMAX, HIT

      REAL*4    EQ2MAX, DQ2

      INTEGER*4 I, ILAST, IDUM
      
      DBE=E
      EQ2MAX=FNNCQ2MAX(E,ITARG)
      DQ2=EQ2MAX/50.
      BMAX=0.E+0
      BLAST=0.E+0
C
C --  SEARCH PEAK ROUGHLY
C      
      DO 10 I=0,50
         DBQ2=FLOAT(I)*DQ2
         SIGMA=DNNCELSQ2(DBE,IPAR,ITARG,DBQ2)
C         write(6,*) i,dbq2,sigma
         IF(SIGMA.GT.BMAX) BMAX=SIGMA
         IF (SIGMA.LT.BLAST) GOTO 20
         BLAST=SIGMA
 10   CONTINUE
 20   CONTINUE
      ILAST=I-1
C
C --  SEARCH PEAK CLOSELY
C      
      DBQ2=FLOAT(ILAST)*DQ2+DQ2*0.02
      SIGMA=DNNCELSQ2(DBE,IPAR,ITARG,DBQ2)
      IF (SIGMA.GT.BLAST) THEN
C
C --     FORWARD
C         
         BMAX=SIGMA
         BLAST=SIGMA
         DO 30 I=2,50
            DBQ2=FLOAT(ILAST)*DQ2+FLOAT(I)*DQ2*0.02
            SIGMA=DNNCELSQ2(DBE,IPAR,ITARG,DBQ2)
C           write(6,*) i,dbq2,sigma
            IF(SIGMA.GT.BMAX) BMAX=SIGMA
            IF (SIGMA.LT.BLAST) GOTO 100
            BLAST=SIGMA
 30      CONTINUE
      ELSE IF (ILAST.NE.0) THEN
C
C --     BACKWORD
C         
         DO 40 I=1,50
            DBQ2=FLOAT(ILAST)*DQ2-FLOAT(I)*DQ2*0.02
            SIGMA=DNNCELSQ2(DBE,IPAR,ITARG,DBQ2)
C           write(6,*) i,dbq2,sigma
            IF(SIGMA.GT.BMAX) BMAX=SIGMA
            IF (SIGMA.LT.BLAST) GOTO 100
            BLAST=SIGMA
 40      CONTINUE
      ENDIF
 100  CONTINUE
C
C --  THROW DICE
C      
 200  DBQ2=RLU(IDUM)*EQ2MAX
      HIT=RLU(IDUM)*BMAX
      SIGMA=DNNCELSQ2(DBE,IPAR,ITARG,DBQ2)
      IF(SIGMA.LE.HIT)GO TO 200
      RNNCELQ2=DBQ2
      RETURN
      END
