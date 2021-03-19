************************************************************************
*     -----------------------
      FUNCTION RNELQ2(IPAR,E)
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
*       RNELQ2 : random Q**2
*
*     (Creation Date and Author)
*       19??.??.?? ; M.NAKAHATA
*       1995.02.01 ; K. KANEYUKI  FOR S.K.
*                    RNELQ2 -> RNELQ2 ( SAME )
*       1995.02.03 ; K. KANEYUKI  BUG FIX
*     
************************************************************************
      EXTERNAL FNQ2MAX, DNELSQ2, FNQ2MIN
      REAL*8   DNELSQ2, DBQ2, DBE
      REAL EQ2MIN, EQ2MAX, DQ2
      REAL BMAX, BLAST, SIGMA
      INTEGER I

      DBE=E
      EQ2MIN=FNQ2MIN(E,IPAR)
      EQ2MAX=FNQ2MAX(E,IPAR)
      DQ2=(EQ2MAX-EQ2MIN)/50.
      BMAX=0.E+0
      BLAST=0.E+0
C
C --  SEARCH PEAK ROUGHLY
C      
C Let's start at the minimum Q2 instead of 0 (eq2min)
C
      DO 10 I=0,50
         DBQ2=EQ2MIN+FLOAT(I)*DQ2
         SIGMA=DNELSQ2(DBE,IPAR,DBQ2)
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
      DBQ2=EQ2MIN+FLOAT(ILAST)*DQ2+DQ2*0.02
      SIGMA=DNELSQ2(DBE,IPAR,DBQ2)
      IF (SIGMA.GT.BLAST) THEN
C
C --     FORWARD
C         
         BMAX=SIGMA
         BLAST=SIGMA
         DO 30 I=2,50
            DBQ2=EQ2MIN+FLOAT(ILAST)*DQ2+FLOAT(I)*DQ2*0.02
            SIGMA=DNELSQ2(DBE,IPAR,DBQ2)
C           write(6,*) i,dbq2,sigma
            IF(SIGMA.GT.BMAX) BMAX=SIGMA
            IF (SIGMA.LT.BLAST) GOTO 100
            BLAST=SIGMA
 30      CONTINUE
      ELSE IF (ILAST.NE.0) THEN
C
C --     BACKWARD
C         
         DO 40 I=1,50
            DBQ2=EQ2MIN+FLOAT(ILAST)*DQ2-FLOAT(I)*DQ2*0.02
            SIGMA=DNELSQ2(DBE,IPAR,DBQ2)
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
 200  DBQ2=EQ2MIN+RLU(DUM)*(EQ2MAX-EQ2MIN)
      HIT=RLU(DUM)*BMAX
      SIGMA=DNELSQ2(DBE,IPAR,DBQ2)
      IF(SIGMA.LE.HIT)GO TO 200
      RNELQ2=DBQ2
      RETURN
      END
