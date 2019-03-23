*     -----------------------------------------------------
      FUNCTION EFGAMMPT(E)
*     -----------------------------------------------------
*     
*     ( purpose )
*       photo single pion production total cross-section
*        (gamma + p -> n + pi+)
*
*     ( input ) 
*       E       : energy of gamma(MeV)
*
*     ( output ) 
*       EFGAMMPT: total cross-section (mb)
*
      
      IMPLICIT NONE
      REAL EFGAMMPT
      REAL E

      REAL CRSTBL(17),ETABLE(17)
      INTEGER I
      REAL CRSSCT

      DATA ETABLE/250,270.,290.,310.,330.,350.,370.,390.,410.,
     $            430,450.,470.,515.,565.,590.,640.,690./
      DATA CRSTBL/156.5,203.2,231.3,251.9,241.7,203.7,178.8,
     $            139.1,124.9,115.3,109.4,103.3, 93.7, 95.0, 
     $             93.6,99.7,111.7/

      IF ((E.LT.ETABLE(1)).OR.(E.GT.ETABLE(16))) THEN
         EFGAMMPT=0.0
         RETURN
      ENDIF

      DO 100 I=1,16
         IF (E.LT.ETABLE(I)) GOTO 110
 100  CONTINUE

 110  CONTINUE
      CRSSCT= CRSTBL(I-1)
     $       + (CRSTBL(I)-CRSTBL(I-1))/(ETABLE(I)-ETABLE(I-1))
     $        *(E-ETABLE(I-1))

      EFGAMMPT=CRSSCT

      RETURN
      END
      
