*     -----------------------------------------------------
      FUNCTION EFGAMMNT(E)      
*     -----------------------------------------------------
*     
*     ( purpose )
*       photo single pion production total cross-section
*        (gamma + n -> p + pi-)
*     ( input ) 
*       E       : energy of gamma(MeV)
*
*     ( output ) 
*       EFGAMMNT: total cross-section (mb)


      IMPLICIT NONE
      REAL EFGAMMNT
      REAL E

      REAL CRSTBL(16),ETABLE(16)
      INTEGER I
      REAL CRSSCT

      DATA ETABLE/270.,290.,310.,330.,350.,370.,390.,410.,430.,
     $            450.,470.,515.,565.,590.,640.,665./
      DATA CRSTBL/242.7,273.4,289.8,274.5,222.8,202.7,158.9,
     $            146.4,143.4,126.7,123.3,105.1,107.8,105.2,
     $            104.2,109.0/

      IF ((E.LT.ETABLE(1)).OR.(E.GT.ETABLE(16))) THEN
         EFGAMMNT=0.0
         RETURN
      ENDIF

      DO 100 I=1,16
         IF (E.LT.ETABLE(I)) GOTO 110
 100  CONTINUE

 110  CONTINUE
      CRSSCT= CRSTBL(I-1)
     $       + (CRSTBL(I)-CRSTBL(I-1))/(ETABLE(I)-ETABLE(I-1))
     $        *(E-ETABLE(I-1))

      EFGAMMNT=CRSSCT

      RETURN
      END
      
