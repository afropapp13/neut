*     -----------------------------------------------------
      FUNCTION EFFRMGAS(DUM1,DUM2,R)
*     -----------------------------------------------------
*
*     ( purpose )
*       GET FERMI SURFACE MOMENTUM FOR NUCLEI
*
*     ( input )
*         DUM1   : # OF NUCLEONS  (DUMMY)
*         DUM2   : # OF PROTONS OR NEUTRONS  (DUMMY)
*         R      : radius(fm)
*     ( output )
*       EFFRMGAS : Fermi surface momentum

#include<efpion.h>

C     INTEGRAL ABRHO FROM 0.0 TO CC2 -> 100.18
C     (4 / 3)* PI * (R0^3 * A)       -> 115.75

C     M(PI) = 140MEV = 0.709 FM^-1
C     0.709 ** 3 = 0.356
C     0.356 * 0.48 = 0.1709

      DOUBLE PRECISION EFABRHO
      EXTERNAL EFABRHO
      REAL DENSITY

      REAL  PI
      PARAMETER (PI=3.141593)

      R0=1.2

      DENSITY = EFABRHO(R)
      DENSITY = DENSITY * 0.1709
      DUMMY=(ZZ/AN*3*PI*PI*DENSITY)**0.33333
      EFFRMGAS=DUMMY*197.33
      RETURN
      END
