*     -----------------------------------------------------
      FUNCTION EFABRHO(R)      
*     -----------------------------------------------------
*     
*     ( purpose )
*       FUNCTION TO GIVE NUCLEAR DENSITY (UNIT : AT CENTER OF NULEI -- 1)
*
*     ( input ) 
*       r       : radius(fm)
*
*     ( output ) 
*       efabrho : normalized nucelar density
*
*     (Creation Date and Author)
*     2007.11.05 ; G.Mitsuka - change EFABRHO formula to 
*                  3-parameter fermi model
C

#include<efpion.h>

      X=(R-C)/AF
C      EFABRHO=1./(1.+EXP(X))/DFACT
      EFABRHO=(1.+WPARM*(R/C)**2)/(1.+EXP(X))/DFACT
      IF(EFABRHO.GE.1.) EFABRHO=0.9999999
      IF(EFABRHO.LE.0.) EFABRHO=1.E-5
      RETURN
      END
