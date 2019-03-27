*     -----------------------------------------------------
      FUNCTION EFRHO(R)
*     -----------------------------------------------------
*
*     ( purpose )
*       FUNCTION TO GIVE NUCLEAR DENSITY*R**2
*
*     ( input )
*         R     : radius(fm)
*
*     ( output )
*       efrho   : nuclear density*R**2
C

C
C     ++
C
#include<efpion.h>
C
      X=(R-C)/AF
      EFRHO=1./(1.+EXP(X))
      EFRHO=EFRHO*R**2
      RETURN
      END
