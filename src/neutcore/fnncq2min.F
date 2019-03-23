************************************************************************
*     ------------------------
      FUNCTION FNNCQ2MIN(E,INUC)
*     ------------------------
*
*     (Purpose)
*       Calc. minimum value of Q**2
*
*     (Input)
*       E      : NEUTRINO ENERGY ( GEV )
*       INUC   : Target TYPE
*                  2212 : Proton
*                  2112 : Neutron
*     
*     (Output)
*       FNNCQ2MIN : minimum value of Q**2
*
*     (Creation Date and Author)
*       2011.04.03 ; Y.Hayato
*     
************************************************************************

      IMPLICIT NONE

      REAL*4    FNNCQ2MIN
      REAL*4    E
      INTEGER*4 INUC

      REAL*4 AM,PM
      REAL*4 A,B,C,TESTV
C     
      AM=0
      CALL MCMASSGV(INUC,PM)
C     
      A=1.+2.*E/PM
      B=2.*E**2-(1.+E/PM)*AM**2
      C=AM**4
C**********************************
C     TEST
C**********************************
      TESTV=B**2-A*C
      IF(TESTV.LT.0.) THEN
         PRINT *,' TESTV< 0 IN Q2MIN ********'
         FNNCQ2MIN=0.
         RETURN
      END IF
C**********************************
      FNNCQ2MIN=(B-SQRT(B**2-A*C))/A
      RETURN
      END
