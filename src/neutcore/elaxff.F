************************************************************************
*     -------------------------
      FUNCTION ELAXFF(Q2, nc, sgn)
*     -------------------------
*
*     (Purpose)
*       Return axial form factor value given chosen form
*
*     (Input)
*       Q2     : FOUR-MOM TRANSFER (GeV^2)
*       nc     : neutral current flag
*                 true  : NC
*                 false : CC
*       sgn    : PDG Sign (Neutrino/Antineutrino)
*                +1 : Neutrino
*                -1 : Anti-neutrino
*     
*     (Output)
*       ELAXFF : Elastic Form Factor (Dimensionless)
*
*     (Creation Date and Author)
*     2016.09.09 ; P. Stowell
*     
*     To add a new fa form, just add a cc file copying the fa2COMP format.
*     Then enable it here, as well as inside qedifcrs.c if you want RW
************************************************************************
      IMPLICIT NONE
#include "../neutcore/neutmodel.h"

      LOGICAL nc
      INTEGER sgn
      
      double precision Q2
      
C     Read in FA Functions
      double precision fabbba07
      external fabbba07

      double precision fa2COMP
      external fa2COMP

      double precision fa3COMP
      external fa3COMP
      
      double precision faZEXP
      external faZEXP
      
     
C     Setup Variables      
      REAL*8 ELAXFF
      REAL*8 FA
      REAL*4 gA
      REAL*8 MA
      
C      Dipole FA
      if (MDLQEAF.eq.1 .or. MDLQEAF.eq.0) then
         gA = -1.2673
C         MA = 1000 * dble(XMAQE)                 -2017/03/28 Y. Hyaato
C         FA = gA / (( 1 -( Q2 /( MA**2 ) ) )**2) -2017/03/28 Y. Hyaato
         if (nc .EQV. .TRUE.) then
            FA = gA / (( 1 +( Q2 /( XMANCEL**2 ) ) )**2)
         else
            FA = gA / (( 1 +( Q2 /( XMAQE**2 ) ) )**2)
         endif

C     BBBA07 FF         
      else if (MDLQEAF.eq.2) then
         FA = fabbba07( Q2 )
C     2 Component FF
      else if (MDLQEAF.eq.3) then
         FA = fa2COMP( Q2 )
C     3 Component FF
      else if (MDLQEAF.eq.4) then
         FA = fa3COMP( Q2 )
C     ZEXP FF
      else if (MDLQEAF.eq.5) then
         FA = faZEXP( Q2 )
      endif

C     Apply NC Scaling to FA
      if (nc .EQV. .TRUE.) then
         FA = sgn*(0.5 * FA)
      endif

C     RETURN FA
      ELAXFF = FA
      
      RETURN
      END
