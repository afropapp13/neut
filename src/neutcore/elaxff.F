************************************************************************
*     -------------------------
      DOUBLE PRECISION FUNCTION ELAXFF(Q2_GeV, isNC, nuPDGSign)
*     -------------------------
*
*     (Purpose)
*       Return axial form factor value given chosen form
*
*     (Input)
*       Q2_GeV    : FOUR-MOM TRANSFER (GeV^2)
*       isNC      : neutral current flag
*                   true  : NC
*                   false : CC
*       nuPDGSign : PDG Sign (Neutrino/Antineutrino)
*                  +1 : Neutrino
*                  -1 : Anti-neutrino
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
#include "neutmodel.h"

      LOGICAL isNC
      INTEGER nuPDGSign

      double precision Q2_GeV

C     Read in FA Functions
      double precision fabbba07
      external fabbba07

      double precision fa2COMP
      external fa2COMP

      double precision fa3COMP
      external fa3COMP

      double precision faZEXP
      external faZEXP

      double precision gA
      parameter (gA = -1.2673)


C     Setup Variables
      DOUBLE PRECISION FA
      DOUBLE PRECISION MA

C      Dipole FA
      if (MDLQEAF.eq.1) then
         if (isNC .EQV. .TRUE.) then
            FA = gA / (( 1 +( Q2_GeV /( XMANCEL**2 ) ) )**2)
         else
            FA = gA / (( 1 +( Q2_GeV /( XMAQE**2 ) ) )**2)
C             WRITE(*,*) "Dipole: MAQE = ", XMAQE, FA
         endif

C     BBBA07 FF
      else if (MDLQEAF.eq.2) then
         FA = fabbba07( Q2_GeV )
C     2 Component FF
      else if (MDLQEAF.eq.3) then
         FA = fa2COMP( Q2_GeV )
C     3 Component FF
      else if (MDLQEAF.eq.4) then
         FA = fa3COMP( Q2_GeV )
C     ZEXP FF
      else if (MDLQEAF.eq.5) then
         FA = faZEXP( Q2_GeV )
      else
        write(*,*) 'INVALID MDLQEAF specified:',MDLQEAF,' should be:'
        write(*,*) '    1: Dipole'
        write(*,*) '    2: BBBA07'
        write(*,*) '    3: 2-component'
        write(*,*) '    4: 3-component'
        write(*,*) '    5: Z-expansion'
        stop
      endif

C     Apply NC Scaling to FA
      if (isNC .EQV. .TRUE.) then
         FA = nuPDGSign*(0.5 * FA)
      endif

C     RETURN FA
      ELAXFF = FA
C       WRITE(*,*) "ELAXFF: Q2_GeV", ELAXFF, Q2_GeV
C       WRITE(*,*) "ELAXFF: MDLQEAF, ISNC",  MDLQEAF, isNC

      RETURN
      END
