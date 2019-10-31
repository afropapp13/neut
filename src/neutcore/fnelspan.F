************************************************************************
*     ------------------------------
      FUNCTION FNELSPAN(E,IPAR,INUC)
*     ------------------------------
*
*     (Purpose)
*       Cross section of neutral current elastic Int. on Bounded Nucleons
*
*     (Input)
*       E      : NEUTRINO ENERGY ( GEV )
*       IPAR   : NEUTRINO TYPE
*                  12 : NEU-E
*                 -12 : NEU-E-BAR
*                  14 : NEU-MU
*                 -14 : NEU-MU-BAR
*                  16 : NEU-TAU
*                 -16 : NEU-TAU-BAR
*       INUC   : =2212 PROTON  =2112 NEUTRON
*
*     (Output)
*       FNELSPAN : CROSS SECTION ( 10**-38cm^2 )
*
*     (Creation Date and Author)
*       1988.10.08 ; T.KAJITA
*       1989.06.10 ; T.KAJITA  SMALL MODIFICATION
*         ***** AT PRESENT THE CROSS SECTIONS ARE NOT EXACT.
*         ***** CROSS SECTIONS ARE SIMPLY SCALED FROM THAT OF CC ELASTIC
*         ***** FOR THESE SCALING DATA FROM K.ABE ET AL., P.R.L.56(1986)
*               1107 ARE USED
*         ***** SIGMA(NEU    N)=1.5*SIGMA(NEU    P)           (10-JUN-89)
*               SIGMA(NEUBAR N)=1.0*SIGMA(NEUBAR P)
*               ARE USED FROM C.H.ALBRIGHT ET AL.,P.R.D14(1976)1780
*       1995.02.16 ; K. KANEYUKI  FOR S.K.
*                    ELSPAN -> FNELSPAN
*
*       2007.11.06 ; G.Mitsuka, FNELSPAU->FNELSPAU_SELECT
*       2011.03.09 ; Y.Hayato , Use FNELSPAN_SELECT
*
************************************************************************
      IMPLICIT NONE

#include "neutmodel.h"
#include "necard.h"

      REAL*4   FNELSPAN, E
      INTEGER  IPAR, INUC

      REAL*4   FNELSPAU_SELECT, FNELSPAN_INTPOL, FNELSPAN_SF,
     $     FNELSPAN_EFFSF
      EXTERNAL FNELSPAU_SELECT, FNELSPAN_INTPOL, FNELSPAN_SF,
     $     FNELSPAN_EFFSF

C -- for SF nuclei tests
      logical HaveSpecFuncForNuclei, HaveEffSpecFuncForNuclei
      external HaveSpecFuncForNuclei, HaveEffSpecFuncForNuclei

C
C ++ CROSS SECTION RATIO OF NC AND CC (NC/CC) FOR NEUTRINO AND ANTI-NEU
C
      REAL*4 RNEUNC, RNEBNC

      DATA RNEUNC/0.153/
      DATA RNEBNC/0.218/

      INTEGER  IPTEMP, NCELMODE
      REAL     CRS
C
      FNELSPAN=0.

      if ((numbndp.eq.0).and.(numbndn.eq.0)) then
         return
      endif

      IF (INUC.NE.2212 .AND. INUC.NE.2112) RETURN
      IF (IABS(IPAR).NE.12 .AND. IABS(IPAR).NE.14 .AND.
     &   IABS(IPAR).NE.16) RETURN

      IF (IPAR.GT.0) IPTEMP=12
      IF (IPAR.LT.0) IPTEMP=-12

      NCELMODE = mod(MDLQE,100)/10
      if (mod(MDLQE, 1000)/100 .EQ. 4 ) then
        NCELMODE = mod(mod(MDLQE, 1000),100)/10
      else if (mod(MDLQE, 1000)/100 .EQ. 7
     $       .OR. mod(MDLQE, 1000)/100 .EQ. 6) then
        NCELMODE = mod(mod(MDLQE, 1000),100)/10
      endif

      if ((mod(MDLQE, 1000)/100 .EQ. 4) .and.
     &     HaveSpecFuncForNuclei() ) then
         FNELSPAN=FNELSPAN_SF(E,IPTEMP,INUC)
      else if ( ((mod(MDLQE, 1000)/100 .EQ. 6 ).or.
     &           (mod(MDLQE, 1000)/100 .EQ. 7 )     ).and.
     &           HaveEffSpecFuncForNuclei() ) then
         FNELSPAN=FNELSPAN_EFFSF(E,IPTEMP,INUC)
      else if (NCELMODE .eq. 0) then
         CRS=FNELSPAU_SELECT(E,IPTEMP)

         IF (IPAR.GT.0)  THEN
            FNELSPAN=CRS*RNEUNC
            IF(INUC.EQ.2112) FNELSPAN=FNELSPAN*1.5
         END IF
         IF (IPAR.LT.0) FNELSPAN=CRS*RNEBNC
         return
      else if (NCELMODE .ge. 1) then
         FNELSPAN=FNELSPAN_INTPOL(E,IPTEMP,INUC)
      endif

      RETURN
      END
