************************************************************************
*     --------------------------------
      REAL FUNCTION FNELSPAU_SELECT(E,IPAR)
*     --------------------------------
*
*     (Purpose)
*       Select FNELSPAU or FNELSPAU_INTPOL following
*       the Fermi momentum set in the neut.card
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
*
*     (Output)
*       FNELSPAU_SELECT : CROSS SECTION ( 10**-38cm^2 )
*
*     (Creation Date and Author)
*       2007.11.06 ; G. Mitsuka
*       2009.01.26 ; Y. Hayato 
*                    Even if PFSURF is 0.225(GeV/c), 
*                           fnelspau_intpol is used.
C       2013.04.13 ; Y. Hayato
C                    Support spectral function function 
C       2013.12.30   A. Redij
C		     Support for RPA correction
*
************************************************************************

      implicit none
#include "neutparams.h"
#include "neutmodel.h"
#include "necard.h"

      real E
      integer IPAR
      real     fnelspau_intpol,fnelspau_sf,fnels_erpascl
      external fnelspau_intpol,fnelspau_sf,fnels_erpascl

      if ((mod(MDLQE, 1000)/400 .EQ. 1) 
     & .AND. (numbndp .eq. 6 .or. 
     &    numbndp .eq. 8 .or. numbndp .eq. 26)) then
         FNELSPAU_SELECT = fnelspau_sf(E, IPAR)
      else
         FNELSPAU_SELECT = fnelspau_intpol(E, IPAR)
      endif
      
      if (mod(MDLQE, 10000)/1000 .EQ. 1)then
         FNELSPAU_SELECT =FNELSPAU_SELECT*fnels_erpascl(E,IPAR)
      endif
   
      
      return
      end
