************************************************************************
*     --------------------------------
      REAL FUNCTION FNMEC_SELECT(E,IPAR)
*     --------------------------------
*
*     (Purpose)
*       To Select  MEC model. 
*
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
*       FNMEC_SELECT : CROSS SECTION ( 10**-38cm^2 )
*
*     (Creation Date and Author)
*       2013.02.25 ; A.Redij; Derived from FNELSPAU_SELECT for npnh
*                    Here we may opt between different MEC models. 
*		     At present only Nieves 2p2h models is added.
*
************************************************************************

      implicit none
#include "neutparams.h"
#include "neutmodel.h"
      real E
      integer IPAR
      real  fnmec_nieves 
      external  fnmec_nieves

      integer iloadmecsel
      save iloadmecsel
      data iloadmecsel/0/


         if (iloadmecsel.eq.0) then
            print*, ' MEC cross section is used '    
         endif 
         
         FNMEC_SELECT = fnmec_nieves(E, IPAR)/1000.

         iloadmecsel = 1
      
      return
      end
