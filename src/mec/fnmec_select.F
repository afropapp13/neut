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
      integer*4 IPAR
      real*4    fnmec_nieves, fn_2p2h_ht
      external  fnmec_nieves, fn_2p2h_ht

      integer*4 mode
      
      if (MDL2P2H.eq.1) then 
         FNMEC_SELECT = fnmec_nieves(E, IPAR)/1000.
      else if (MDL2P2H.eq.2) then
         FNMEC_SELECT = fn_2p2h_ht(e,ipar,mode)
      else
         write(*,*) 'Unknown MDL2P2H value ',MDL2P2H
         stop
      endif

      return
      end
