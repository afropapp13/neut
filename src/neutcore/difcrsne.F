************************************************************************
*     ---------------------------------
      REAL FUNCTION DIFCRSNE(E,IPAR,IC)
*     ---------------------------------
*
*     (Purpose)
*       Return cross section of Diffractive pion production
*
*     (Input)
*       E     : Incident Neutrino Energy ( GeV )
*       IPAR  : Incident Neutrino ID ( NOT USED )
*       IC    : Interaction mode
*                1 : C.C.
*                0 : N.C.
*     (Output)
*       DIFCRSNE :  Cross Section / 16O Nucleous ( 10**(-38) cm**2 )
*      
*     (Creation Date and Author)
*       2006.01.08  G.Mitsuka
*
*     (Comment)
*       Original difcrsne.F was renamed to difcrsnerein.F,
*       which returns cross section of Rein & Sehgal model.
*       Current difcrsne.F can select each models;
*       Rein & Sehgal or Kartvtsev et al.
*
************************************************************************
      implicit none

#include "necard.h"

      integer ipar,ic
      real e
      real     difcrsnerein
      external difcrsnerein

      if(NEDIFPI.eq.0) then     ! Select Rein model
         DIFCRSNE = DIFCRSNEREIN(E,IPAR,IC)
      else if (NEDIFPI.eq.1) then ! no W-cut applied, but Rein still model
         DIFCRSNE = DIFCRSNEREIN(E,IPAR,IC)
      else
         write(*,*) 'undefined model of diffractive pi ( ',
     $        nedifpi ,') was specified.'
         stop
      endif

      return
      end
