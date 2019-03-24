************************************************************************
*     ---------------------------------
      REAL FUNCTION COHCRSNE(E,IPAR,IC)
*     ---------------------------------
*
*     (Purpose)
*       Return cross section of Coherent pion production
*
*     (Input)
*       E     : Incident Neutrino Energy ( GeV )
*       IPAR  : Incident Neutrino ID ( NOT USED )
*       IC    : Interaction mode
*                1 : C.C.
*                0 : N.C.
*     (Output)
*       COHCRSNE :  Cross Section / 16O Nucleous ( 10**(-38) cm**2 )
*      
*     (Creation Date and Author)
*       2006.01.08  G.Mitsuka
*
*     (Comment)
*       Original cohcrsne.F was renamed to cohcrsnerein.F,
*       which returns cross section of Rein & Sehgal model.
*       Current cohcrsne.F can select each models;
*       Rein & Sehgal or Kartvtsev et al.
*
************************************************************************
      implicit none

#include "necard.h"

      integer ipar,ic
      real e
      real     cohcrsnerein,cohcrsnekart,cohcrsneberger
      external cohcrsnerein,cohcrsnekart,cohcrsneberger

      if ((numbndp+numbndn).le.1) then
         cohcrsne = 0.
         return
      endif

      if(NECOHEPI.eq.0) then     ! Select Rein & Sehgal model
         COHCRSNE = COHCRSNEREIN(E,IPAR,IC)
      elseif(NECOHEPI.eq.1) then ! Select Kartavtsev et al. model
         COHCRSNE = COHCRSNEKART(E,IPAR,IC)
      elseif(NECOHEPI.eq.2) then ! Select Berger and Sehgal model
         COHCRSNE = COHCRSNEBERGER(E,IPAR,IC)
      else
         write(*,*) 'undefined model of coherent pi ( ',
     $        cohcrsne ,') was specified.'
         stop
      endif

      return
      end
