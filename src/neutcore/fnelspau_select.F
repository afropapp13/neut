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
      integer IPAR, mode
      real fnelscrs
      external fnelscrs

      real     fnelspau_intpol,fnelspau_sf,fnelspau_effsf,fnels_erpascl,
     $         fn_1p1h_n
      external fnelspau_intpol,fnelspau_sf,fnelspau_effsf,fnels_erpascl,
     $         fn_1p1h_n

C -- currently dummy argument
      if ( ipar.gt.0 ) then
         mode = 1
      else
         mode = -1
      endif

C     Free Nucleon Calculation
      if ((numbndp + numbndn).eq.1) then
         FNELSPAU_SELECT = fnelscrs(E, IPAR)
         return
      endif


C     Nuclear Calculations

      if (MDLQEAF.gt.2) then
         write(*,*) "ERROR: Axial FF other than Dipole and BBBA07"
         write(*,*) "Cannot generate on Nuclear Target"
         write(*,*) "Generte on dipole and use RW tools instead"
         stop 9
      endif
      
C
C--- Nieves model includes RPA correction ( = 2XXX )
C
      if (mod(MDLQE, 10000)/1000 .EQ. 2)then
         FNELSPAU_SELECT = fn_1p1h_n(e,ipar,mode)
         goto 100
      endif

      if ( (mod(MDLQE, 1000)/100 .EQ. 4) .and.
     &     ( (numbndp .eq. 6) .or. 
     &       (numbndp .eq. 8) .or.
     $       (numbndp .eq. 26))) then
         FNELSPAU_SELECT = fnelspau_sf(E, IPAR)
      else if (  ((mod(MDLQE, 1000)/100 .EQ. 6 ).or.
     &            (mod(MDLQE, 1000)/100 .EQ. 7 )    ).and.
     &           ((NUMATOM .eq.  2) .or. (NUMATOM .eq.  3) .or.
     &            (NUMATOM .eq.  4) .or. (NUMATOM .eq. 12) .or.
     &            (NUMATOM .eq. 16) .or. (NUMATOM .eq. 20) .or.
     &            (NUMATOM .eq. 27) .or. (NUMATOM .eq. 40) .or.
     &            (NUMATOM .eq. 56) .or. (NUMATOM .eq. 63) .or.
     &            (NUMATOM .eq. 64) .or. (NUMATOM .eq.208))) then
         FNELSPAU_SELECT = fnelspau_effsf(E, IPAR)
      else
         FNELSPAU_SELECT = fnelspau_intpol(E, IPAR)
      endif
      
      if (mod(MDLQE, 10000)/1000 .EQ. 1)then
         FNELSPAU_SELECT =FNELSPAU_SELECT*fnels_erpascl(E,IPAR)
      endif
   
 100  continue
      return
      end
