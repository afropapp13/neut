************************************************************************
*     -----------------
      subroutine necrsmode(ipar,mode)
*     -----------------
*
*     (Purpose)
*        Modify cross section scaling factors to select mode
*
*     (Input)
*        IPAR   : Neutrino type
*                   12  nu_e
*                  -12  nu_e_bar
*                   14  nu_mu
*                  -14  nu_mu_bar
*                   16  nu_tau
*                  -16  nu_tau
*        MODE   : Neut mode (0: total)
*
*     (Creation Date and Author)
*       2010.11.03 ; P. de Perio
*                    WARNING: This must be kept in sync with:
*                             - neutcore/nemodsel.F
*                             - neutcore/fntotpau.F
*       2013.04.11 ; A. Redij - Add npnh mode
************************************************************************

      IMPLICIT NONE

#include "necard.h"

      integer ipar, mode, i
      logical first


      if (mode.ne.0) then
         NEMODFLG=-1
         do i=1,30
            CRSNEUT(I) = 0
            CRSNEUTB(I)= 0
         end do

C     Neutrino
         if (ipar.gt.0) then
            if (mode.eq.1) then 
               CRSNEUT(1) = 1.
            else if (mode.eq.2)then
               CRSNEUT(28) = 1.
            else if (mode.eq.11) then
               CRSNEUT(2) = 1.
            else if (mode.eq.12) then
               CRSNEUT(3) = 1.
            else if (mode.eq.13) then
               CRSNEUT(4) = 1.
            else if (mode.eq.15) then
               CRSNEUT(29) = 1.
            else if (mode.eq.16) then
               CRSNEUT(14) = 1.
            else if (mode.eq.17) then
               CRSNEUT(25) = 1.
            else if (mode.eq.21) then
               CRSNEUT(5) = 1.
            else if (mode.eq.22) then
               CRSNEUT(16) = 1.
            else if (mode.eq.23) then
               CRSNEUT(19) = 1.
            else if (mode.eq.26) then
               CRSNEUT(23) = 1.
            else if (mode.eq.31) then
               CRSNEUT(6) = 1.
            else if (mode.eq.32) then
               CRSNEUT(7) = 1.
            else if (mode.eq.33) then
               CRSNEUT(8) = 1.
            else if (mode.eq.34) then
               CRSNEUT(9) = 1.
            else if (mode.eq.35) then
               CRSNEUT(30) = 1.
            else if (mode.eq.36) then
               CRSNEUT(15) = 1.
            else if (mode.eq.38) then
               CRSNEUT(26) = 1.
            else if (mode.eq.39) then
               CRSNEUT(27) = 1.
            else if (mode.eq.41) then
               CRSNEUT(10) = 1.
            else if (mode.eq.42) then
               CRSNEUT(17) = 1.
            else if (mode.eq.43) then
               CRSNEUT(18) = 1.
            else if (mode.eq.44) then
               CRSNEUT(20) = 1.
            else if (mode.eq.45) then
               CRSNEUT(21) = 1.
            else if (mode.eq.46) then
               CRSNEUT(24) = 1.
            else if (mode.eq.51) then
               CRSNEUT(11) = 1.
               CRSNEUT(12) = 1.
            else if (mode.eq.52) then
               CRSNEUT(13) = 1.
            else
               write(*,*) "neutmodesel: Unknown mode=",mode
               stop
            endif

C     Anti-neutrino
         else 
            if (mode.eq.1) then 
               CRSNEUTB(1) = 1.
               CRSNEUTB(11) = 1.
            else if (mode.eq.2) then
               CRSNEUTB(28)=1.
            else if (mode.eq.11) then
               CRSNEUTB(2) = 1.
            else if (mode.eq.12) then
               CRSNEUTB(3) = 1.
            else if (mode.eq.13) then
               CRSNEUTB(4) = 1.
            else if (mode.eq.15) then
               CRSNEUTB(29) = 1.
            else if (mode.eq.16) then
               CRSNEUTB(15) = 1.
            else if (mode.eq.17) then
               CRSNEUTB(25) = 1.
            else if (mode.eq.21) then
               CRSNEUTB(5) = 1.
            else if (mode.eq.22) then
               CRSNEUTB(17) = 1.
            else if (mode.eq.23) then
               CRSNEUTB(20) = 1.
            else if (mode.eq.26) then
               CRSNEUTB(23) = 1.
            else if (mode.eq.31) then
               CRSNEUTB(6) = 1.
            else if (mode.eq.32) then
               CRSNEUTB(7) = 1.
            else if (mode.eq.33) then
               CRSNEUTB(8) = 1.
            else if (mode.eq.34) then
               CRSNEUTB(9) = 1.
            else if (mode.eq.35) then
               CRSNEUTB(30) = 1.
            else if (mode.eq.36) then
               CRSNEUTB(16) = 1.
            else if (mode.eq.38) then
               CRSNEUTB(26) = 1.
            else if (mode.eq.39) then
               CRSNEUTB(27) = 1.
            else if (mode.eq.41) then
               CRSNEUTB(10) = 1.
            else if (mode.eq.42) then
               CRSNEUTB(18) = 1.
            else if (mode.eq.43) then
               CRSNEUTB(19) = 1.
            else if (mode.eq.44) then
               CRSNEUTB(21) = 1.
            else if (mode.eq.45) then
               CRSNEUTB(22) = 1.
            else if (mode.eq.46) then
               CRSNEUTB(24) = 1.
            else if (mode.eq.51) then
               CRSNEUTB(12) = 1.
               CRSNEUTB(13) = 1.
            else if (mode.eq.52) then
               CRSNEUTB(14) = 1.
            else
               write(*,*) "neutmodesel: Unknown mode=",mode
               stop
            endif
         endif
      else
         NEMODFLG=0
         do i=1,30
            CRSNEUT(I) = 1.
            CRSNEUTB(I)= 1.
         end do
      endif

      END

