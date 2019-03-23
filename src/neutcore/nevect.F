**********************************************************************
*     ------------------------------
      SUBROUTINE NEVECT(IP,POS,PMOM,IERR)
*     ------------------------------
*     
*     ( purpose )
*       VECTOR GENERATION ( NEUTRINO INTERACTION )
*
*     ( input ) 
*       IP      : PARICLE CODE
*       POS(3)  : VERTEX POSITION
*       PMOM(3) : MOMENTUM ( MEV/C)
*
*     ( output ) 
*       COMMON VCWORK,VCVRTX
*
*     ( creation date and author )
*       1995.02.25 ; First version by K.Kaneyuki
*       1996.01.11 ; K.Kaneyuki  NENEFO16 = 1 -> no nuclear effect in O16
*       1998.03.07 ; add nek0(k0 -> KL or KS) by J.K.
*       2004.07.06 ; save interaction position in O16 by M.Ishitsuka
*       2007.11.05 ; G.Mitsuka de-excited gamma is considered only fot O16
*                              this treatment is tentatively
*       2009.12.09 ; R.Tacik, P. de Perio
*                    - Add nucleon ejection by pi/Delta absorption
*
**********************************************************************
      IMPLICIT NONE

#include "mcgenpar.h"
#include "necard.h"
#include "fsihist.h"

      INTEGER*4 IP

      INTEGER*4 IERR

      REAL POS(3),PMOM(3)

      real*4    EVPIPROB
      external  EVPIPROB
C
      IERR = 1
C
      if (QUIET.le.1) then
         write(*,*) '-----------------------------------------------------'
      endif
C     
C --  CLEAR COMMONS
C
      CALL VCCLCM
      CALL NECLPOS
C     
C --  MAKE VECTOR
C      
      CALL NEVENT(IP,PMOM,IERR)

      IF (IERR.ne.0) THEN
         RETURN
      ENDIF
C
C --  COPY PARAMETER FROM NEWOK TO VCWORK, VCVRTX
C
      CALL NECPNEWK(POS)
C
C --  TRACE IN O16
C
      CALL EFCLFSI
      IF (NENEFO16.EQ.0) CALL NEFFECT
C
C --  EXECUTE TAU DECAY
C
      CALL TAUVECT
C
C --  CONSIDER NULCEAR EXCITATION OF GAMMA
C
      if(NUMATOM.eq.16) CALL NEMKGAM ! only 16O
c      CALL NEMKGAM
C
C --  translate K0 -> K LONG or K SHORT 
C
      CALL NEK0
C
C --  decay omega meson 
C
      CALL NEOMEGA
C
C --  nucleon rescattering
C
      call nrintr
C      
C --  nucleon ejection following pion/delta absorption
C
      call absneject
C
C -- Calculate FSIPROB after all possible instances of EFTRACE
C
      FSIPROB=EVPIPROB()
     
      RETURN
      END
