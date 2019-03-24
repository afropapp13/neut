************************************************************************
*     -----------------
      SUBROUTINE NESETFGPARAMS
*     -----------------
*
*     (Purpose)
*       Set Fermi gas parameters based on IZ (atomic number) 
*       and IA ( mass number )
*
*     (Input)
*       IZ : Atomic number
*       IA : Mass number
*
*     (Output)
*
*     (Creation Date and Author)
*       2008.08.15 ; Y.Hayato 
*       2008.11.03 ; R. Tacik  modified to return ok for iz=1
C                                SUBROUTINE NESETFGPARAMS(IZ, IA)
*       2009.04.11 ; Y.Hayato  modifyed to copy IZ  or IA from common.
*                             ( NUMFREP is required to be given from card
C
*                               to support H2O or CH2                     )
**************************************************************************
      IMPLICIT NONE

      INTEGER*4 IZ,IA

#include "necard.h"
#include "neutparams.h"      

C      IF (IZ.gt.1) then
C         NUMBNDP = IZ
C         NUMBNDN = IA-IZ
C         NUMFREP = 0
C         NUMATOM = IA
C      else if (iz.eq.1) then
C         NUMBNDP = 0
C         NUMBNDN = 0
C         NUMFREP = 1
C         NUMATOM = 1
C         return
C      endif         

      IA = NUMBNDP + NUMBNDN
      if (IA.eq.0 .and. NUMFREP.eq.1 .and. NUMATOM.eq.1) then
          return
      else if (IA.ne.NUMATOM) then
         write(*,*) 'NUMATOM(',NUMATOM,') is inconsistent with NUMBNDP(',
     $              NUMBNDP,')+NUMBNDN(',NUMBNDN,')'
         stop
      endif
      IZ = NUMBNDP

      if (NUMBNDP.eq.1 .and. NUMBNDN.eq.0) then !1p
         PFSURF = 0.001
         PFMAX  = 0.001
         VNUINI =-0.001
         VNUFIN = 0.0
      else if (NUMBNDP.eq.0 .and. NUMBNDN.eq.1) then !1n
         PFSURF = 0.000
         PFMAX  = 0.000
         VNUINI =-0.001
         VNUFIN = 0.0
      else if (IZ.EQ.5) then        ! Boron
         PFSURF = 0.200
         PFMAX  = 0.200
         VNUINI =-0.020
         VNUFIN = 0.0
      else if (IZ.EQ.6) then   ! Carbon
         PFSURF = 0.217
         PFMAX  = 0.217
         VNUINI =-0.025
         VNUFIN = 0.0
      else if (IZ.EQ.7) then   ! Nitrogen
         PFSURF = 0.221
         PFMAX  = 0.221
         VNUINI =-0.026
         VNUFIN = 0.0
      else if (IZ.EQ.8) then   ! Oxygen
         PFSURF = 0.225
         PFMAX  = 0.225
         VNUINI =-0.027
         VNUFIN = 0.0
      else if (IZ.EQ.9) then   ! Flourine
         PFSURF = 0.226
         PFMAX  = 0.226
         VNUINI =-0.027
         VNUFIN = 0.0
      else if (IZ.EQ.11) then   ! Sodium
         PFSURF = 0.229
         PFMAX  = 0.229
         VNUINI =-0.027
         VNUFIN = 0.0
      else if (IZ.EQ.13) then   ! Aluminum
         PFSURF = 0.233
         PFMAX  = 0.233
         VNUINI =-0.028
         VNUFIN = 0.0
      else if (IZ.EQ.14) then   ! Silicon
         PFSURF = 0.234
         PFMAX  = 0.234
         VNUINI =-0.028
         VNUFIN = 0.0
      else if (IZ.EQ.16) then   ! Sulfur
         PFSURF = 0.235
         PFMAX  = 0.235
         VNUINI =-0.029
         VNUFIN = 0.0
      else if (IZ.EQ.17) then   ! Clorine
         PFSURF = 0.236
         PFMAX  = 0.236
         VNUINI =-0.029
         VNUFIN = 0.0
      else if (IZ.EQ.18) then  ! Argon
         PFSURF = 0.237
         PFMAX  = 0.237
         VNUINI =-0.030
         VNUFIN = 0.0
      else if (IZ.EQ.19) then  ! Potassium
         PFSURF = 0.239
         PFMAX  = 0.239
         VNUINI =-0.030
         VNUFIN = 0.0
      else if (IZ.EQ.20) then  ! Calcium
         PFSURF = 0.241
         PFMAX  = 0.241
         VNUINI =-0.031
         VNUFIN = 0.0
      else if (IZ.EQ.22) then  ! Titanium
         PFSURF = 0.243
         PFMAX  = 0.243
         VNUINI =-0.031
         VNUFIN = 0.0
      else if (IZ.EQ.25) then  ! Manganese
         PFSURF = 0.250
         PFMAX  = 0.250
         VNUINI =-0.033
         VNUFIN = 0.0
      else if (IZ.EQ.26) then  ! Iron
         PFSURF = 0.250
         PFMAX  = 0.250
         VNUINI =-0.033
         VNUFIN = 0.0
      else if (IZ.EQ.27) then  ! Cobalt
         PFSURF = 0.250
         PFMAX  = 0.250
         VNUINI =-0.035
         VNUFIN = 0.0
      else if (IZ.EQ.28) then  ! Nickel
         PFSURF = 0.250
         PFMAX  = 0.250
         VNUINI =-0.035
         VNUFIN = 0.0
      else if (IZ.EQ.29) then  ! Copper
         PFSURF = 0.250
         PFMAX  = 0.250
         VNUINI =-0.035
         VNUFIN = 0.0
      else if (IZ.EQ.30) then  ! Zinc
         PFSURF = 0.250
         PFMAX  = 0.250
         VNUINI =-0.035
         VNUFIN = 0.0
      else if (IZ.EQ.38) then  ! Strontium
         PFSURF = 0.250
         PFMAX  = 0.250
         VNUINI =-0.035
         VNUFIN = 0.0
      else if (IZ.EQ.40) then  ! Zirconium
         PFSURF = 0.250
         PFMAX  = 0.250
         VNUINI =-0.035
         VNUFIN = 0.0
      else if (IZ.EQ.41) then  ! Niobium
         PFSURF = 0.250
         PFMAX  = 0.250
         VNUINI =-0.035
         VNUFIN = 0.0
      else if (IZ.EQ.45) then  ! Rhodium
         PFSURF = 0.250
         PFMAX  = 0.250
         VNUINI =-0.035
         VNUFIN = 0.0
      else if (IZ.EQ.50) then  ! Tin
         PFSURF = 0.250
         PFMAX  = 0.250
         VNUINI =-0.035
         VNUFIN = 0.0
      else if (IZ.EQ.67) then  ! Holmium
         PFSURF = 0.250
         PFMAX  = 0.250
         VNUINI =-0.035
         VNUFIN = 0.
      else if (IZ.EQ.79) then  ! Gold
         PFSURF = 0.250
         PFMAX  = 0.250
         VNUINI =-0.035
         VNUFIN = 0.0
       else if (IZ.EQ.82) then  ! Lead
         PFSURF = 0.245
         PFMAX  = 0.245
         VNUINI =-0.044
         VNUFIN = 0.0
       else if (IZ.EQ.83) then  ! Bismuth
         PFSURF = 0.245
         PFMAX  = 0.245
         VNUINI =-0.044
         VNUFIN = 0.0
      else if (NUMATOM .EQ. 2 .OR. NUMATOM .EQ. 3 .OR.
     $        NUMATOM .EQ. 4 .OR. NUMATOM .EQ. 20) then
         write(*,*) "Parameters not set in nesetfgparams.F 
     $- only valid for Effective SF"
      else
         print*, 'This target, IZ=',IZ,' is not supported yet, see nesetfgparams.F'
         stop
      end if

      call nesettarg
      call nrsettarg

      end
