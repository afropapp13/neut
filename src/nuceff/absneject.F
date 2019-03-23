**********************************************************************
*     ------------------
      SUBROUTINE absneject
*     ------------------
*     
*     ( purpose )
*       Eject nucleons after absorption of Delta or pion
*
*     ( input ) 
*       COMMON VCWORK
*
*     ( output ) 
*       COMMON VCWORK
*
*     ( creation date and author )
*        R.Tacik/P.de Perio
*        
**********************************************************************
      IMPLICIT NONE

#include "necard.h"
#include "vcwork.h"

      INTEGER*4 i,j, ipi(50), npi, maxAbs

C Check whether use or not
      if (neabspiemit.eq.0) then
         return
      endif

C     DELTA ABSORPTION (assuming 4th slot only)
      if ((ipvc(4).eq.2224).or.(ipvc(4).eq.2214).or.
     &     (ipvc(4).eq.2114).or.(ipvc(4).eq.1114)) then
         call absdelta
      end if	

C     Note: ipi(50) must account for largest nucleus (so far 208Pb)
      maxAbs = NUMATOM/5

C     PION ABSORPTION (assuming >5th slot)
      npi = 0
      do i = 5, nvc
         if ( ((ipvc(i).eq.211).or.(ipvc(i).eq.-211).or.
     &	      (ipvc(i).eq.111)).and.
     &        (icrnvc(i).eq.0).and.
     &        (iflgvc(i).eq.3) ) then
            
            if (npi.LE.maxAbs) then
               npi = npi + 1
C              Store slot of absorbed pion
               ipi(npi) = i
            else if (npi.GT.maxAbs) then
               write(*,*) "absneject, WARNING: >",maxAbs,"PI ABSORBED"
            end if

         end if
      end do
      
      
C     Only absorb up to maxAbs pions:
C     Assume chance of absorbing more than maxAbs on one nucleus
C     is negligible and/or nucleus is blown apart.
      if (npi.ge.1) then
         do i = 1, npi

C           Erase nucleons ejected via pion intranuclear 
C           scattering since this is included in data used 
C           for the multiplicity calculations in abspi()
            do j = 6, nvc
               if ( iorgvc(j).eq.ipi(i) .and. 
     &              (ipvc(j).eq.2112 .or. ipvc(j).eq.2212) ) then
                  icrnvc(j) = 0
                  iflgvc(j) = -3
               end if
            end do


            if (i.le.maxAbs) then
               call abspi(ipi(i))
            end if
         end do
      end if
      
      RETURN
      END
