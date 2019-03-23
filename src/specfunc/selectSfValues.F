***********************************************************************
*  subroutine that selects p,E values from the Sf at random  **********
***********************************************************************

      subroutine selectsfValues(p)
        implicit none

c ----- include spectral function arrays, variables etc
#include "specfunc.h"
#include "../neutcore/necard.h"
#include "../neutcore/neutmodel.h"
        logical RFG
        double precision theta, phi, mag, y
        double precision val, maxValueE, maxValuep, interpVal
        double precision p(4)
        double precision modMomentum
        real dummy, RLU

        external modMomentum, RLU

        integer Pbin, Ebin

        integer nuc

c ----- work out which nucleus to choose from, and set the fermiMomentum accordingly
        if (NUMBNDP .EQ. 6) then
          nuc = 1
          fermiMomentum = 209
        else if (NUMBNDP .EQ. 8 )then
          nuc = 2
          fermiMomentum = 209
        else if (NUMBNDP .EQ. 26 ) then
          nuc = 3
          fermiMomentum = 209
        endif

c ----- if fermiMomentum set in card, then use this value (default set to -1)
        if (PFSF .ge. 0) then
          fermiMomentum = PFSF
        endif


        theta = acos(dble(RLU(dummy))*2 - 1)
        phi = 2 * 3.14159265358 * dble(RLU(dummy))


c RFG implementation for comparison
c        if (RFG .EQV. .TRUE.) then
c       select y at random, then mag follows p^2 probability
c          y = rand(0)
c          mag = 220 * y**(1./3.)
c          ETilde = 34 + 939.565 - sqrt(939.565**2 + mag**2)
c
c          
c          p(1) = 0
c          p(2) = mag*sin(theta)*cos(phi)
c          p(3) = mag*sin(theta)*sin(phi)
c          p(4) = mag*cos(theta)
c          return
c        endif

****** ----- SF selection starts here ----- ******
        maxValuep = 0
        maxValueE = 0
c ----- find highest point of function to throw against in 1D

        do pbin = 1, nBinsP(nuc)
          maxValuep = max(maxValuep, sFArray1D(nuc,pbin))
        enddo

c ----- pick values of momentum and see if they fall under the 1D curve
c ----- goto 11 occurs for failure to select a momentum value
11      continue

        mag = RLU(dummy) * p_Array(nuc,nBinsP(nuc))
        val = RLU(dummy) * maxValuep

c        if (mag .GT. 200) then
c          val = 0.5 * val
c        endif
        
        pBin = (mag) / (p_Array(nuc,2) - p_Array(nuc,1))
        if (pBin .GE. 200) then
          pBin = 199
        endif
        if (pBin .EQ. 0) then
          goto 11
        endif

c ----- interpolate between p-bins to get probability of momentum p
        interpVal = sFArray1D(nuc,pBin) +
     &  (mag - p_Array(nuc,pBin))*(sFArray1D(nuc,pBin+1)-
     & sFArray1D(nuc,pBin))/(p_Array(nuc,pBin+1)-p_Array(nuc,pBin))

        if (val .GT. interpVal) then
c ------- momentum value not selected
          goto 11
        endif
        
        do Ebin = 1, nBinsE(nuc)
          maxValueE = max(maxValueE, sFArray(nuc,pbin,Ebin) )
          maxValueE = max(maxValueE, sFArray(nuc,pbin+1,Ebin) )
        enddo

12      continue
c ----- pick Energies at random and see if they fall under the sf surface for p = mag
c ----- goto 12 occurs if the selected point is above the sf surface


        ETilde = RLU(dummy) * E_Array(nuc,nBinsE(nuc))
        val = RLU(dummy) * maxValueE

c        if (ETilde .GT. 100) then
c          val = 0.5 * val
c        endif

        EBin = (ETilde) / (E_Array(nuc,2) - E_Array(nuc,1))

        if (EBin .GE. 200) then
          EBin = 199
        endif
        if (EBin .EQ. 0) then
          goto 12
        endif
      
c interpolate between bins
      interpVal =
     +(p_Array(nuc,pBin+1)-mag)*(E_Array(nuc,EBin+1)-ETilde)*
     + sfArray(nuc,pBin, EBin)
     ++
     +(mag-p_Array(nuc,pBin))*(E_Array(nuc,EBin+1)-ETilde)*
     + sfArray(nuc,pBin+1, EBin)
     ++
     +(p_Array(nuc,pBin+1)-mag)*(ETilde-E_Array(nuc,EBin))*
     + sfArray(nuc,pBin, EBin+1)
     ++
     +(mag-p_Array(nuc,pBin))*(ETilde-E_Array(nuc,EBin))*
     + sfArray(nuc,pBin+1, EBin+1)
      
      interpVal = 
     +interpVal/((p_Array(nuc,2)-p_Array(nuc,1))*
     + (E_Array(nuc,2)-E_Array(nuc,1)))


c check if value falls above or below maxAllowed - if above, reselect
        if( val .GT. interpVal ) then
c ------- energy value not selected
          goto 12
        endif
        
        p(1) = 0
        p(2) = mag*sin(theta)*cos(phi)
        p(3) = mag*sin(theta)*sin(phi)
        p(4) = mag*cos(theta)

      return
      end
