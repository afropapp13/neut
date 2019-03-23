***********************************************************************
*  subroutine that selects p,E values from the Sf at random  **********
***********************************************************************

      subroutine selectsfValues(p)
        implicit none

c ----- include spectral function arrays, variables etc
        include "specfunc.h"
        logical RFG
        double precision theta, phi, mag, y
        double precision val, maxvalue, interpVal
        double precision p(4)
        double precision modMomentum
        real dummy, RLU

        external modMomentum, RLU

        integer Pbin, Ebin
        
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
c         p(2) = mag*sin(theta)*cos(phi)
c          p(3) = mag*sin(theta)*sin(phi)
c          p(4) = mag*cos(theta)
c          return
c        endif

****** ----- SF selection starts here ----- ******

c ----- find highest point of function to throw against
        maxvalue = 0

        do pbin = 1, nBinsP
        do Ebin = 1, nBinsE
          maxvalue = max(maxvalue, sFArray(pbin,Ebin) )
        enddo
        enddo

c ----- pick points at random and see if they fall under the sf surface
c ----- goto 11 occurs if the selected point is above the sf surface
11      continue
        
        mag = rand(0) * p_Array(nBinsP)
        ETilde = rand(0) * E_Array(nBinsE)
        val = rand(0) * maxvalue

        pBin = (mag + p_Array(1)) / (p_Array(2) - p_Array(1))

        EBin = (ETilde + E_Array(1)) / (E_Array(2) - E_Array(1))

        if (EBin .GE. 200) then
          EBin = 199
        endif
        if (pBin .GE. 200) then
          pBin = 199
        endif

      
c interpolate between bins
      interpVal =
     +(p_Array(pBin+1)-mag)*(E_Array(EBin+1)-ETilde)*sfArray(pBin, EBin)
     ++
     +(mag-p_Array(pBin))*(E_Array(EBin+1)-ETilde)*sfArray(pBin+1, EBin)
     ++
     +(p_Array(pBin+1)-mag)*(ETilde-E_Array(EBin))*sfArray(pBin, EBin+1)
     ++
     +(mag-p_Array(pBin))*(ETilde-E_Array(EBin))*sfArray(pBin+1, EBin+1)
      
      interpVal = 
     +interpVal/((p_Array(2)-p_Array(1))*(E_Array(2)-E_Array(1)))


c check if value falls above or below maxAllowed - if above, reselect
        if( val .GT. interpVal ) then
          goto 11
        endif
        
        p(1) = 0
        p(2) = mag*sin(theta)*cos(phi)
        p(3) = mag*sin(theta)*sin(phi)
        p(4) = mag*cos(theta)

      return
      end
