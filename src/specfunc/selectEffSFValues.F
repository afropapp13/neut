***********************************************************************
* subroutine that selects p values from the Effective SF distribution *
***********************************************************************

      subroutine selectEffSFValues(p)
        implicit none

c ----- include spectral function arrays, variables etc
#include "effspecfunc.h"
#include "necard.h"
#include "neutmodel.h"
        logical RFG
        double precision theta, phi, mag, y
        double precision val, maxValueE, maxValuep, interpVal
        double precision p(4)
        double precision modMomentum
        real dummy, RLU

        external modMomentum, RLU

        integer Pbin 

        integer nuc

c ----- work out which nucleus to choose from, and set the fermiMomentum accordingly
C     CW: I think I need to ignore Fermi momentum for the moment
        fermiMomentum = 0

        if (NUMATOM .EQ. 2) then
          nuc = 1
          fermiMomentum = 100
          Etilde = 0.13
        else if (NUMATOM .EQ. 3 )then
          nuc = 2
          fermiMomentum = 115
          Etilde = 5.3
        else if (NUMATOM .EQ. 4 )then
          nuc = 3
          fermiMomentum = 190
          Etilde = 14.0
        else if (NUMATOM .EQ. 12 )then
          nuc = 4
          fermiMomentum = 228
          Etilde = 12.5
        else if (NUMATOM .EQ. 16 )then
          nuc = 5
          fermiMomentum = 228
          Etilde = 12.5
        else if (NUMATOM .EQ. 20 )then
          nuc = 6
          fermiMomentum = 230
          Etilde = 16.6        
       else if (NUMATOM .EQ. 27 )then
          nuc = 7
          fermiMomentum = 236
          Etilde = 12.5
       else if (NUMATOM .EQ. 40 )then
          nuc = 8
          fermiMomentum = 241
          Etilde = 20.6
       else if (NUMATOM .EQ. 56 )then
          nuc = 9
          fermiMomentum = 241
          Etilde = 15.1
       else if (NUMATOM .EQ. 63 )then
          nuc = 10
          fermiMomentum = 245
          Etilde = 18.8
       else if (NUMATOM .EQ. 64 )then
          nuc = 11
          fermiMomentum = 245
          Etilde = 18.8
       else if (NUMATOM .EQ. 208 )then
          nuc = 12
          fermiMomentum = 245
          Etilde = 18.8
       endif

C     CW: randomly choose direction
        theta = acos(dble(RLU(dummy))*2 - 1)
        phi = 2 * 3.14159265358 * dble(RLU(dummy))

****** ----- SF selection starts here ----- ******
        maxValuep = 0
c ----- find highest point of function to throw against in 1D

        do pbin = 1, nBinsP(nuc)
          maxValuep = max(maxValuep, SFArray(nuc,pbin))
        enddo

c ----- pick values of momentum and see if they fall under the 1D curve
c ----- goto 11 occurs for failure to select a momentum value
11      continue

C     CW: only need to pick this value, remember momentum < 0.65 GeV by construction.
        mag = RLU(dummy) * p_Array(nuc, nBinsP(nuc))
        val = RLU(dummy) * maxValuep

        pBin = (mag) / (p_Array(nuc,2) - p_Array(nuc,1))
        if (pBin .GE. 651) then
          pBin = 650
        endif
        if (pBin .EQ. 0) then
          goto 11
        endif

c ----- interpolate between p-bins to get probability of momentum p
C     CW: only use this for the moment
C     CW: is this interpolation correct?
        interpVal = SFArray(nuc,pBin) +
     &       (mag - p_Array(nuc,pBin))*(SFArray(nuc,pBin+1)-
     &       SFArray(nuc,pBin))/(p_Array(nuc,pBin+1)-p_Array(nuc,pBin))
        
        if (val .GT. interpVal) then
c ------- momentum value not selected
          goto 11
        endif
        
        p(1) = 0
        p(2) = mag*sin(theta)*cos(phi)
        p(3) = mag*sin(theta)*sin(phi)
        p(4) = mag*cos(theta)

      return
      end
