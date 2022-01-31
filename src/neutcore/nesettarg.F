**********************************************************
*     -------------------------
      subroutine nesettarg
*     -------------------------
*     
*     ( purpose )
*       Set parameters related to target nucleon
*
*     ( input ) Note:NUMATOM is from necard.h
*       NUMATOM : atomic number
*
*     ( output ) 
*       
*
*     ( creation date and author )
*      2007.11.05 ; G. Mitsuka
*     2008.11.07 ; R. Tacik    added additional nuclei
*     2016.08.05 ; P. Stowell  added free nucleons
*     
*     ( comment )
*      Charge density parameters refer
*      "Atomic Data and Nuclear Data Tables,
*       Vol.14, No.5,6, Nov/Dec 1974, Jager and Vries"
**********************************************************

      implicit none
#include "necard.h"
#include "efpion.h"

      if (NUMBNDP.eq.1 .and. NUMBNDN.eq.0) then !p
         RMSRAD= 0.01
         C     = 0.01
         CNN   = 0.01
         AF    = 0.01
         WPARM = 0.00
         DFACT = 1.0
      else if (NUMBNDP.eq.0 .and. NUMBNDN.eq.1) then !n
         RMSRAD= 0.01
         C     = 0.01
         CNN   = 0.01
         AF    = 0.01
         WPARM = 0.00
         DFACT = 1.0
      else if (NUMBNDP.eq.5) then ! 11B (estimate)
         RMSRAD= 2.42
         C     = 2.00
         CNN   = 2.00
         AF    = 0.51
         WPARM = 0.00
         DFACT = 1.0
      else if (NUMBNDP.eq.6) then ! 12C
         RMSRAD= 2.455
         C     = 2.355
         CNN   = 2.355
         AF    = 0.5224
         WPARM =-0.149
         DFACT = 1.0
      else if (NUMBNDP.eq.7) then ! 14N (estimate)
         RMSRAD=  2.524 
         C     =  2.570 
         CNN   =  2.570 
         AF    =  0.5052
         WPARM = -0.180 
         DFACT =  1.0         
      else if (NUMBNDP.eq.8) then ! 16O
         RMSRAD= 2.730
         C     = 2.69
         CNN   = 2.69
         AF    = 0.40961
         WPARM = 0.
         DFACT = 0.9985962
      else if (NUMBNDP.eq.9) then ! 19F
         RMSRAD= 2.900
         C     = 2.59
         CNN   = 2.59
         AF    = 0.564
         WPARM = 0.
         DFACT = 1.0
      else if (NUMBNDP.eq.11) then ! 23Na (estimate)
         RMSRAD= 2.94
         C     = 2.70
         CNN   = 2.70
         AF    = 0.56
         WPARM = 0.
         DFACT = 1.0
      else if (NUMBNDP.eq.13) then ! 27Al
         RMSRAD= 3.05
         C     = 2.84
         CNN   = 2.84
         AF    = 0.569
         WPARM = 0.
         DFACT = 1.0
      else if (NUMBNDP.eq.14) then ! 28Si (2pF)
         RMSRAD= 3.15
         C     = 3.14
         CNN   = 3.14
         AF    = 0.537
         WPARM = 0.
         DFACT = 1.0
      else if (NUMBNDP.eq.16) then ! 32S
         RMSRAD= 3.243
         C     = 3.503
         CNN   = 3.503
         AF    = 0.633
         WPARM = -0.250
         DFACT = 1.0   
      else if (NUMBNDP.eq.17) then ! 35Cl
         RMSRAD= 3.388
         C     = 3.476
         CNN   = 3.476
         AF    = 0.599
         WPARM =-0.10
         DFACT = 1.0
      else if (NUMBNDP.eq.18) then ! 40Ar
         RMSRAD= 3.48
         C     = 3.73
         CNN   = 3.73
         AF    = 0.62
         WPARM =-0.19
         DFACT = 1.0
      else if (NUMBNDP.eq.19) then ! 39K
         RMSRAD= 3.408
         C     = 3.743
         CNN   = 3.743
         AF    = 0.585
         WPARM =-0.201
         DFACT = 1.0
      else if (NUMBNDP.eq.20) then ! 40Ca
         RMSRAD= 3.482
         C     = 3.766
         CNN   = 3.766
         AF    = 0.586
         WPARM =-0.161
         DFACT = 1.0
      else if (NUMBNDP.eq.22) then ! 48Ti
         RMSRAD= 3.59
         C     = 3.75
         CNN   = 3.75
         AF    = 0.567
         WPARM = 0.0
         DFACT = 1.0
      else if (NUMBNDP.eq.25) then ! 55Mn (2pf)
         RMSRAD= 3.68
         C     = 3.89
         CNN   = 3.89
         AF    = 0.567
         WPARM = 0.
         DFACT = 1.0
      else if (NUMBNDP.eq.26) then ! 56Fe
         RMSRAD= 3.787
         C     = 3.971
         CNN   = 3.971
         AF    = 0.5935
         WPARM = 0.
         DFACT = 1.0
      else if (NUMBNDP.eq.27) then ! 59Co
         RMSRAD= 3.80
         C     = 4.08
         CNN   = 4.08
         AF    = 0.569
         WPARM = 0.
         DFACT = 1.0
      else if (NUMBNDP.eq.28) then ! 58Ni
         RMSRAD= 3.823
         C     = 4.14
         CNN   = 4.14
         AF    = 0.56
         WPARM = 0.
         DFACT = 1.0  
      else if (NUMBNDP.eq.29) then ! 63Cu
         RMSRAD= 3.925
         C     = 4.214
         CNN   = 4.214
         AF    = 0.586
         WPARM = 0.
         DFACT = 1.0
      else if (NUMBNDP.eq.30) then ! 64Zn
         RMSRAD= 3.965
         C     = 4.285
         CNN   = 4.285
         AF    = 0.584
         WPARM = 0.
         DFACT = 1.0
      else if (NUMBNDP.eq.38) then ! 88Sr
         RMSRAD= 4.17
         C     = 4.83
         CNN   = 4.83
         AF    = 0.496
         WPARM = 0.
         DFACT = 1.0
      else if (NUMBNDP.eq.40) then ! 90Zr (Estimate)
         RMSRAD= 4.274
         C     = 4.8
         CNN   = 4.8
         AF    = 0.57
         WPARM = 0.
         DFACT = 1.0   
      else if (NUMBNDP.eq.41) then ! 93Nb
         RMSRAD= 4.31
         C     = 4.87
         CNN   = 4.87
         AF    = 0.573
         WPARM = 0.
         DFACT = 1.0   
      else if (NUMBNDP.eq.45) then ! 103Rh (Estimate from 110Cd)
         RMSRAD= 4.6
         C     = 5.3
         CNN   = 5.3
         AF    = 0.55
         WPARM = 0.
         DFACT = 1.0   
      else if (NUMBNDP.eq.50) then ! 120Sn
         RMSRAD= 4.64
         C     = 5.315
         CNN   = 5.315
         AF    = 0.576
         WPARM = 0.
         DFACT = 1.0   
      else if (NUMBNDP.eq.67) then ! 165Ho
         RMSRAD= 5.23
         C     = 6.18
         CNN   = 6.18
         AF    = 0.57
         WPARM = 0.
         DFACT = 1.0 
      else if (NUMBNDP.eq.74) then ! 184W
         RMSRAD= 5.42
         C     = 6.51
         CNN   = 6.51
         AF    = 0.535
         WPARM = 0.
         DFACT = 1.0 
      else if (NUMBNDP.eq.79) then ! 197Au
         RMSRAD= 5.33
         C     = 6.38
         CNN   = 6.38
         AF    = 0.535
         WPARM = 0.
         DFACT = 1.0  
      else if (NUMBNDP.eq.82) then ! 208Pb
         RMSRAD= 5.521
         C     = 6.624
         CNN   = 6.624
         AF    = 0.549
         WPARM = 0.
         DFACT = 1.0
      else if (NUMBNDP.eq.83) then ! 209Bi
         RMSRAD= 5.51
         C     = 6.75
         CNN   = 6.75
         AF    = 0.468
         WPARM = 0.
         DFACT = 1.0
      else if (NUMATOM .EQ. 2 .OR. NUMATOM .EQ. 3 .OR.
     $        NUMATOM .EQ. 4 .OR. NUMATOM .EQ. 20) then
         RMSRAD= 2.455
         C     = 2.355
         CNN   = 2.355
         AF    = 0.5224
         WPARM =-0.149
         DFACT = 1.0
         write(*,*) "Parameters not set in nesettarg.F 
     $        - only valid for Effective SF"
      else
         print*, 'This target is not supported yet, see nesettarg.F'
         stop
      end if

      AN = NUMATOM
      ZZ = NUMBNDP
      cc2= 2.5*C

      A=1.80
      ANN=1.80
      DR=0.05

      return
      end
