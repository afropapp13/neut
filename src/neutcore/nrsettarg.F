**********************************************************
*     -------------------------
      subroutine nrsettarg
*     -------------------------
*     
*     ( purpose )
*       Set parameters related to target nucleon
*
*     ( input )
*       
*
*     ( output ) 
*       
*
*     ( creation date and author )
*      2007.11.05 ; G. Mitsuka
*      2008.11.07 ; R. Tacik    added additional nuclei
*     
*     ( comment )
*      Charge density parameters refer
*      "Atomic Data and Nuclear Data Tables,
*       Vol.14, No.5,6, Nov/Dec 1974, Jager and Vries"
*
*      Parameters for 16O is not same as ones in nesettarg.F
*
*     2010.11.01 ; P. de Perio - Change 16O oxygen to be same
*                                as ones in nesettarg.F
**********************************************************

      implicit none
#include "necard.h"
#include "neutparams.h"
#include "nrnuclparam.h"

      APN = float(NUMATOM)
      AP  = float(NUMBNDP)

      if (NUMBNDP.eq.5) then ! 11B (estimate)
         NRRMSRAD= 2.42
         NRC     = 2.00
         NRCNN   = 2.00
         NRAF    = 0.51
         NRWPARM = 0.00
         NRDFACT = 1.0
      else if (NUMBNDP.eq.6) then ! 12C
         NRRMSRAD= 2.455
         NRC     = 2.355
         NRCNN   = 2.355
         NRAF    = 0.5224
         NRWPARM =-0.149
         NRDFACT = 1.0
      else if (NUMBNDP.eq.7) then ! 14N
         NRRMSRAD= 2.524
         NRC     = 2.570
         NRCNN   = 2.570
         NRAF    = 0.5052
         NRWPARM =-0.180
         NRDFACT = 1.0
      else if (NUMBNDP.eq.8) then ! 16O
         NRRMSRAD= 2.730
         NRC     = 2.69
         NRCNN   = 2.69
         NRAF    = 0.40961
         NRWPARM = 0.
         NRDFACT = 0.9985962
      else if (NUMBNDP.eq.9) then ! 19F
         NRRMSRAD= 2.900
         NRC     = 2.59
         NRCNN   = 2.59
         NRAF    = 0.564
         NRWPARM = 0.
         NRDFACT = 1.0
      else if (NUMBNDP.eq.11) then ! 23Na (estimate)
         NRRMSRAD= 2.94
         NRC     = 2.70
         NRCNN   = 2.70
         NRAF    = 0.56
         NRWPARM = 0.
         NRDFACT = 1.0
      else if (NUMBNDP.eq.13) then ! 27Al
         NRRMSRAD= 3.05
         NRC     = 2.84
         NRCNN   = 2.84
         NRAF    = 0.569
         NRWPARM = 0.
         NRDFACT = 1.0
      else if (NUMBNDP.eq.14) then ! 28Si (2pF)
         NRRMSRAD= 3.15
         NRC     = 3.14
         NRCNN   = 3.14
         NRAF    = 0.537
         NRWPARM = 0.
         NRDFACT = 1.0
      else if (NUMBNDP.eq.16) then ! 32S
         NRRMSRAD= 3.243 
         NRC     = 3.503 
         NRCNN   = 3.503 
         NRAF    = 0.633 
         NRWPARM = -0.250
         NRDFACT = 1.0   
      else if (NUMBNDP.eq.17) then ! 35Cl
         NRRMSRAD= 3.388
         NRC     = 3.476
         NRCNN   = 3.476
         NRAF    = 0.599
         NRWPARM =-0.10
         NRDFACT = 1.0
      else if (NUMBNDP.eq.18) then ! 40Ar
         NRRMSRAD= 3.48
         NRC     = 3.73
         NRCNN   = 3.73
         NRAF    = 0.62
         NRWPARM =-0.19
         NRDFACT = 1.0
      else if (NUMBNDP.eq.20) then ! 40Ca
         NRRMSRAD= 3.482
         NRC     = 3.766
         NRCNN   = 3.766
         NRAF    = 0.586
         NRWPARM =-0.161
         NRDFACT = 1.0
      else if (NUMBNDP.eq.22) then ! 48Ti
         NRRMSRAD= 3.59
         NRC     = 3.75
         NRCNN   = 3.75
         NRAF    = 0.567
         NRWPARM = 0.0
         NRDFACT = 1.0
      else if (NUMBNDP.eq.26) then ! 56Fe
         NRRMSRAD= 3.787
         NRC     = 3.971
         NRCNN   = 3.971
         NRAF    = 0.5935
         NRWPARM = 0.
         NRDFACT = 1.0
      else if (NUMBNDP.eq.27) then ! 59Co
         NRRMSRAD= 3.80
         NRC     = 4.08
         NRCNN   = 4.08
         NRAF    = 0.569
         NRWPARM = 0.
         NRDFACT = 1.0
      else if (NUMBNDP.eq.28) then ! 58Ni
         NRRMSRAD= 3.823
         NRC     = 4.14 
         NRCNN   = 4.14 
         NRAF    = 0.56 
         NRWPARM = 0.   
         NRDFACT = 1.0  
      else if (NUMBNDP.eq.29) then ! 63Cu
         NRRMSRAD= 3.925
         NRC     = 4.214
         NRCNN   = 4.214
         NRAF    = 0.586
         NRWPARM = 0.
         NRDFACT = 1.0
      else if (NUMBNDP.eq.30) then ! 64Zn
         NRRMSRAD= 3.965
         NRC     = 4.285
         NRCNN   = 4.285
         NRAF    = 0.584
         NRWPARM = 0.
         NRDFACT = 1.0
      else if (NUMBNDP.eq.38) then ! 88Sr
         NRRMSRAD= 4.17
         NRC     = 4.83
         NRCNN   = 4.83
         NRAF    = 0.496
         NRWPARM = 0.
         NRDFACT = 1.0
      else if (NUMBNDP.eq.40) then ! 90Zr (Estimate)
         NRRMSRAD= 4.274 
         NRC     = 4.8   
         NRCNN   = 4.8   
         NRAF    = 0.57  
         NRWPARM = 0.    
         NRDFACT = 1.0   
      else if (NUMBNDP.eq.41) then ! 93Nb
         NRRMSRAD= 4.31  
         NRC     = 4.87  
         NRCNN   = 4.87  
         NRAF    = 0.573 
         NRWPARM = 0.    
         NRDFACT = 1.0 
      else if (NUMBNDP.eq.45) then ! 103Rh (Estimate from 110Cd)
         NRRMSRAD= 4.6   
         NRC     = 5.3   
         NRCNN   = 5.3   
         NRAF    = 0.55  
         NRWPARM = 0.    
         NRDFACT = 1.0  
      else if (NUMBNDP.eq.50) then ! 120Sn
         NRRMSRAD= 4.64  
         NRC     = 5.315 
         NRCNN   = 5.315 
         NRAF    = 0.576 
         NRWPARM = 0.    
         NRDFACT = 1.0   
      else if (NUMBNDP.eq.67) then ! 165Ho
         NRRMSRAD= 5.23 
         NRC     = 6.18 
         NRCNN   = 6.18 
         NRAF    = 0.57 
         NRWPARM = 0.   
         NRDFACT = 1.0 
      else if (NUMBNDP.eq.79) then ! 197Au
         NRRMSRAD= 5.33  
         NRC     = 6.38  
         NRCNN   = 6.38  
         NRAF    = 0.535 
         NRWPARM = 0.    
         NRDFACT = 1.0  
      else if (NUMBNDP.eq.82) then ! 208Pb
         NRRMSRAD= 5.521
         NRC     = 6.624
         NRCNN   = 6.624
         NRAF    = 0.549
         NRWPARM = 0.
         NRDFACT = 1.
      else if (NUMBNDP.eq.83) then ! 209Bi
         NRRMSRAD= 5.51
         NRC     = 6.75
         NRCNN   = 6.75
         NRAF    = 0.468
         NRWPARM = 0.
         NRDFACT = 1.0
      else
         print*, 'This target is not supported yet, see nrsettarg.F'
         stop
      end if

      NRPFSURF = PFSURF

      call nrinit(apn,1)

      return
      end
