C     
C     TAUMC Common block for the properties of taus which have decayed 
C     in Neut.  Included energy, polarization decay mode
C     
C     First Added 2004/04/04 C.W. Walter
      
#include "taumc_com.h"
C-- now in taumc_com.h
C      COMMON / TAUMC  / TAUMOM, TAUDCY, TAUDIR, POLTAU
C      
C      REAL    TAUMOM            ! Momentum of Tau in GeV/c
C      REAL    TAUDIR(3)         ! Direction of Tau
C      INTEGER TAUDCY            ! Decay mode
C      REAL    POLTAU(3)         ! Polarization of Tau
C-- end here      
C
C     The info below is for using this common in a PAW ntuple bank.

      CHARACTER*60 TAUMCTAGS(1)
      CHARACTER*60 TAUMCTAG
      EQUIVALENCE (TAUMCTAG, TAUMCTAGS(1))
C     ***********************************************************************
      DATA TAUMCTAGS/
     &     'TAUMOM:R, TAUDCY:I, TAUDIR(3):R, POLTAU(3):R'/
