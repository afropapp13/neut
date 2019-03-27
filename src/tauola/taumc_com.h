C     
C     TAUMC Common block for the properties of taus which have decayed 
C     in Neut.  Included energy, polarization decay mode
C     
C     First Added 2004/04/04 C.W. Walter
      
#ifndef TAUMC_COM_INCLUDED
  
      COMMON / TAUMC  / TAUMOM, TAUDCY, TAUDIR, POLTAU
      
      REAL    TAUMOM            ! Momentum of Tau in GeV/c
      REAL    TAUDIR(3)         ! Direction of Tau
      INTEGER TAUDCY            ! Decay mode
      REAL    POLTAU(3)         ! Polarization of Tau
      
#define TAUMC_COM_INCLUDED
#endif
