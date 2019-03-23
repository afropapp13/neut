C     This routine has the common blocks for the hepevt record as
C     used in Tauola.

C     This is the cersion used in Tauola 2.7
C
C

C     this is the hepevt class in old style. No d_h_ class pre-name
      INTEGER NMXHEP
      PARAMETER (NMXHEP=4000)
      REAL*8  phep,  vhep ! to be real*4/ *8  depending on host
      INTEGER nevhep,nhep,isthep,idhep,jmohep,
     $        jdahep
      COMMON /hepevt/
     $      nevhep,               ! serial number
     $      nhep,                 ! number of particles
     $      isthep(nmxhep),   ! status code
     $      idhep(nmxhep),    ! particle ident KF
     $      jmohep(2,nmxhep), ! parent particles
     $      jdahep(2,nmxhep), ! childreen particles
     $      phep(5,nmxhep),   ! four-momentum, mass [GeV]
     $      vhep(4,nmxhep)    ! vertex [mm]
* ----------------------------------------------------------------------
      LOGICAL qedrad
      COMMON /phoqed/ 
     $     qedrad(nmxhep)    ! Photos flag
* ----------------------------------------------------------------------
      SAVE hepevt,phoqed


C     This was the version used in tauola version 2.6 and previous versions
C
C      INTEGER NMXHEP, NEVHEP, NHEP, ISTHEP, IDHEP, JMOHEP, JDAHEP
C      REAL    PHEP, VHEP
C      PARAMETER (NMXHEP=2000)
C      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
C     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
