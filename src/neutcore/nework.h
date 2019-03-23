C-------------------------------------------------------
C
C     INCLUDE FILE FOR NEUTRINO INTERACTION ( nework.h )
C
C     WORK AREA OF NEUTRINO INTERACTION
C  
C     CAUTION **** UNIT IS GEV ****
C
C-------------------------------------------------------
      INTEGER MAXNE, MODENE, NUMNE, IPNE, IORGNE, IFLGNE, ICRNNE
      REAL    PNE
      PARAMETER(MAXNE=100)
C
C     COMMON /NEWORK/
C  
C     MODENE      : MODE OF INTERACTION
C     NUMNE       : # OF PARTICLE
C     IPNE(I)     : PARTICLE CODE OF I-TH PARTICLE
C     PNE(3,I)    : MOMENTUM OF I-TH PARTICLE ( GEV/C )
C     IORGNE(I)   : ID OF ORIGIN PARTICLE
C     IFLGNE(I)   : FLAG OF FINAL STATE 
C                  -1 : initial particle
C                   0 : DETERMINED LATER PROCEDURE  
C                   1 : DECAY TO OTHER PARTICLE
C                   2 : ESCAPE FROM DETECTOR
C                   3 : ABSORPTION
C                   4 : CHARGE EXCHANGE
C                   5 : STOP AND NOT CONSIDER IN M.C. 
C                   6 : E.M. SHOWER
C     ICRNNE(I)  : FLAG OF TO CHASE OR NOT
C                   0 : DO NOT CHASE
C                   1 : CHASE
C
      COMMON /NEWORK/ MODENE,NUMNE,IPNE(MAXNE),
     $    PNE(3,MAXNE),IORGNE(MAXNE),IFLGNE(MAXNE),
     $    ICRNNE(MAXNE)

C-------------------------------------------------------------
