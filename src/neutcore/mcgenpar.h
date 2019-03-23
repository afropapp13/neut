C-------------------------------------------------------
C
C     INCLUDE FILE FOR MC ( mcgenpar.h )
C
C     FLAG OF NUCLEAR INTERATION
C
C-------------------------------------------------------
C
C     COMMON /NEGENPAR/
C
C     NPARENT : PARENT NUCLEON ( NOT YET SUPPORTED )
C               0 : H2O ( DEFAULT )
C               1 : H2
C               2 : O
C
C     NUCFLG1 : FLAG OF NUCLEAR INTERATION IN O16
C               0 : CONSIDER NUCLEAR INT. ( DEFAULT )
C               1 : IGNORE NUCLEAR INT.
C
C     NUCFLG2 : FLAG OF NUCLEAR INTERATION IN WATER
C               0 : CONSIDER NUCLEAR INT. ( DEFAULT )
C               1 : IGNORE NUCLEAR INT.
C
      INTEGER NPAERNT, NUCFLG1, NUCFLG2
      COMMON /NEGENPAR/ NPAERNT, NUCFLG1, NUCFLG2


