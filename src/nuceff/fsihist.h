C-------------------------------------------------------------
C
C     INCLUDE FILE FOR NUCEFF VERTICES ( fsihist.h )
C
C     WORK AREA OF M.C. VECTOR GENERATOR to store FSI vertices
C
C-------------------------------------------------------------
      INTEGER MAXVERT
      PARAMETER(MAXVERT=100)

      INTEGER MAXVCVERT
      PARAMETER(MAXVCVERT=300)
C
C     NVERT       : # OF VERTICES
C     POSVERT(3,I): POSITION OF I-TH VERTEX
C     IFLGVERT(I) : INTERACTION TYPE OF I-TH VERTEX
C                     (*10 FOR HI-NRG)
C                     (*100 for SKDETSIM non-NEFFECT interactions e.g. elastic SGPIEL;
C                          +0 Free Hydrogen, +1 Oxygen)
C                  -1 : ESCAPE
C                   0 : INITIAL (or unmatched parent vertex if I>1)
C                   1 :
C                   2 :
C                   3 : ABSORPTION
C                   4 : CHARGE EXCHANGE
C                   5 :
C                   6 :
C                   7 : HADRON PRODUCTION (hi-nrg only, i.e. 70)
C                   8 : QUASI-ELASTIC SCATTER
C                   9 : FORWARD (ELASTIC-LIKE) SCATTER
C
C     NVCVERT     : # OF INTERMEDIATE PARTICLES
C     DIRVERT(3,J): DIRECTION OF J-TH PARTICLE
C     ABSPVERT(J) : ABSOLUTE MOM. OF J-TH PART. in lab frame (MeV/c)
C     ABSTPVERT(J): ABSOLUTE MOM. OF J-TH PART. in nucleon rest frame
C     IPVERT(J)   : PARTICLE CODE OF J-TH PARTICLE
C     IVERTI(J)   : INDEX OF INITIAL VERTEX OF J-TH PARTICLE
C     IVERTF(J)   : INDEX OF FINAL VERTEX OF J-TH PARTICLE
C
      INTEGER NVERT, IFLGVERT, NVCVERT, IPVERT, IVERTI, IVERTF
      REAL    POSVERT, DIRVERT, ABSPVERT, ABSTPVERT, FSIPROB
      COMMON /FSIHIST/ NVERT,POSVERT(3,MAXVERT),
     $    IFLGVERT(MAXVERT), NVCVERT, DIRVERT(3,MAXVCVERT),
     $    ABSPVERT(MAXVCVERT), ABSTPVERT(MAXVCVERT),
     $    IPVERT(MAXVCVERT),IVERTI(MAXVCVERT),IVERTF(MAXVCVERT),
     $    FSIPROB

