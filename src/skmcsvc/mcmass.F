C     
      SUBROUTINE MCMASS(ID,XMASS)
C     
C     ( purpose )
C       SET MASS ON GIVEN ID
C
C     ( input ) 
C       ID
C
C     ( output ) 
C       XMASS ( MEV/C**2 )
C
C     ( creation date and author )
C       1992.6.5 ; First version by K.Kaneyuki
C       1998.03.15 ;Add eta',lambda,etc by J.K.
C
C     ( comment )
C       PARTICLE       ID
C
C       GAMMA          22
C       Z              23      
C       W              24
C     
C       E-             11
C       MU-            13
C       TAU-           15
C       NEU-E          12
C       NEU-MU         14
C       NEU-TAU        16
C
C       PI+           211
C       PI0           111
C       ETA           221
C       ETA'          331
C       RO+(770)      213
C       RO0(770)      113
C       OMEGA(783)    223
C       a1(1260)    20113,20213
C       K+            321
C       K0            311
C       K0S           310
C       K0L           130
C       K*+(892)      323
C       K*0(892)      313
C  
C       P            2212
C       N            2112
C       LAMBDA       3122
C       
C     Deuteron        100045
C     Tritium         100046
C     Alpha           100047
C     Geantino        100048
C     He3             100049
C     O16             100069
C
C
      XMASS=0.0
      ID2=ABS(ID)
      IF (ID2.EQ.  11) XMASS=0.511
      IF (ID2.EQ.  12) XMASS=0.0
      IF (ID2.EQ.  13) XMASS=105.658
      IF (ID2.EQ.  14) XMASS=0.0
      IF (ID2.EQ.  15) XMASS=1776.99
      IF (ID2.EQ.  16) XMASS=0.0

      IF (ID2.EQ.  22) XMASS=0.0
      IF (ID2.EQ.  23) XMASS=91.173E+3
      IF (ID2.EQ.  24) XMASS=80.22E+3
      
      IF (ID2.EQ. 211) XMASS=139.568
      IF (ID2.EQ. 111) XMASS=134.973
      IF (ID2.EQ. 221) XMASS=548.8
      IF (ID2.EQ. 213) XMASS=770.0
      IF (ID2.EQ. 113) XMASS=770.0
      IF (ID2.EQ. 223) XMASS=782.0
      IF (ID2.EQ.20113) XMASS=1260.0
      IF (ID2.EQ.20213) XMASS=1260.0
      IF (ID2.EQ. 321) XMASS=493.646
      IF (ID2.EQ. 311) XMASS=497.671
      IF (ID2.EQ. 310) XMASS=497.671
      IF (ID2.EQ. 130) XMASS=497.671
      IF (ID2.EQ. 323) XMASS=892.1
      IF (ID2.EQ. 313) XMASS=892.1
      IF (ID2.EQ. 331) XMASS=957.77

      IF (ID2.EQ.2212) XMASS=938.272
      IF (ID2.EQ.2112) XMASS=939.566
      IF (ID2.EQ.3122) XMASS=1115.68

      IF (ID2.EQ.100045) XMASS=1875.61
      IF (ID2.EQ.100046) XMASS=2808.0
      IF (ID2.EQ.100047) XMASS=3727.42
      IF (ID2.EQ.100048) XMASS=0.0
      IF (ID2.EQ.100049) XMASS=2814.48
      IF (ID2.EQ.100069) XMASS=14899.1
     
C     CW: add some nuclear masses for the effective SF model
C     1u c^{2} = 931.49 MeV used for conversion
      IF (ID2.EQ.1000010020) XMASS=1876.12
      IF (ID2.EQ.1000020030) XMASS=2809.40
      IF (ID2.EQ.1000020040) XMASS=3728.38
      IF (ID2.EQ.1000060120) XMASS=11177.88
      IF (ID2.EQ.1000080160) XMASS=14899.10
      IF (ID2.EQ.1000100200) XMASS=18622.76
      IF (ID2.EQ.1000130270) XMASS=25133.03
      IF (ID2.EQ.1000180400) XMASS=37224.56
      IF (ID2.EQ.1000260560) XMASS=52102.83
      IF (ID2.EQ.1000290630) XMASS=58618.29
      IF (ID2.EQ.1000300640) XMASS=59549.36
      IF (ID2.EQ.1000822080) XMASS=193728.17

      RETURN
      END


