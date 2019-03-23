***********************************************************************
*     -------------------
      SUBROUTINE TAINIDCY
*     -------------------
*     
*     ( purpose )
*       set decay for tauola
*
*     ( input ) 
*       none
*
*     ( output ) 
*       COMMON  /TAUBRA/,/DECPAR/,/TAUKLE/,/JAKI/
*
*     ( creation date and author )
*       1995.10.24 ; First version by K.Kaneyuki
*
*     ( comment )
*
***********************************************************************
#include "taucom.h"
#include "tacard.h"
C
C    conventions of particles names
C                 K-,P-,K+,  K0,P-,KB,  K-,P0,K0
C                  3, 1,-3  , 4, 1,-4  , 3, 2, 4  ,
C                 P0,P0,K-,  K-,P-,P+,  P-,KB,P0
C                  2, 2, 3  , 3, 1,-1  , 1,-4, 2  ,
C                 ET,P-,P0   P-,P0,GM
C                  9, 1, 2  , 1, 2, 8
C
      DIMENSION NOPIK(6,NMODE),NPIK(NMODE)
CAM   outgoing multiplicity and flavors of multi-pion /multi-K modes    
      DATA   NPIK  /                4,                    4,  
     1                              5,                    5,
     2                              6,                    6,
     3                              3,                    3,            
     4                              3,                    3,            
     5                              3,                    3,            
     6                              3,                    3,  
     7                              2                         /         
      DATA  NOPIK / -1,-1, 1, 2, 0, 0,     2, 2, 2,-1, 0, 0,  
     1              -1,-1, 1, 2, 2, 0,    -1,-1,-1, 1, 1, 0,  
     2              -1,-1,-1, 1, 1, 2,    -1,-1, 1, 2, 2, 2, 
     3              -3,-1, 3, 0, 0, 0,    -4,-1, 4, 0, 0, 0,  
     4              -3, 2,-4, 0, 0, 0,     2, 2,-3, 0, 0, 0,  
     5              -3,-1, 1, 0, 0, 0,    -1, 4, 2, 0, 0, 0,  
     6               9,-1, 2, 0, 0, 0,    -1, 2, 8, 0, 0, 0,
     7              -3, 4, 0, 0, 0, 0                         /

C      REAL BRANCH(22)
C      DATA BRANCH/17.93, 17.58, 11.6, 24.0, 14.1,
C     &             0.67,  1.2 , 

C     
C     COMMON /TAUBRA/,/DECOMP/
C

      NCHAN = NMODE + 7
      DO 10 I = 1,22
         IF (I.LE.NCHAN) THEN
            JLIST(I) = I
            IF (ITAMODE.EQ.0) THEN

C     Taken from ATLAS software 
C     http://atlas.web.cern.ch/Atlas/GROUPS/SOFTWARE/OO/dist/
C     5.3.0/Generators/Tauola_i/Tauola_i-00-00-10/src/tauini.F
C
C     Branching rations acc to:
C     C. Caso et al, The European Physical Journal C3 (1998) 1 
C     and 1999 off-year partial update for the 2000 edition available on 
C     the PDG WWW pages (URL: http://pdg.lbl.gov/). 
               IF(I.EQ. 1) GAMPRT(I) = 0.17810 ! e nu nu~
               IF(I.EQ. 2) GAMPRT(I) = 0.17370 ! mu nu nu~
               IF(I.EQ. 3) GAMPRT(I) = 0.11080 ! pi- nu
               IF(I.EQ. 4) GAMPRT(I) = 0.25840 ! rho 
               IF(I.EQ. 5) GAMPRT(I) = 0.18250 ! A1
               IF(I.EQ. 6) GAMPRT(I) = 0.00710 ! K-
               IF(I.EQ. 7) GAMPRT(I) = 0.01290 ! K*-
               IF(I.EQ. 8) GAMPRT(I) = 0.04350 ! 2PI-  PI+  PI0
               IF(I.EQ. 9) GAMPRT(I) = 0.01110 ! PI-       3PI0
               IF(I.EQ.10) GAMPRT(I) = 0.00110 ! 2PI-  PI+ 2PI0
               IF(I.EQ.11) GAMPRT(I) = 0.00075 ! 3PI-, 2PI+
               IF(I.EQ.12) GAMPRT(I) = 0.00022 ! 3PI-, 2PI+,  PI0
               IF(I.EQ.13) GAMPRT(I) = 0.00140 ! 2PI-,  PI+, 3PI0
               IF(I.EQ.14) GAMPRT(I) = 0.00160 ! K-, PI-,  K+ 
               IF(I.EQ.15) GAMPRT(I) = 0.00121 ! K0, PI-, K0B
               IF(I.EQ.16) GAMPRT(I) = 0.00150 ! K-,  K0, PI0
               IF(I.EQ.17) GAMPRT(I) = 0.00080 ! 2PI0 K-
               IF(I.EQ.18) GAMPRT(I) = 0.00180 ! K-, PI-, PI+ 
               IF(I.EQ.19) GAMPRT(I) = 0.00390 ! PI-, K0B, PI0
               IF(I.EQ.20) GAMPRT(I) = 0.00170 ! ETA, PI-, PI0
               IF(I.EQ.21) GAMPRT(I) = 0.00150 ! PI-, PI0, GAM
               IF(I.EQ.22) GAMPRT(I) = 0.00159 ! K-,  K0
C
C              ORIGINAL VALUES FOR THIS FILE
C
C              Note: these aren't as different as those above
C              as you think.  The overall normalization 
C              isn't important.  These are normalized such that 
C              e and mu decays have branching ratios of 1 and are 
C              equal.
C
C               IF(I.EQ. 1) GAMPRT(I) = 1.0000
C               IF(I.EQ. 2) GAMPRT(I) = 1.0000
C               IF(I.EQ. 3) GAMPRT(I) = 0.605
C               IF(I.EQ. 4) GAMPRT(I) = 1.326
C               IF(I.EQ. 5) GAMPRT(I) = 0.719
C               IF(I.EQ. 6) GAMPRT(I) = 0.040
C               IF(I.EQ. 7) GAMPRT(I) = 0.070
C               IF(I.EQ. 8) GAMPRT(I) = 0.0486+0.03308
C               IF(I.EQ. 9) GAMPRT(I) = 0.0157
C               IF(I.EQ.10) GAMPRT(I) = 0.07
C               IF(I.EQ.11) GAMPRT(I) = 0.029
C               IF(I.EQ.12) GAMPRT(I) = 0.004
C               IF(I.EQ.13) GAMPRT(I) = 0.004
C               IF(I.EQ.14) GAMPRT(I) = 0.00606
C               IF(I.EQ.15) GAMPRT(I) = 0.00606
C               IF(I.EQ.16) GAMPRT(I) = 0.00053
C               IF(I.EQ.17) GAMPRT(I) = 0.00599
C               IF(I.EQ.18) GAMPRT(I) = 0.03147
C               IF(I.EQ.19) GAMPRT(I) = 0.03237
C               IF(I.EQ.20) GAMPRT(I) = 0.01082
C               IF(I.EQ.21) GAMPRT(I) = 0.003
C               IF(I.EQ.22) GAMPRT(I) = 0.018
            ELSE IF (ITAMODE.EQ.-1) THEN
               GAMPRT(I) = TAUBRA(I)
            ELSE IF (ITAMODE.GE.1 .AND. ITAMODE.LE.22) THEN
               IF (I.EQ.ITAMODE) GAMPRT(I) = 1.0
            ELSE
               WRITE(6,*) ' TAINIDCY : ITAMODE IS BAD !!! '
               WRITE(6,*) ' TAINIDCY : CHECK YOUR CARD INPUT '
               STOP
            ENDIF
            IF(I.EQ. 8) NAMES(I-7)='  TAU-  --> 2PI-,  PI0,  PI+   '
            IF(I.EQ. 9) NAMES(I-7)='  TAU-  --> 3PI0,        PI-   '
            IF(I.EQ.10) NAMES(I-7)='  TAU-  --> 2PI-,  PI+, 2PI0   '
            IF(I.EQ.11) NAMES(I-7)='  TAU-  --> 3PI-, 2PI+,        '
            IF(I.EQ.12) NAMES(I-7)='  TAU-  --> 3PI-, 2PI+,  PI0   '
            IF(I.EQ.13) NAMES(I-7)='  TAU-  --> 2PI-,  PI+, 3PI0   '
            IF(I.EQ.14) NAMES(I-7)='  TAU-  -->  K-, PI-,  K+      '
            IF(I.EQ.15) NAMES(I-7)='  TAU-  -->  K0, PI-, K0B      '
            IF(I.EQ.16) NAMES(I-7)='  TAU-  -->  K-,  K0, PI0      '
            IF(I.EQ.17) NAMES(I-7)='  TAU-  --> PI0, PI0,  K-      '
            IF(I.EQ.18) NAMES(I-7)='  TAU-  -->  K-, PI-, PI+      '
            IF(I.EQ.19) NAMES(I-7)='  TAU-  --> PI-, K0B, PI0      '
            IF(I.EQ.20) NAMES(I-7)='  TAU-  --> ETA, PI-, PI0      '
            IF(I.EQ.21) NAMES(I-7)='  TAU-  --> PI-, PI0, GAM      '
            IF(I.EQ.22) NAMES(I-7)='  TAU-  -->  K-,  K0           '
         ELSE
            JLIST(I) = 0
            GAMPRT(I) = 0.
         ENDIF
 10   CONTINUE

      DO 20 I=1,NMODE
        MULPIK(I)=NPIK(I)
        DO 30 J=1,MULPIK(I)
           IDFFIN(J,I)=NOPIK(J,I)
 30     CONTINUE
 20   CONTINUE

C
C     COMMON /TAUKLE/
C      
      BRA1=0.5
      BRK0=0.5
      BRK0B=0.5
      BRKS=0.6667

C
C     COMMON /JAKI/
C      
      JAK1=0
      JAK2=0
      
      RETURN
      END
      
