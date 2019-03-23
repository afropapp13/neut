*     -----------------------------------------------------
      FUNCTION EFQEPATHTUN(PPI,RIN)      
*     -----------------------------------------------------
*     
*     ( purpose )
*       GIVE INELASTIC SCATTER MEAN FREE PATH **(-1)
*
*     ( input ) 
*       PPI   : PION MOMENUTM(MeV/c)
*       r     : radius
*
*     ( output ) 
*       EFQEPATHTUN : REACTION MEAN FREE PATH **(-1)  UNIT (/FERMI)
*
*     ( Creation Date and Author )
*       2007.11.05 ; G.Mitsuka - add support other target than 16O
*                                using efconv2roxy function
C     
C       2010.06.28 ; P. de Perio - Tuned to pion scattering data
C                                - Extended to 10 GeV/c

      IMPLICIT NONE

#include "necard.h"
#include "efpion.h"

      REAL*4 EFQEPATHTUN

      REAL*4 PPI
      REAL*4 DIFFP,DIFFR
      REAL*4 PATH1,PATH2, PATH3,RATIO

      integer*4 i,j

      REAL*4 REFF
      
      REAL*4    PATH(8,19),PPITBL(19)
      REAL*4    OPATH(18)
      DATA PPITBL/  0.,110.,180.,200.,225.,250.,280.,300.,325.,350.,
     $            375.,400.,425.,500.,500.1,850.,1000.,1400.,10000./
      real efconv2roxy,Rorg,Rorgcut,RMSRADOXY,RIN
      external efconv2roxy
      parameter (RMSRADOXY = 2.730)

      DATA PATH/0,0,0,0,0,0,0,0,

     $  3.66E-02,3.49E-02,2.87E-02,1.60E-02,4.77E-03,1.03E-03,
     $     1.52E-04,1.50E-05,
     $  2.71E-01,2.58E-01,2.12E-01,1.19E-01,3.52E-02,7.64E-03,
     $     1.12E-03,1.11E-04,
     $  4.36E-01,4.17E-01,3.42E-01,1.91E-01,5.72E-02,1.09E-02,
     $     1.29E-03,1.04E-04,
     $  6.41E-01,6.24E-01,5.29E-01,2.99E-01,9.01E-02,1.51E-02,
     $     1.47E-03,1.01E-04,
     $  6.76E-01,6.66E-01,5.86E-01,3.45E-01,1.02E-01,1.51E-02,
     $     1.37E-03,8.88E-05,
     $  6.92E-01,6.88E-01,6.26E-01,3.80E-01,1.10E-01,1.51E-02,
     $     1.22E-03,7.40E-05,

      
     $  7.16E-01,7.13E-01,6.53E-01,3.88E-01,1.03E-01,1.22E-02,
     $     8.78E-04,5.15E-05,
     $  7.43E-01,7.39E-01,6.75E-01,3.85E-01,8.92E-02,9.08E-03,
     $     6.00E-04,3.39E-05,
     $  7.74E-01,7.69E-01,6.97E-01,3.77E-01,7.68E-02,6.85E-03,
     $     4.28E-04,2.35E-05,
     $  7.97E-01,7.91E-01,7.11E-01,3.66E-01,6.62E-02,5.34E-03,
     $     3.21E-04,1.74E-05,
     $  8.14E-01,8.07E-01,7.21E-01,3.55E-01,5.79E-02,4.33E-03,
     $     2.54E-04,1.35E-05,
     $  8.31E-01,8.23E-01,7.30E-01,3.46E-01,5.17E-02,3.64E-03,
     $     2.10E-04,1.11E-05,
     $  8.66E-01,8.58E-01,7.62E-01,3.60E-01,5.39E-02,3.79E-03,
     $     2.19E-04,1.16E-05,


     $  7.94E-01,7.87E-01,6.98E-01,3.31E-01,4.94E-02,3.48E-03,2.01E-04,
     $     1.06E-05,
     $  1.27,1.26,1.12,5.29E-01,7.91E-02,5.57E-03,3.22E-04,1.70E-05,
     $  1.84,1.83,1.62,7.67E-01,1.15E-01,8.08E-03,4.67E-04,2.46E-05,
     $  2.58,2.56,2.27,1.07,1.60E-01,1.13E-02,6.53E-04,3.45E-05,
     $  0,0,0,0,0,0,0,0   /

      DATA OPATH/.881,.752,.666,.612,.578,.570,.571,.580,.590,
     $           .600,.600,.600,.610,.610,.619,.638,.667,.696/
      
      EFQEPATHTUN = 0.

      REFF=RIN
      if(NUMATOM.ne.16) REFF=efconv2roxy(RIN)

      IF (PPI.GT.425 .AND. PPI.LE.500) goto 100
      IF (PPI.GT.10000) RETURN

      Rorgcut = 7.*(RMSRAD/RMSRADOXY)
      IF (RIN.GT.Rorgcut) RETURN

      DO 10 I=1,19
         IF (PPI.LT.PPITBL(I)) GOTO 15
 10   CONTINUE
 15   CONTINUE

      DO 20 J=1,15
         IF (REFF.LT.(J-1)) GOTO 25
 20   CONTINUE
 25   CONTINUE

      DIFFP = (PPI-PPITBL(I-1))/(PPITBL(I)-PPITBL(I-1))
      DIFFR = (REFF-(J-2))
      PATH1 = PATH(J-1,I-1)+(PATH(J-1,I)-PATH(J-1,I-1))*DIFFP
      PATH2 = PATH(J,I-1)+(PATH(J,I)-PATH(J,I-1))*DIFFP

      EFQEPATHTUN = (PATH1+(PATH2-PATH1)*DIFFR)

      RETURN
         

C     Weird extrapolation 
 100  I=INT((PPI-425.)/25.)+1

      DO 120 J=1,15
         IF (REFF.LT.(J-1)) GOTO 125
 120  CONTINUE
 125  CONTINUE

      DIFFR = (REFF-(J-2))
      PATH1 = PATH(J-1,13)
      PATH2 = PATH(J,13)
      PATH3 = PATH(1,13)

      RATIO = (PATH1+(PATH2-PATH1)*DIFFR)/PATH3 * 0.66

      DIFFP = (PPI-(I-1)*25.-425.)/25.
      PATH1 = OPATH(I)
      PATH2 = OPATH(I+1)

      EFQEPATHTUN = RATIO*(PATH1 + (PATH2-PATH1)*DIFFP)
      EFQEPATHTUN = EFQEPATHTUN * (PPI*2.500E-3 + 0.35000)

      RETURN
      END

