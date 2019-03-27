**********************************************************************
*     -----------------------------------------------------
      FUNCTION EFABPATHTUN(PPI,RIN)
*     -----------------------------------------------------
*
*     ( purpose )
*       GIVE ABSORPTION MEAN FREE PATH **(-1)
*
*     ( input )
*       PPI      : Momentum of pion(MeV/c)
*       R        : radius(fm)
*
*     ( output )
*       efabpathtun : absorption mean free path**(-1) (fm ** (-1))
*
*     ( Creation Date and Author )
*       2007.11.05 ; G.Mitsuka - add support other target than 16O
*                                using efconv2roxy function
C
C       2010.06.28 ; P. de Perio - Tuned to pion scattering data
C                                - Extended to 2 GeV/c


      Implicit None
#include "necard.h"
#include "efpion.h"
      real*4   efabpathtun

      real*4   ppi,path1,path2,diffp,diffr,rin

      real*4   reff

      integer*4 i,j

      REAL*4    PATH(8,15),PPITBL(15),PATHA(8,9),PATHB(8,6)

      EQUIVALENCE (PATH(1,1),PATHA(1,1))
      EQUIVALENCE (PATH(1,10),PATHB(1,1))

      DATA PPITBL/  0.,50.,100.,130.,180.,200.,225.,250.,275.,300.,325.,
     $            350.,375.,400.,500./

      DATA PATHA/
     $  8.49E-02,7.78E-02,5.49E-02,2.02E-02,5.06E-03,8.13E-04,7.20E-05,
     $     3.50E-06,
     $  4.43E-02,4.05E-02,2.86E-02,1.05E-02,2.64E-03,4.24E-04,3.76E-05,
     $     1.82E-06,
     $  6.53E-02,5.98E-02,4.23E-02,1.55E-02,3.90E-03,6.26E-04,5.54E-05,
     $     2.69E-06,
     $  1.04E-01,9.52E-02,6.72E-02,2.47E-02,6.20E-03,9.95E-04,
     $     8.82E-05,4.28E-06,
     $  2.34E-01,2.14E-01,1.52E-01,5.59E-02,1.40E-02,2.24E-03,
     $     1.99E-04,9.68E-06,
     $  3.47E-01,3.17E-01,2.20E-01,8.31E-02,1.85E-02,2.25E-03,
     $     1.35E-04,4.06E-06,
     $  5.23E-01,4.82E-01,3.38E-01,1.27E-01,2.42E-02,2.11E-03,
     $     7.91E-05,1.38E-06,
     $  6.08E-01,5.68E-01,4.06E-01,1.51E-01,2.46E-02,1.55E-03,
     $     3.76E-05,4.17E-07,
     $  6.36E-01,6.02E-01,4.37E-01,1.57E-01,2.13E-02,9.59E-04,
     $     1.52E-05,1.21E-07/

      DATA PATHB/
     $  4.59E-01,4.37E-01,3.18E-01,1.06E-01,1.15E-02,3.58E-04,
     $     3.87E-06,2.41E-08,
     $  2.42E-01,2.32E-01,1.69E-01,5.11E-02,4.30E-03,9.33E-05,
     $     7.36E-07,3.99E-09,
     $  8.54E-02,8.21E-02,5.95E-02,1.64E-02,1.08E-03,1.68E-05,
     $     1.06E-07,5.30E-10,
     $  6.41E-02,6.16E-02,4.46E-02,1.23E-02,8.10E-04,1.26E-05,
     $     7.96E-08,3.98E-10,
     $  4.27E-02,4.10E-02,2.98E-02,8.22E-03,5.40E-04,8.42E-06,
     $     5.30E-08,2.65E-10,
     $  0,0,0,0,0,0,0,0/

      real efconv2roxy,Rorg,Rorgcut,RMSRADOXY
      external efconv2roxy
      parameter (RMSRADOXY = 2.730)

      EFABPATHTUN=0.

      IF (PPI.GE.500) RETURN

      Reff=Rin
C     R is vonverted to that in 16O
      if(NUMATOM.ne.16) Reff=efconv2roxy(Rin)

      Rorgcut = 7.*(RMSRAD/RMSRADOXY)
      IF (RIN.GT.Rorgcut) RETURN

      DO 10 I=1,15
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

      EFABPATHTUN = (PATH1+(PATH2-PATH1)*DIFFR)




      RETURN
      END
