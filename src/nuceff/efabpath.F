**********************************************************************
*     -----------------------------------------------------
      FUNCTION EFABPATH(PPI,RIN)
*     -----------------------------------------------------
*
*     ( purpose )
*       GIVE ABSORPTION MEAN FREE PATH **(-1)
*
*     ( input )
*       PPI      : momentum of pion(MeV/c)
*       R        : radius(fm)
*
*     ( output )
*       efabpath : absorption mean free path**(-1) (fm ** (-1))
*
*     ( Creation Date and Author )
*       2007.11.05 ; G.Mitsuka - add support other target than 16O
*                                using efconv2roxy function
C


      Implicit None
#include "necard.h"
#include "efpion.h"
      real*4   efabpath

      real*4   ppi,path1,path2,diffp,diffr,rin

      real*4   reff

      integer*4 i,j

      REAL*4    PATH(8,14),PPITBL(14),PATHA(8,7),PATHB(8,7)

      EQUIVALENCE (PATH(1,1),PATHA(1,1))
      EQUIVALENCE (PATH(1,8),PATHB(1,1))

      DATA PPITBL/  0.,130.,180.,200.,225.,250.,275.,300.,325.,
     $            350.,375.,400.,425.,480./

      DATA PATHA/
     $  0.119,0.109,0.077,2.83E-02,0.71E-02,1.14E-03,1.01E-04,4.90E-06,
     $  0.119,0.109,0.077,2.83E-02,0.71E-02,1.14E-03,1.01E-04,4.90E-06,
     $  0.237,0.217,0.154,5.66E-02,1.42E-02,2.27E-03,2.01E-04,9.79E-06,
     $  0.335,0.306,0.213,8.03E-02,1.79E-02,2.17E-03,1.30E-04,3.92E-06,
     $  0.479,0.441,0.309,0.116,2.22E-02,1.93E-03,7.24E-05,1.26E-06,
     $  0.587,0.549,0.392,0.146,2.38E-02,1.50E-03,3.63E-05,4.03E-07,
     $  0.601,0.569,0.413,0.148,2.01E-02,9.06E-04,1.44E-05,1.14E-07/

      DATA PATHB/
     $  0.537,0.511,0.372,0.124,1.34E-02,4.18E-04,4.52E-06,2.82E-08,
     $  0.444,0.425,0.309,9.37E-02,7.88E-03,1.71E-04,1.35E-06,7.32E-09,
     $  0.356,0.342,0.248,6.85E-02,4.50E-03,7.02E-05,4.42E-07,2.21E-09,
     $  0.282,0.272,0.197,5.00E-02,2.60E-03,3.07E-05,1.67E-07,7.90E-10,
     $  0.223,0.215,0.157,3.68E-02,1.54E-03,1.46E-05,7.20E-08,3.28E-10,
     $  0.175,0.170,0.124,2.73E-02,9.42E-04,7.50E-06,3.50E-08,1.55E-10,
     $  0.,0.,0.,0.,0.,0.,0.,0./

      real efconv2roxy,Rorg,Rorgcut,RMSRADOXY
      external efconv2roxy
      parameter (RMSRADOXY = 2.730)

      EFABPATH=0.

      IF (PPI.GE.480) RETURN

      Reff=Rin
C     R is vonverted to that in 16O
      if(NUMATOM.ne.16) Reff=efconv2roxy(Rin)

      Rorgcut = 7.*(RMSRAD/RMSRADOXY)
      IF (RIN.GT.Rorgcut) RETURN

      DO 10 I=1,14
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

      EFABPATH = (PATH1+(PATH2-PATH1)*DIFFR)
      RETURN
      END
