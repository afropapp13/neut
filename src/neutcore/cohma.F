************************************************************************
*     ----------------------------
      FUNCTION COHMA(E,X,Y,XMLEP)
*     ----------------------------
*
*     (Purpose)
*       RETURN DIFFERENTIAL CROSSSECTION WHICH DEPENDS X,Y,Z,PHI.
*
*
*     (Input)
*       E     : INCIDENT NEUTRINO ENERGY ( GEV )
*       COMMON /nework/
*
*     (Output)
*       DCOHNE : DIFFERENTIAL CROSSECTION
*     
*     (Creation Date and Author)
*       1994.08.30 ; G.DION & K.Okumura
*       1996.02.05 ; K.Kaneyuki   DCOH -> DCOHNE
*       2007.01.08 ; G.Mitsuka,  Add lepton mass correction,
*                                rename from DCOHNE to COHMA
*       2007.11.05 ; G.Mitsuka,  support target nucleus besides 16O
*                                Atomic number is given by "NUMATOM"
*       2010.10.29 ; R.Terri, Simplified for reweighting program
*
*     (Comment)
*       Based on the function dcohnerein.F (make sure things are synced)
*
************************************************************************

#include "necard.h"
#include "nework.h"
#include "neutmodel.h"
#include "pdf804/w50516.inc"

C     COMMON /NVECT/ IMOD,LVECT,IP(20),PIN(3,20)
C
      REAL*8 SOLN,corr
      real QSMX,QSMN
      real XMAX,XMIN
      real PRMASS,PIMASS
      real Y2MPE2

C     NEED Proton & pion Mass here
      PRMASS = .9383
C     PIMASS TAKEN FROM RSCONS.h
      PIMASS = .14

      COHMA = 0.
      corr = 1.
      SOLN = 0.

      if(Y.lt.PIMASS/E) goto 100
      if(Y.gt.1./(1.+PRMASS*X/(2.*E))) goto 100

      QSMX=4.*E*E*(1.-Y)
      QSMN=0.
      if(QSMX.lt.0.) goto 100

      XMAX=QSMX/(2.*PRMASS*Y*E)
      XMIN=QSMN/(2.*PRMASS*Y*E)
      if(X.gt.XMAX) goto 100
      if(X.lt.XMIN) goto 100


C  ((OUTGOING PION MOMENTUM) / (INCIDENT NEUTRINO ENERGY))**2
C      Y2MPE2=Y**2-(PIMASS/E)**2

      SOLN = (1.-Y)/(1.+2.*PRMASS*E*X*Y/XMACOH**2)**2

      COHMA = SOLN

 100  continue

      RETURN
      END
