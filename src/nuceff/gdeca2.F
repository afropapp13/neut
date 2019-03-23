*
* $Id: gdeca2.F,v 1.1.1.1 2016/11/04 20:31:13 mdunkman Exp $
*
* $Log: gdeca2.F,v $
* Revision 1.1.1.1  2016/11/04 20:31:13  mdunkman
* Import of NEUT 5.3.3 v1r27 for use with T2KReWeight v1r27 (Winter OA 2017)
*
* Revision 1.1.1.1  1995/10/24 10:21:23  cernlib
* Geant
*
*
C#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GDECA2(XM0,XM1,XM2,PCM)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Simulates two body decay process with isotropic angular       *
C.    *  distribution in CMS.                                          *
C.    *                                                                *
C.    *    ==>Called by : GDECAY                                       *
C.    *       Author    G.Patrick  *********                           *
C.    *       Modified by P.de Perio for NEUT                          *
C.    ******************************************************************
C.
C#include "geant321/gconsp.inc"
      DIMENSION PCM(4,2)
C      DIMENSION RNDM(2)
C.
C.    ------------------------------------------------------------------
C.
C             Generate first decay product in CMS.
C
      E1=(XM0*XM0+XM1*XM1-XM2*XM2)/(2.*XM0)
      P1=SQRT(ABS((E1-XM1)*(E1+XM1)))
C
C             Isotropic decay angular distribution.
C
      COSTH=2.*RLU(IDUM)-1.
      IF(ABS(COSTH).GE.1.0) THEN
         COSTH=SIGN(1.,COSTH)
         SINTH=0.
      ELSE
         SINTH=SQRT((1.-COSTH)*(1.+COSTH))
      ENDIF
      PHI=6.238318*RLU(IDUM)
C
C             Polar co-ordinates to momentum components.
C
      PCM(1,1)=P1*SINTH*COS(PHI)
      PCM(2,1)=P1*SINTH*SIN(PHI)
      PCM(3,1)=P1*COSTH
      PCM(4,1)=E1
C
C             Generate second decay product.
C
      PCM(1,2)=-PCM(1,1)
      PCM(2,2)=-PCM(2,1)
      PCM(3,2)=-PCM(3,1)
      PCM(4,2)=SQRT(PCM(1,2)**2+PCM(2,2)**2+PCM(3,2)**2+XM2*XM2)
C
      END
