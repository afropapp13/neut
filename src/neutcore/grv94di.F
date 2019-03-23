*
* $Id: grv94di.F,v 1.1 2007-01-20 07:09:16 skrep Exp $
*
* $Log: not supported by cvs2svn $
* Revision 1.3  2006/08/10 08:03:12  skrep
* Add Bodek Correction for GRV94di and GRV98_lo except for charm production
* (Aug'06 G.Mitsuka)
*
* Revision 1.2  2003/10/21 10:05:09  skrep
* Merging between main truck and ap02b branch applied
*
* Revision 1.1.2.1  2003/04/09 15:28:44  saji
* fnapicrs.F:
* 	bug of p/n ratio was fixed.
* grv94di.F nejet.F pdfset.F qgdisg.F:
* 	GRV94 in PDFLIB(not PYTHIA) package is
* 	used for caliculation of neutrino kinematics
* 	(in old version, CTEQ was used).
*
* 	Bodek modification is correctly applied
* 	in low Q2(Q2<Q2MIN) region.
*
* Revision 1.1  2002/08/29 13:46:22  saji
* Add modifications and corrections to GRV94 PDF according to
* Bodek's paper[hep-ex/0203009].
* If no need this effect or need original one,
* you should do 'imake_boot -DNO_GRV94CORR'
*
* Revision 1.1.1.2  1996/10/30 08:28:34  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:22  plothow
* Version 7.01
*
*
#include "pdf804/pilot.h"
*
*...NLO PARAMETRIZATION (DIS) :
*
       SUBROUTINE GRV94DI (ZX,ZQ,ZUV,ZDV,ZUB,ZDB,ZSB,ZCB,ZBB,ZTB,ZGL)
       IMPLICIT REAL (A - Y)
C+SEQ, IMPDP.
#include "necard.h"
#include "pdf804/w50516.inc"
#include "pdf804/expdp.inc"
     +      ZX,ZQ,ZUV,ZDV,ZDEL,ZUDB,ZUB,ZDB,ZSB,ZCB,ZBB,ZTB,ZGL
       REAL X, Q, DURATIO, RUV, RDV
C#ifndef NO_GRV94CORR
*--------> correction by Bodek[hep-ex/0203009](C.Saji Aug. 02) 
*--------> d/u-ratio correction is added(G.Mitsuka, Aug'06)
       if(NEBODEK.eq.1) then
          if(first) print *, '>>> GRV94DI with Bodek corr.'
          X = ZX*(ZQ*ZQ+0.624)/(ZQ*ZQ+1.735*ZX)
          CORR_BODEK=ZQ*ZQ/(ZQ*ZQ+0.188)
          Q2=ZQ*ZQ
          IF(Q2.lt.0.24) then
             Q2=0.24  
             ZQ=sqrt(Q2) 
          endif 
          DURATIO = -0.0161 + 0.0549*ZX + 0.355*ZX**2 - 0.193*ZX**3
C#else
       else if(NEBODEK.eq.0) then
          if(first) print *, '>>> GRV94DI original'
          X = ZX
C#endif
       endif
*<-------
       Q = ZQ
       MU2  = 0.34
       LAM2 = 0.248 * 0.248
       Q2 = Q*Q
       ARG1 = MAX(MU2,Q2)/LAM2
       ARG2 = MU2/LAM2
       S  =  LOG ( LOG( ARG1  ) /  LOG( ARG2   ))
       DS =  SQRT (S)
       S2 = S * S
       S3 = S2 * S
*...UV :
       NU  =  2.484 + 0.116 * S + 0.093 * S2
       AKU =  0.563 - 0.025 * S
       BKU =  0.054 + 0.154 * S
       AU  = -0.326 - 0.058 * S - 0.135 * S2
       BU  = -3.322 + 8.259 * S - 3.119 * S2 + 0.291 * S3
       CU  =  11.52 - 12.99 * S + 3.161 * S2
       DU  =  2.808 + 1.400 * S - 0.557 * S2 + 0.119 * S3
       UV  = GR94FV (X, NU, AKU, BKU, AU, BU, CU, DU)
       ZUV = UV
*...DV :
       ND  =  0.156 - 0.017 * S
       AKD =  0.299 - 0.022 * S
       BKD =  0.259 - 0.015 * S
       AD  =  3.445 + 1.278 * S + 0.326 * S2
       BD  = -6.934 + 37.45 * S - 18.95 * S2 + 1.463 * S3
       CD  =  55.45 - 69.92 * S + 20.78 * S2
       DD  =  3.577 + 1.441 * S - 0.683 * S2 + 0.179 * S3
       DV  = GR94FV (X, ND, AKD, BKD, AD, BD, CD, DD)
       ZDV = DV
*...  DEL :
       NE  =  0.099 + 0.019 * S + 0.002 * S2
       AKE =  0.419 - 0.013 * S
       BKE =  1.064 - 0.038 * S
       AE  = -44.00 + 98.70 * S - 14.79 * S2
       BE  =  28.59 - 40.94 * S - 13.66 * S2 + 2.523 * S3
       CE  =  84.57 - 108.8 * S + 31.52 * S2
       DE  =  7.469 + 2.480 * S - 0.866 * S2
       DEL = GR94FV (X, NE, AKE, BKE, AE, BE, CE, DE)
       ZDEL= DEL
*...UDB :
       ALX =  1.215
       BEX =  0.466
       AKX =  0.326 + 0.150 * S
       BKX =  0.956 + 0.405 * S
       AGX =  0.272
       BGX =  3.794 - 2.359 * DS
       CX  =  2.014
       DX  =  7.941 + 0.534 * DS - 0.940 * S + 0.410 * S2
       EX  =  3.049 + 1.597 * S
       ESX =  4.396 - 4.594 * DS + 3.268 * S
       UDB = GR94FW(X, S, ALX, BEX, AKX, BKX, AGX, BGX, CX, DX, EX, ESX)
       ZUDB= UDB
*...UB
       ZUB = (ZUDB - ZDEL)/2.D0
*...DB
       ZDB = (ZUDB + ZDEL)/2.D0
*...SB :
       ALS =  0.175
       BES =  0.344
       AKS =  1.415 - 0.641 * DS
       AS  =  0.580 - 9.763 * DS + 6.795 * S  - 0.558 * S2
       BS  =  5.617 + 5.709 * DS - 3.972 * S
       DST =  13.78 - 9.581 * S  + 5.370 * S2 - 0.996 * S3
       EST =  4.546 + 0.372 * S2
       ESS =  5.053 - 1.070 * S  + 0.805 * S2
       SB  = GR94FWS (X, S, ALS, BES, AKS, AS, BS, DST, EST, ESS)
       ZSB = SB
*...GL :
       ALG =  1.258
       BEG =  1.846
       AKG =  2.423
       BKG =  2.427 + 1.311 * S  - 0.153 * S2
       AG  =  25.09 - 7.935 * S
       BG  = -14.84 - 124.3 * DS + 72.18 * S
       CG  =  590.3 - 173.8 * S
       DG  =  5.196 + 1.857 * S
       EG  = -1.648 + 3.988 * S  - 0.432 * S2
       ESG =  3.232 - 0.542 * S
       GL  = GR94FW (X, S, ALG, BEG, AKG, BKG, AG, BG, CG, DG, EG, ESG)
       ZGL = GL
*...
       ZCB = 0.D0
       ZBB = 0.D0
       ZTB = 0.D0
*...
C#ifndef NO_GRV94CORR
*--------> correction by Bodek(C.Saji Aug. 02)
*--------> d/u ratio correction is added
       if(NEBODEK.eq.1) then
          RUV=ZUV/(1. + DURATIO*ZUV/(ZUV+ZDV))
          RDV=(ZDV+ZUV*DURATIO)/(1. + DURATIO*ZUV/(ZUV+ZDV))
          ZUV=RUV*CORR_BODEK
          ZDV=RDV*CORR_BODEK
          ZUB=ZUB*CORR_BODEK
          ZDB=ZDB*CORR_BODEK
          ZSB=ZSB*CORR_BODEK
C#else
       else if(NEBODEK.eq.0) then
          ZUV=ZUV
          ZDV=ZDV
          ZUB=ZUB
          ZDB=ZDB
          ZSB=ZSB
C#endif
       endif
*      ZCB=ZCB
*      ZBB=ZBB
*      ZTB=ZTB
       RETURN
       END
