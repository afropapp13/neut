************************************************************************
*     ------------------------
      subroutine dumpcrs(IPAR,E)
*     ------------------------
*
*     (Purpose)
*       Calc. TOTAL CROSS SECTION OF NEUTRINO
*
*     (Input)
*       IPAR   : NEUTRINO TYPE
*                  12 : NEU-E
*                 -12 : NEU-E-BAR            
*                  14 : NEU-MU
*                 -14 : NEU-MU-BAR            
*                  16 : NEU-TAU
*                 -16 : NEU-TAU-BAR
*       E      : NEUTRINO ENERGY ( GEV )
*     
*     (Output)
*       FNTOTPAU : CROSS SECTION ( 10^-38cm^2 )
*     
************************************************************************
C
C    MEANING OF DATA
C    RATN   : NEUTRON RATIO IN TARGET
C    RATP   : PROTON RATIO IN TARGET
C    RATFRE : FREE PROTON RATIO IN TARGET
C     
      IMPLICIT NONE

      integer*4 IPAR
      real*4    E

      real*4 apicrscc, apicrsnc

      real*4 xccqen
      real*4 xccqef,xccqep

      real*4 xncelf,xncelp,xnceln

      real*4 xccspi1,xccspi2,xccspi3
      real*4 xncspi1,xncspi2,xncspi3,xncspi4

      real*4 xcccoh,xnccoh

      real*4 xccmpi1,xccmpi2
      real*4 xncmpi1,xncmpi2

      real*4 xccseta,xncseta1,xncseta2
      real*4 xccsk,xncsk1,xncsk2
      real*4 xccsg,xncsg1,xncsg2

      real*4 fnapicrs,fnelscrs,rsspicrs,fnelspau_select
      real*4 fnelscrn,fnelspan,rssetcrs,rsskcrs,rssgcrs
      real*4 cohcrsne

#include "necard.h"

      real*4 etmp,r1320cc,r1320nc

      real*4 ratn,ratp,ratfre
      REAL MOLECULE
      
#ifdef UPMU_RCK
      DATA RATN/0.5/
      DATA RATP/0.5/
      DATA RATFRE/0./
      MOLECULE=18.0
#else
C     NUMBNDN, NUMBNDN etc. are set in necard (2007 Nov G.Mitsuka)
      RATN  =float(NUMBNDN)/float(NUMBNDN+NUMBNDP+NUMFREP)
      RATP  =float(NUMBNDP+NUMFREP)/float(NUMBNDN+NUMBNDP+NUMFREP)
      RATFRE=float(NUMFREP)/float(NUMBNDN+NUMBNDP+NUMFREP)
      MOLECULE=float(NUMBNDN+NUMBNDP+NUMFREP)
#endif

C     --ENERGY UPPER LIMIT FOR 1PI

      IF(E.LE.1000.) THEN
         ETMP=E
      ELSE
         ETMP=999.999
      ENDIF

C     
      IF (IPAR.GT.0) THEN
C     
C     -- FOR NEUTRINO
C     
#ifdef UPMU_RCK
         APICRSCC=FNAPICRS(E,IPAR,1,RATP,R1320CC)
         APICRSNC=FNAPICRS(E,IPAR,0,RATN,R1320NC)
#else
         APICRSCC=FNAPICRS(E,IPAR,1,-1.,R1320CC)
         APICRSNC=FNAPICRS(E,IPAR,0,-1.,R1320NC)
#endif
         XCCQEN=FNELSPAU_SELECT(E,IPAR)

         XCCSPI1=RSSPICRS(ETMP,IPAR,1)
         XCCSPI2=RSSPICRS(ETMP,IPAR,2)
         XCCSPI3=RSSPICRS(ETMP,IPAR,3)

         XCCMPI1=APICRSCC*R1320CC

         XNCSPI1=RSSPICRS(ETMP,IPAR,11)
         XNCSPI2=RSSPICRS(ETMP,IPAR,12)
         XNCSPI3=RSSPICRS(ETMP,IPAR,13)
         XNCSPI4=RSSPICRS(ETMP,IPAR,14)

         XNCMPI1=APICRSNC*R1320NC

         XNCELF=FNELSCRN(E,IPAR,2212)
         XNCELP=FNELSPAN(E,IPAR,2212)
         XNCELN=FNELSPAN(E,IPAR,2112)

         XCCCOH=COHCRSNE(E,IPAR,1) /MOLECULE
         XNCCOH=COHCRSNE(E,IPAR,0) /MOLECULE
         
         XCCSETA=RSSETCRS(ETMP,IPAR,1)
         XNCSETA1=RSSETCRS(ETMP,IPAR,11)
         XNCSETA2=RSSETCRS(ETMP,IPAR,12)

         XCCSK=RSSKCRS(ETMP,IPAR,1)
         XNCSK1=RSSKCRS(ETMP,IPAR,11)
         XNCSK2=RSSKCRS(ETMP,IPAR,12)
         
         XCCMPI2=APICRSCC*(1.-R1320CC)
         XNCMPI2=APICRSNC*(1.-R1320NC)

         XCCSG=RSSGCRS(E,IPAR, 1)
         XNCSG1=RSSGCRS(E,IPAR,11)
         XNCSG2=RSSGCRS(E,IPAR,12)
         
         write(*,10) ipar,E,
     $        0.00,xccqen,
     $        xccspi1,xccspi2,xccspi3,
     $        xcccoh,
     $        xccmpi1,xccmpi2,
     $        xccseta,xccsk,xccsg,
     $        xncelf,xncelp,xnceln,
     $        xncspi1,xncspi2,xncspi3,xncspi4,
     $        xnccoh,xncmpi1,xncmpi2,
     $        xncseta1,xncseta2,xncsk1,xncsk2,xncsg1,xncsg2

 10      format(I4,F8.3,F8.5,F8.5,F8.5,F8.5,F8.5,F8.5,F8.5,F12.5,
     $          F8.5,F8.5,F8.5,F8.5,F8.5,F8.5,F8.5,F8.5,F8.5,F8.5,
     $          F8.5,F8.5,F12.5,F8.5,F8.5,F8.5,F8.5,F8.5,F8.5)

      ELSE
C     
C     -- FOR ANTI NEUTRINO
C     
#ifdef UPMU_RCK
         APICRSCC=FNAPICRS(E,IPAR,1,RATP,R1320CC)
         APICRSNC=FNAPICRS(E,IPAR,0,RATN,R1320NC)
#else
         APICRSCC=FNAPICRS(E,IPAR,1,-1.,R1320CC)
         APICRSNC=FNAPICRS(E,IPAR,0,-1.,R1320NC)
#endif
         XCCQEF=FNELSCRS(E,IPAR)
         
         XCCSPI1=RSSPICRS(ETMP,IPAR,4)
         XCCSPI2=RSSPICRS(ETMP,IPAR,5)
         XCCSPI3=RSSPICRS(ETMP,IPAR,6)
         
         XCCMPI1=APICRSCC*R1320CC

         XNCSPI1=RSSPICRS(ETMP,IPAR,15)
         XNCSPI2=RSSPICRS(ETMP,IPAR,16)
         XNCSPI3=RSSPICRS(ETMP,IPAR,17)
         XNCSPI4=RSSPICRS(ETMP,IPAR,18)*RATP*CRSNEUTB(9)

         XNCMPI1=APICRSNC*R1320NC
         
         XCCQEP=FNELSPAU_SELECT(E,IPAR)
         
         XNCELF=FNELSCRN(E,IPAR,2212)
         XNCELP=FNELSPAN(E,IPAR,2212)
         XNCELN=FNELSPAN(E,IPAR,2112)
         
         XCCCOH=COHCRSNE(E,IPAR,1) /MOLECULE
         XNCCOH=COHCRSNE(E,IPAR,0) /MOLECULE

         XCCSETA=RSSETCRS(ETMP,IPAR,2)
         XNCSETA1=RSSETCRS(ETMP,IPAR,13)
         XNCSETA2=RSSETCRS(ETMP,IPAR,14)

         XCCSK=RSSKCRS(ETMP,IPAR,2)
         XNCSK1=RSSKCRS(ETMP,IPAR,13)
         XNCSK2=RSSKCRS(ETMP,IPAR,14)
         
         XCCMPI2=APICRSCC*(1.-R1320CC)
         XNCMPI2=APICRSNC*(1.-R1320NC)

         XCCSG=RSSGCRS(E,IPAR, 2)
         XNCSG1=RSSGCRS(E,IPAR,13)
         XNCSG2=RSSGCRS(E,IPAR,14)
         
         write(*,20)  ipar,E,
     $        xccqef,xccqep,
     $        xccspi1,xccspi2,xccspi3,
     $        xcccoh,
     $        xccmpi1,xccmpi2,
     $        xccseta,xccsk,xccsg,
     $        xncelf,xncelp,xnceln,
     $        xncspi1,xncspi2,xncspi3,xncspi4,
     $        xnccoh,
     $        xncmpi1,xncmpi2,
     $        xncseta1,xncseta2,xncsk1,xncsk2,xncsg1,xncsg2

 20      format(I4,F8.3,F8.5,F8.5,F8.5,F8.5,F8.5,F8.5,F8.5,F12.5,
     $          F8.5,F8.5,F8.5,F8.5,F8.5,F8.5,F8.5,F8.5,F8.5,F8.5,
     $          F8.5,F8.5,F12.5,F8.5,F8.5,F8.5,F8.5,F8.5,F8.5)

      ENDIF
      RETURN
      END
