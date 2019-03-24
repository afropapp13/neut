      FUNCTION fnrawcrs( IPAR, MODE, E )

      IMPLICIT NONE

      REAL*4 fnrawcrs

      INTEGER*4 IPAR,MODE
      REAL*4 E

#include "necard.h"

C-- external cross-section functions

      real*4   fnelscrs,fnelscrn,fnelspau_select,fnelspan,cohcrsne
      external fnelscrs,fnelscrn,fnelspau_select,fnelspan,cohcrsne
      real*4   rsspicrs,rssetcrs,rsskcrs,rssgcrs,fnapicrs,fnmec_select
      external rsspicrs,rssetcrs,rsskcrs,rssgcrs,fnapicrs,fnmec_select
      real*4   difcrsnerein
      external difcrsnerein



C-- temporary variables for DIS
      REAL*4 APICRSCC,APICRSNC,R1320CC,R1320NC

C--
      REAL*4 RATN, RATP, RATFRE, H2O

      RATN  =float(NUMBNDN)/float(NUMBNDN+NUMBNDP+NUMFREP)
      RATP  =float(NUMBNDP+NUMFREP)/float(NUMBNDN+NUMBNDP+NUMFREP)
      RATFRE=float(NUMFREP)/float(NUMBNDN+NUMBNDP+NUMFREP)
      H2O   =float(NUMBNDN+NUMBNDP+NUMFREP)

      APICRSCC=FNAPICRS(E,IPAR,1,-1.,R1320CC)
      APICRSNC=FNAPICRS(E,IPAR,0,-1.,R1320NC)
C--
      fnrawcrs = -1.

      if (ipar.gt.0) then

         if (mode.eq.1) then
            fnrawcrs=FNELSPAU_SELECT(E,IPAR)
         endif
         
         if (mode.eq.2) then
            fnrawcrs=FNMEC_SELECT(E,IPAR)
         endif

         if (mode.eq.11) then
            fnrawcrs=RSSPICRS(E,IPAR,1)
         endif
         
         if (mode.eq.12) then
            fnrawcrs=RSSPICRS(E,IPAR,2)
         endif
         
         if (mode.eq.13) then
            fnrawcrs=RSSPICRS(E,IPAR,3)
         endif
         
         if (mode.eq.21) then
            fnrawcrs=APICRSCC*R1320CC
         endif
         
         if (mode.eq.31) then
            fnrawcrs=RSSPICRS(E,IPAR,11)
         endif
         
         if (mode.eq.32) then
            fnrawcrs=RSSPICRS(E,IPAR,12)
         endif
         
         if (mode.eq.33) then
            fnrawcrs=RSSPICRS(E,IPAR,13)
         endif
         
         if (mode.eq.34) then
            fnrawcrs=RSSPICRS(E,IPAR,14)
         endif
         
         if (mode.eq.41) then
            fnrawcrs=APICRSNC*R1320NC
         endif
         
         if (mode.eq.51) then
            fnrawcrs = (
     $           FNELSCRN(E,IPAR,2212)*RATFRE +
     $           FNELSPAN(E,IPAR,2212)*(RATP-RATFRE))/RATP
         endif
         
         if (mode.eq.52) then
            fnrawcrs=FNELSPAN(E,IPAR,2112)
         endif
         
         if (mode.eq.16) then
            fnrawcrs=COHCRSNE(E,IPAR,1)
         endif
         
         if (mode.eq.36) then
            fnrawcrs=COHCRSNE(E,IPAR,0)
         endif
         
         if (mode.eq.15) then
            fnrawcrs=DIFCRSNEREIN(E,IPAR,1)
         endif
         
         if (mode.eq.35) then
            fnrawcrs=DIFCRSNEREIN(E,IPAR,0)
         endif
         
         if (mode.eq.22) then
            fnrawcrs=RSSETCRS(E,IPAR,1)
         endif
         
         if (mode.eq.42) then
            fnrawcrs=RSSETCRS(E,IPAR,11)
         endif
         
         if (mode.eq.43) then
            fnrawcrs=RSSETCRS(E,IPAR,12)
         endif
         
         if (mode.eq.23) then
            fnrawcrs=RSSKCRS(E,IPAR,1)
         endif
         
         if (mode.eq.44) then
            fnrawcrs=RSSKCRS(E,IPAR,11)
         endif
         
         if (mode.eq.45) then
            fnrawcrs=RSSKCRS(E,IPAR,12)
         endif
         
         if (mode.eq.26) then
            fnrawcrs=APICRSCC*(1.-R1320CC)
         endif
         
         if (mode.eq.46) then
            fnrawcrs=APICRSNC*(1.-R1320NC)
         endif
         
         if (mode.eq.17) then
            fnrawcrs=RSSGCRS(E,IPAR, 1)
         endif
         
         if (mode.eq.38) then
            fnrawcrs=RSSGCRS(E,IPAR,11)
         endif
         
         if (mode.eq.39) then
            fnrawcrs=RSSGCRS(E,IPAR,12)
         endif

      endif
      
C---------------------------------------------------------------------
      if (ipar.lt.0) then

         if (mode.eq.-1) then
            fnrawcrs =(
     $           FNELSCRS(E,IPAR)*RATFRE +
     $           FNELSPAU_SELECT(E,IPAR)*(RATP-RATFRE))/RATP
         endif

         if (mode.eq.-2) then
            fnrawcrs = FNMEC_SELECT(E,IPAR)
         endif

         if (mode.eq.-11) then
            fnrawcrs = RSSPICRS(E,IPAR,4)
         endif
            
         if (mode.eq.-12) then
            fnrawcrs = RSSPICRS(E,IPAR,5)
         endif            

         if (mode.eq.-13) then
            fnrawcrs = RSSPICRS(E,IPAR,6)
         endif
         
         if (mode.eq.-21) then
            fnrawcrs = APICRSCC*R1320CC
         endif

         if (mode.eq.-31) then
            fnrawcrs = RSSPICRS(E,IPAR,15)
         endif

         if (mode.eq.-32) then
            fnrawcrs = RSSPICRS(E,IPAR,16)
         endif

         if (mode.eq.-33) then
            fnrawcrs = RSSPICRS(E,IPAR,17)
         endif

         if (mode.eq.-34) then
            fnrawcrs = RSSPICRS(E,IPAR,18)
         endif

         if (mode.eq.-41) then            
            fnrawcrs = APICRSNC*R1320NC
         endif

         if (mode.eq.-51) then                        
             fnrawcrs = (
     $           FNELSCRN(E,IPAR,2212)*RATFRE +
     $           FNELSPAN(E,IPAR,2212)*(RATP-RATFRE))/RATP
         endif

         if (mode.eq.-52) then
            fnrawcrs = FNELSPAN(E,IPAR,2112)
         endif
           
         if (mode.eq.-16) then
            fnrawcrs = COHCRSNE(E,IPAR,1) /H2O
         endif
         
         if (mode.eq.-36) then
            fnrawcrs = COHCRSNE(E,IPAR,0) /H2O
         endif

         if (mode.eq.-15) then
            fnrawcrs=DIFCRSNEREIN(E,IPAR,1)
         endif
         
         if (mode.eq.-35) then
            fnrawcrs=DIFCRSNEREIN(E,IPAR,0)
         endif
         
         if (mode.eq.-22) then
            fnrawcrs = RSSETCRS(E,IPAR,2)
         endif

         if (mode.eq.-42) then
            fnrawcrs = RSSETCRS(E,IPAR,13)
         endif

         if (mode.eq.-43) then
            fnrawcrs = RSSETCRS(E,IPAR,14)
         endif

         if (mode.eq.-23) then
            fnrawcrs = RSSKCRS(E,IPAR,2)
         endif

         if (mode.eq.-44) then
            fnrawcrs = RSSKCRS(E,IPAR,13)
         endif

         if (mode.eq.-45) then
            fnrawcrs = RSSKCRS(E,IPAR,14)
         endif

         if (mode.eq.-26) then
            fnrawcrs = APICRSCC*(1.-R1320CC)
         endif

         if (mode.eq.-46) then
            fnrawcrs = APICRSNC*(1.-R1320NC)
         endif

         if (mode.eq.-17) then
            fnrawcrs = RSSGCRS(E,IPAR, 2)
         endif

         if (mode.eq.-38) then
            fnrawcrs = RSSGCRS(E,IPAR,13)
         endif

         if (mode.eq.-39) then
            fnrawcrs = RSSGCRS(E,IPAR,14)
         endif

      endif

      END
