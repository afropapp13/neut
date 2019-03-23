************************************************************************
*     -------------------------
      REAL FUNCTION FNELSPAU_INTPOL(E,IPAR)
*     -------------------------
*
*     (Purpose)
*       TOTAL CROSS SECTION OF ELASTIC SCATTERING
*       FERMI MOMENTUM AND PAULI PRINCIPLE IS CONSIDERED.
*
*     (Input)
*       E      : NEUTRINO ENERGY ( GEV )
*       IPAR   : NEUTRINO TYPE
*                  12 : NEU-E
*                 -12 : NEU-E-BAR            
*                  14 : NEU-MU
*                 -14 : NEU-MU-BAR            
*                  16 : NEU-TAU
*                 -16 : NEU-TAU-BAR
*     
*     (Output)
*       FNELSPAU_INTPOL : CROSS SECTION ( 10**-38cm^2 )
*
*     (Creation Date and Author)
*       19??.??.?? ; M.Nakahata
*       1987.06.22 ; N.SATO     FOR TAU
*       1988.06.25 ; T.KAJITA   MA 0.95 -> 1.01
*       1988.06.28 ; M.TAKITA   UNIFIED
*       1995.02.16 ; K. KANEYUKI  FOR S.K.
*                    ELSPAU -> FNELSPAU
*                    TEMPORARY VERSION ( FNNUCL = 0. )
*                    CROSS SECTIONS OF NEU_TAU_BAR ARE LARGE AT
*                    HIGH ENERGY REGION, BUT WE ADD NO CORRECTION.
*       2000.09.05 ; K. KANEYUKI FIX FNNUCL, UPDATE CF
*     
C
C       2002.09.17 ; Y.Hayato ; 
C                      Cross-section was calculated with programs in ver 99b.
C                             modifications : dnels.F  
C                                                 line 65 was commented out
C                                                 line 57 do not subtract V
C                                                 fixed potential
C                                             dnelsq2.F
C                                                 MA=1.1 
C
C       2007.11,05 ; G.Mitsuka - cross section for each energy and
C                      target is set after loading cross section table
C                      And rename to FNELSPAU_INTPOL
C       2009.01.19 ; Y.Hayato  - PF table change from XX0 to XX5
************************************************************************
      IMPLICIT NONE

#include "neutparams.h"      
#include "qecrs.h"
#include "neutmodel.h"
#include "necard.h"
#include "neutfilepath.h"

      REAL E
      INTEGER IPAR

      INTEGER LUN
      DATA LUN/80/

      INTEGER I,II,J,JJ,IDUM,ILOADQECRS,ILOADSFQECRS,IQEMA
      REAL RDUM
      CHARACTER*80 DUMSTR
      integer*4      LENSTR
      external       lenstr
      integer*4      lenpath
      SAVE ILOADQECRS
      DATA ILOADQECRS/0/
      SAVE ILOADSFQECRS
      DATA ILOADSFQECRS/0/

      INTEGER*4 IDTAG
      SAVE      IDTAG
      DATA      IDTAG/0/

      INTEGER*4 IDVECFORM

      REAL FERMISURF(9)
      DATA FERMISURF/0.195, 0.205, 0.215, 0.225, 0.235, 
     &               0.245, 0.255, 0.265, 0.275/

      REAL CRSLE1(210),CRSLE2(210),CRSLE3(210),CRSLE4(210)
      REAL CRSLE5(210),CRSLE6(210),CRSLE7(210),CRSLE8(210)
      REAL CRSLE9(210)
      REAL CRSRE1(210),CRSRE2(210),CRSRE3(210),CRSRE4(210)
      REAL CRSRE5(210),CRSRE6(210),CRSRE7(210),CRSRE8(210)
      REAL CRSRE9(210)
      REAL CRSLM1(210),CRSLM2(210),CRSLM3(210),CRSLM4(210)
      REAL CRSLM5(210),CRSLM6(210),CRSLM7(210),CRSLM8(210)
      REAL CRSLM9(210)
      REAL CRSRM1(210),CRSRM2(210),CRSRM3(210),CRSRM4(210)
      REAL CRSRM5(210),CRSRM6(210),CRSRM7(210),CRSRM8(210)
      REAL CRSRM9(210)
      REAL CRSLT1(210),CRSLT2(210),CRSLT3(210),CRSLT4(210)
      REAL CRSLT5(210),CRSLT6(210),CRSLT7(210),CRSLT8(210)
      REAL CRSLT9(210)
      REAL CRSRT1(210),CRSRT2(210),CRSRT3(210),CRSRT4(210)
      REAL CRSRT5(210),CRSRT6(210),CRSRT7(210),CRSRT8(210)
      REAL CRSRT9(210)

C     Enu for e/mu-nu
      REAL EE(210) ,EE1(50) ,EE2(50) ,EE3(50) ,EE4(50), EE5(10)
C     Enu for tau-nu
      REAL ET(125) ,ET1(50) ,ET2(50) ,ET3(25)

      REAL WeightEnu,WeightPF
      REAL CRSLOW,CRSHIG

      INTEGER*4 IDDUM
      REAL*4    PFDUM(9)

      EQUIVALENCE (EE(1),EE1(1))
      EQUIVALENCE (EE(51),EE2(1))
      EQUIVALENCE (EE(101),EE3(1))
      EQUIVALENCE (EE(151),EE4(1))
      EQUIVALENCE (EE(201),EE5(1))

      EQUIVALENCE (ET(1),ET1(1))
      EQUIVALENCE (ET(51),ET2(1))
      EQUIVALENCE (ET(101),ET3(1))

C     
C --  ENERGY OF ELECTRON/MUON NEUTRINO
C
       DATA EE1/0.025,0.075,0.125,0.175,0.225, 
     $          0.275,0.325,0.375,0.425,0.475, 
     $          0.525,0.575,0.625,0.675,0.725, 
     $          0.775,0.825,0.875,0.925,0.975, 
     $          1.025,1.075,1.125,1.175,1.225, 
     $          1.275,1.325,1.375,1.425,1.475, 
     $          1.525,1.575,1.625,1.675,1.725, 
     $          1.775,1.825,1.875,1.925,1.975, 
     $          2.025,2.075,2.125,2.175,2.225, 
     $          2.275,2.325,2.375,2.425,2.475/

      DATA EE2/ 2.525,2.575,2.625,2.675,2.725, 
     $          2.775,2.825,2.875,2.925,2.975, 
     $          3.025,3.075,3.125,3.175,3.225, 
     $          3.275,3.325,3.375,3.425,3.475, 
     $          3.525,3.575,3.625,3.675,3.725, 
     $          3.775,3.825,3.875,3.925,3.975, 
     $          4.025,4.075,4.125,4.175,4.225, 
     $          4.275,4.325,4.375,4.425,4.475, 
     $          4.525,4.575,4.625,4.675,4.725, 
     $          4.775,4.825,4.875,4.925,4.975/

      DATA EE3/ 5.025,5.075,5.125,5.175,5.225, 
     $          5.275,5.325,5.375,5.425,5.475, 
     $          5.525,5.575,5.625,5.675,5.725, 
     $          5.775,5.825,5.875,5.925,5.975, 
     $          6.025,6.075,6.125,6.175,6.225, 
     $          6.275,6.325,6.375,6.425,6.475, 
     $          6.525,6.575,6.625,6.675,6.725, 
     $          6.775,6.825,6.875,6.925,6.975, 
     $          7.025,7.075,7.125,7.175,7.225, 
     $          7.275,7.325,7.375,7.425,7.475/
     
      DATA EE4/ 7.525,7.575,7.625,7.675,7.725, 
     $          7.775,7.825,7.875,7.925,7.975, 
     $          8.025,8.075,8.125,8.175,8.225, 
     $          8.275,8.325,8.375,8.425,8.475, 
     $          8.525,8.575,8.625,8.675,8.725, 
     $          8.775,8.825,8.875,8.925,8.975, 
     $          9.025,9.075,9.125,9.175,9.225, 
     $          9.275,9.325,9.375,9.425,9.475, 
     $          9.525,9.575,9.625,9.675,9.725, 
     $          9.775,9.825,9.875,9.925,9.975/

      DATA EE5/ 11.00,12.00,13.00,14.00,15.00,
     $          16.00,17.00,18.00,19.00,20.00/

C
C --  ENERGY OF TAU NEUTRINO
C      
      DATA ET1 /.270E+01,.280E+01,.290E+01,.300E+01,.310E+01,
     &          .320E+01,.330E+01,.340E+01,.350E+01,.360E+01,
     &          .370E+01,.380E+01,.390E+01,.400E+01,.410E+01,
     &          .420E+01,.430E+01,.440E+01,.450E+01,.460E+01,
     &          .470E+01,.480E+01,.490E+01,.500E+01,.510E+01,
     &          .520E+01,.530E+01,.540E+01,.550E+01,.560E+01,
     &          .570E+01,.580E+01,.590E+01,.600E+01,.610E+01,
     &          .620E+01,.630E+01,.640E+01,.650E+01,.660E+01,
     &          .670E+01,.680E+01,.690E+01,.700E+01,.710E+01,
     &          .720E+01,.730E+01,.740E+01,.750E+01,.760E+01/
      DATA ET2 /.770E+01,.780E+01,.790E+01,.800E+01,.810E+01,
     &          .820E+01,.830E+01,.840E+01,.850E+01,.860E+01,
     &          .870E+01,.880E+01,.890E+01,.900E+01,.910E+01,
     &          .920E+01,.930E+01,.940E+01,.950E+01,.960E+01,
     &          .970E+01,.980E+01,.990E+01,.100E+02,.100E+02,
     &          .110E+02,.120E+02,.130E+02,.140E+02,.150E+02,
     &          .160E+02,.170E+02,.180E+02,.190E+02,.200E+02,
     &          .210E+02,.220E+02,.230E+02,.240E+02,.250E+02,
     &          .260E+02,.270E+02,.280E+02,.290E+02,.300E+02,
     &          .310E+02,.320E+02,.330E+02,.340E+02,.350E+02/
      DATA ET3 /.360E+02,.370E+02,.380E+02,.390E+02,.400E+02,
     &          .410E+02,.420E+02,.430E+02,.440E+02,.450E+02,
     &          .460E+02,.470E+02,.480E+02,.490E+02,.500E+02,
     &          .510E+02,.520E+02,.530E+02,.540E+02,.550E+02,
     &          .560E+02,.570E+02,.580E+02,.590E+02,.600E+02/

      lenpath = lenstr(crstblpath)

      if (MDLQE/100 .EQ. 4 .AND. (numbndp .eq. 6 .or. 
     &    numbndp .eq. 8 .or. numbndp .eq. 26)) then
c      if (MDLQE/100 .EQ. 4) then
        if (iloadsfqecrs.ne.0) goto 300
      else
        if (iloadqecrs.ne.0) goto 100
      endif

      IQEMA = INT(XMAQE*100+0.5)

*******  Spectral function implementation *****
* A.Furmanski - 2013
*
      if (MDLQE/100 .EQ. 4 .AND. (numbndp .eq. 6 .or. 
     &    numbndp .eq. 8 .or. numbndp .eq. 26)) then
c       if (MDLQE/100 .EQ. 4) then
c ----- load xsec table required
        DUMSTR = 'qelSfData/totXsec/'
        if (IPAR .EQ. 14) then
          DUMSTR=DUMSTR(1:18) // '+14'
        else if (IPAR .EQ. -14) then
          DUMSTR=DUMSTR(1:18) // '-14'
        else if (IPAR .EQ. 12) then
          DUMSTR=DUMSTR(1:18) // '+12'
        else if (IPAR .EQ. -12) then
          DUMSTR=DUMSTR(1:18) // '-12'
        else if (IPAR .EQ. 16) then
          DUMSTR=DUMSTR(1:18) // '+16'
        else if (IPAR .EQ. -16) then
          DUMSTR=DUMSTR(1:18) // '-16'
        endif
        
        
        if (NUMBNDP .EQ. 6) then
          DUMSTR=DUMSTR(1:21) // '_1000060120'
        else if (NUMBNDP .EQ. 8) then
          DUMSTR=DUMSTR(1:21) // '_1000080160'
        else if (NUMBNDP .EQ. 26) then
          DUMSTR=DUMSTR(1:21) // '_1000260560'
        endif
        
        DUMSTR=DUMSTR(1:32) // '_cc.csv'
        DUMSTR=CRSTBLPATH(1:lenpath)//DUMSTR
        open(LUN, file=DUMSTR, status="old")

c ----- RDUM here represents the number of entries in the cross-section table
c ----- note - values of EE can be overwritten here if the binning of tables changes in the future
c ----- CRSLE1 is re-used here for storing cross-section values for SF implementation
        read(LUN,*) RDUM
        do i = 1, RDUM
          read(LUN,*) PFDUM(1), CRSLE1(i)
        enddo
        close(LUN)
        write (*,*) "cross section table for sf loaded"
        iloadsfqecrs = 1
300     continue        

c ----- find energy wanted and interpolate to get value of total xsec
        ii = 0
        if (E.le.EE(1)) return ! smaller than energy threshold
        ii = 1
        do i=1,209
          if (E.gt.EE(i) .and. E.le.EE(i+1)) then
            ii = i
            goto 21
          endif
        enddo
21      continue
        WeightEnu = (E-EE(ii))/(EE(ii+1)-EE(ii))    

        if (e.gt.ee(210)) then ! larger than maximum energy - assume linear extrapolation
C          write(*,*) "larger than maximum energy - extrapolating"
          ii = 210
          FNELSPAU_INTPOL = CRSLE1(ii-1)
     $              + (CRSLE1(ii)-CRSLE1(ii-1))*WeightEnu
        else
           FNELSPAU_INTPOL = CRSLE1(ii) 
     $              + (CRSLE1(ii+1)-CRSLE1(ii))*WeightEnu
        endif
c        iloadqecrs = 1
        return
      endif

**************************************

C old code - RFG model implementation
C
C----------Based on Smith-Moniz formula / Code by H.Nakamura
C
      DUMSTR='ccqe_xsec'

      IDVECFORM = mod(MDLQE,10)
      if (IDVECFORM.eq.1) then
         DUMSTR = DUMSTR(1:9) // '_dipole_'
      else if (IDVECFORM.eq.2) then
         DUMSTR = DUMSTR(1:9) // '_bbba05_'
      elseif (IDVECFORM.eq.3) then
         DUMSTR = DUMSTR(1:9) // '_bbba07_'
      else
         write(*,*) 'fnelspau_interpol: ERROR MODEL for QE', MDLQE,
     $        ' IS NOT DEFINED'
         STOP
      endif         

      if (IQEMA.eq.101) then
         DUMSTR = DUMSTR(1:17) // 'ma1.01'
         print*, 'Loading Cross section table for QE : MA=1.01'
         IDTAG=101
      else if (IQEMA.eq.103) then
         DUMSTR = DUMSTR(1:17) // 'ma1.03'
         print*, 'Loading Cross section table for QE : MA=1.03'
         IDTAG=103
      else if (IQEMA.eq.111) then
         DUMSTR = DUMSTR(1:17) // 'ma1.11'
         print*, 'Loading Cross section table for QE : MA=1.11'
         IDTAG=111
      else if (IQEMA.eq.121) then
         DUMSTR = DUMSTR(1:17) // 'ma1.21'
         print*, 'Loading Cross section table for QE : MA=1.21'
         IDTAG=121
      else if (IQEMA.eq.131) then
         DUMSTR = DUMSTR(1:17) // 'ma1.31'
         print*, 'Loading Cross section table for QE : MA=1.31'
         IDTAG=131
      else if (IQEMA.eq.160) then
         DUMSTR = DUMSTR(1:17) // 'ma1.60'
         print*, 'Loading Cross section table for QE : MA=1.60'
         IDTAG=160
      else
         write(*,*) 'fnelspau_interpol: ERROR MA FOR QE', XMAQE,
     $        ' IS NOT DEFINED'
         STOP
      end if


      if (MDLQE/100.eq.1) THEN
         DUMSTR = DUMSTR(1:23) // '_trcor.dat'
         print*, ' with Bodek transverse enhancement'
      else
         DUMSTR = DUMSTR(1:23) // '.dat'
      end if
      
      DUMSTR=CRSTBLPATH(1:lenpath)//DUMSTR
      print*, ' from file: ', DUMSTR(1:lenstr(dumstr))
      open(LUN,file=DUMSTR,form='formatted',status='old')               
      
      read(LUN,*) IDDUM,RDUM,
     $            PFDUM(1),PFDUM(2),PFDUM(3),PFDUM(4),PFDUM(5),
     $            PFDUM(6),PFDUM(7),PFDUM(8),PFDUM(9)

      if (IDDUM.ne.IDTAG) then
         write(*,*) 'This QE cross-section table is not compatible.'
         stop
      endif

      do 200 i = 1, 9
         if (abs(PFDUM(i)-FERMISURF(i)).gt.0.001) then
            write(*,*) 'This QE cross-section table is not compatible.'
            stop
         endif
 200  continue

      do i = 1,1090
         if(i.le.210) then ! nue
            read(LUN,*) IDUM,RDUM,
     &           CRSLE1(i),CRSLE2(i),CRSLE3(i),CRSLE4(i),
     &           CRSLE5(i),CRSLE6(i),CRSLE7(i),CRSLE8(i),CRSLE9(i)
C            print*, 'debug:CRSLE1',i,CRSLE1(i),CRSLE2(i),CRSLE3(i)
            
            if ( abs(RDUM-EE(I)).gt.0.001 ) then
               write(*,*) 'Cross-section table seems to be corrupted.'
               stop
            endif

         elseif(i.gt.210 .and. i.le.420) then ! nue-bar
            read(LUN,*) IDUM,RDUM,
     &           CRSRE1(i-210),CRSRE2(i-210),CRSRE3(i-210),
     $           CRSRE4(i-210),CRSRE5(i-210),CRSRE6(i-210),
     $           CRSRE7(i-210),CRSRE8(i-210),CRSRE9(i-210)
C            print*, 'debug:CRSRE1',i,CRSRE1(i),CRSRE2(i),CRSRE3(i)

            if ( abs(RDUM-EE(I-210)).gt.0.001 ) then
               write(*,*) 'Cross-section table seems to be corrupted.'
               stop
            endif

         elseif(i.gt.420 .and. i.le.630) then ! numu
            read(LUN,*) IDUM,RDUM,
     &           CRSLM1(i-420),CRSLM2(i-420),CRSLM3(i-420),
     $           CRSLM4(i-420),CRSLM5(i-420),CRSLM6(i-420),
     $           CRSLM7(i-420),CRSLM8(i-420),CRSLM9(i-420)
C            print*, 'debug:CRSLM1',i,CRSLM1(i),CRSLM2(i),CRSLM3(i)

            if ( abs(RDUM-EE(I-420)).gt.0.001 ) then
               write(*,*) 'Cross-section table seems to be corrupted.'
               stop
            endif

         elseif(i.gt.630 .and. i.le.840) then ! numu-bar
            read(LUN,*) IDUM,RDUM,
     &           CRSRM1(i-630),CRSRM2(i-630),CRSRM3(i-630),
     $           CRSRM4(i-630),CRSRM5(i-630),CRSRM6(i-630),
     $           CRSRM7(i-630),CRSRM8(i-630),CRSRM9(i-630)
C            print*, 'debug:CRSRM1',i,CRSRM1(i),CRSRM2(i),CRSRM3(i)

            if ( abs(RDUM-EE(I-630)).gt.0.001 ) then
               write(*,*) 'Cross-section table seems to be corrupted.'
               stop
            endif

         elseif(i.gt.840 .and. i.le.965) then ! nutau
            read(LUN,*) IDUM,RDUM,
     &           CRSLT1(i-840),CRSLT2(i-840),CRSLT3(i-840),
     $           CRSLT4(i-840),CRSLT5(i-840),CRSLT6(i-840),
     $           CRSLT7(i-840),CRSLT8(i-840),CRSLT9(i-840)
C            print*, 'debug:CRSLT1',i,CRSLT1(i),CRSLT2(i),CRSLT3(i)

            if ( abs(RDUM-ET(I-840)).gt.0.001 ) then
               write(*,*) 'Cross-section table seems to be corrupted.'
               stop
            endif

         elseif(i.gt.965 .and. i.le.1090) then ! nutau-bar
            read(LUN,*) IDUM,RDUM,
     &           CRSRT1(i-965),CRSRT2(i-965),CRSRT3(i-965),
     $           CRSRT4(i-965),CRSRT5(i-965),CRSRT6(i-965),
     $           CRSRT7(i-965),CRSRT8(i-965),CRSRT9(i-965)
C            print*, 'debug:CRSRT1',i,CRSRT1(i),CRSRT2(i),CRSRT3(i)

            if ( abs(RDUM-ET(I-965)).gt.0.001 ) then
               write(*,*) 'Cross-section table seems to be corrupted.'
               stop
            endif

         endif
      enddo

      do i=1,210
*     nue
         CRSRFGLE(1,i)=CRSLE1(i)
         CRSRFGLE(2,i)=CRSLE2(i)
         CRSRFGLE(3,i)=CRSLE3(i)
         CRSRFGLE(4,i)=CRSLE4(i)
         CRSRFGLE(5,i)=CRSLE5(i)
         CRSRFGLE(6,i)=CRSLE6(i)
         CRSRFGLE(7,i)=CRSLE7(i)
         CRSRFGLE(8,i)=CRSLE8(i)
         CRSRFGLE(9,i)=CRSLE9(i)
*     nue-bar
         CRSRFGRE(1,i)=CRSRE1(i)
         CRSRFGRE(2,i)=CRSRE2(i)
         CRSRFGRE(3,i)=CRSRE3(i)
         CRSRFGRE(4,i)=CRSRE4(i)
         CRSRFGRE(5,i)=CRSRE5(i)
         CRSRFGRE(6,i)=CRSRE6(i)
         CRSRFGRE(7,i)=CRSRE7(i)
         CRSRFGRE(8,i)=CRSRE8(i)
         CRSRFGRE(9,i)=CRSRE9(i)
*     numu
         CRSRFGLM(1,i)=CRSLM1(i)
         CRSRFGLM(2,i)=CRSLM2(i)
         CRSRFGLM(3,i)=CRSLM3(i)
         CRSRFGLM(4,i)=CRSLM4(i)
         CRSRFGLM(5,i)=CRSLM5(i)
         CRSRFGLM(6,i)=CRSLM6(i)
         CRSRFGLM(7,i)=CRSLM7(i)
         CRSRFGLM(8,i)=CRSLM8(i)
         CRSRFGLM(9,i)=CRSLM9(i)
*     numu-bar
         CRSRFGRM(1,i)=CRSRM1(i)
         CRSRFGRM(2,i)=CRSRM2(i)
         CRSRFGRM(3,i)=CRSRM3(i)
         CRSRFGRM(4,i)=CRSRM4(i)
         CRSRFGRM(5,i)=CRSRM5(i)
         CRSRFGRM(6,i)=CRSRM6(i)
         CRSRFGRM(7,i)=CRSRM7(i)
         CRSRFGRM(8,i)=CRSRM8(i)
         CRSRFGRM(9,i)=CRSRM9(i)
      enddo
      do i=1,125
*     nutau
         CRSRFGLT(1,i)=CRSLT1(i)
         CRSRFGLT(2,i)=CRSLT2(i)
         CRSRFGLT(3,i)=CRSLT3(i)
         CRSRFGLT(4,i)=CRSLT4(i)
         CRSRFGLT(5,i)=CRSLT5(i)
         CRSRFGLT(6,i)=CRSLT6(i)
         CRSRFGLT(7,i)=CRSLT7(i)
         CRSRFGLT(8,i)=CRSLT8(i)
         CRSRFGLT(9,i)=CRSLT9(i)
*     nutau-bar
         CRSRFGRT(1,i)=CRSRT1(i)
         CRSRFGRT(2,i)=CRSRT2(i)
         CRSRFGRT(3,i)=CRSRT3(i)
         CRSRFGRT(4,i)=CRSRT4(i)
         CRSRFGRT(5,i)=CRSRT5(i)
         CRSRFGRT(6,i)=CRSRT6(i)
         CRSRFGRT(7,i)=CRSRT7(i)
         CRSRFGRT(8,i)=CRSRT8(i)
         CRSRFGRT(9,i)=CRSRT9(i)
      enddo

      close(LUN)
      iloadqecrs = 1
      write(*,*) 'Completed CCQE cross-section loading.'
 100  continue

      FNELSPAU_INTPOL=0.
      WeightEnu=0.
      WeightPF=0.
      ii=0
      jj=0
      CRSLOW=0.
      CRSHIG=0.

C     set weight for interpolation of fermi surface momentum
      if( (PFSURF.le.FERMISURF(1)) .or. 
     $     (PFSURF.gt.FERMISURF(9))) then
         return
      endif
      do 10 j=1,8
         if( (PFSURF.gt.FERMISURF(j) )
     $        .and. 
     $        (PFSURF.le.FERMISURF(j+1))) then
            jj = j
            goto 10
         endif
 10   continue
      WeightPF 
     $     = (PFSURF-FERMISURF(jj))
     $     /(FERMISURF(jj+1)-FERMISURF(jj))

C
C ++ FOR E-NEUTRINO
C
*     set weight for interpolation of neutrino energy
      if (iabs(ipar).eq.12) then
         if (e.le.ee(1)) return ! smaller than energy threshold
         ii = 1
         do 11 i=1,209
            if (e.gt.ee(i) .and. e.le.ee(i+1)) then
               ii = i
               goto 11
            endif
 11      continue
         WeightEnu = (E-EE(ii))/(EE(ii+1)-EE(ii))

         if (ipar.eq.12) then
            if (e.gt.ee(210)) then ! larger than maximum energy
               ii = 210
               CRSLOW   = CRSRFGLE(jj  ,ii)
               CRSHIG   = CRSRFGLE(jj+1,ii)
               FNELSPAU_INTPOL = CRSLOW + (CRSHIG-CRSLOW)*WeightPF
            else
               CRSLOW   = CRSRFGLE(jj  ,ii) 
     $              + (CRSRFGLE(jj  ,ii+1)-CRSRFGLE(jj  ,ii))*WeightEnu
               CRSHIG   = CRSRFGLE(jj+1,ii) 
     $              + (CRSRFGLE(jj+1,ii+1)-CRSRFGLE(jj+1,ii))*WeightEnu
               FNELSPAU_INTPOL = CRSLOW + (CRSHIG-CRSLOW)*WeightPF
            endif
         else
            if (e.gt.ee(210)) then ! larger than maximum energy
               ii = 210
               CRSLOW   = CRSRFGRE(jj  ,ii)
               CRSHIG   = CRSRFGRE(jj+1,ii)
               FNELSPAU_INTPOL = CRSLOW + (CRSHIG-CRSLOW)*WeightPF
            else
               CRSLOW   = CRSRFGRE(jj  ,ii) 
     $              + (CRSRFGRE(jj  ,ii+1)-CRSRFGRE(jj  ,ii))*WeightEnu
               CRSHIG   = CRSRFGRE(jj+1,ii) 
     $              + (CRSRFGRE(jj+1,ii+1)-CRSRFGRE(jj+1,ii))*WeightEnu
               FNELSPAU_INTPOL = CRSLOW + (CRSHIG-CRSLOW)*WeightPF
            endif
         endif
         return
      endif

C
C ++ FOR MU-NEUTRINO
C
*     set weight for interpolation of neutrino energy
      if (iabs(ipar).eq.14) then
         if (e.le.ee(1)) return ! smaller than energy threshold
         ii = 1
         do 12 i=1,209
            if (e.gt.ee(i) .and. e.le.ee(i+1)) then
               ii = i
               goto 12
            endif
 12      continue
         WeightEnu = (E-EE(ii))/(EE(ii+1)-EE(ii))

         if (ipar.eq.14) then
            if (e.gt.ee(210)) then ! larger than maximum energy
               ii = 210
               CRSLOW   = CRSRFGLM(jj  ,ii)
               CRSHIG   = CRSRFGLM(jj+1,ii)
               FNELSPAU_INTPOL = CRSLOW + (CRSHIG-CRSLOW)*WeightPF
            else
               CRSLOW   = CRSRFGLM(jj  ,ii) 
     $              + (CRSRFGLM(jj  ,ii+1)-CRSRFGLM(jj  ,ii))*WeightEnu
               CRSHIG   = CRSRFGLM(jj+1,ii) 
     $              + (CRSRFGLM(jj+1,ii+1)-CRSRFGLM(jj+1,ii))*WeightEnu
               FNELSPAU_INTPOL = CRSLOW + (CRSHIG-CRSLOW)*WeightPF
            endif
         else
            if (e.gt.ee(210)) then ! larger than maximum energy
               ii = 210
               CRSLOW   = CRSRFGRM(jj  ,ii)
               CRSHIG   = CRSRFGRM(jj+1,ii)
               FNELSPAU_INTPOL = CRSLOW + (CRSHIG-CRSLOW)*WeightPF
            else
               CRSLOW   = CRSRFGRM(jj  ,ii) 
     $              + (CRSRFGRM(jj  ,ii+1)-CRSRFGRM(jj  ,ii))*WeightEnu
               CRSHIG   = CRSRFGRM(jj+1,ii) 
     $              + (CRSRFGRM(jj+1,ii+1)-CRSRFGRM(jj+1,ii))*WeightEnu
               FNELSPAU_INTPOL = CRSLOW + (CRSHIG-CRSLOW)*WeightPF
            endif
         endif
         return
      endif

C
C ++ FOR TAU-NEUTRINO
C
*     set weight for interpolation of neutrino energy
      if (iabs(ipar).eq.16) then
         if (e.le.et(1)) return ! smaller than energy threshold
         ii = 1
         do 13 i=1,124
            if (e.gt.et(i) .and. e.le.et(i+1)) then
               ii = i
               goto 13
            endif
 13      continue
         WeightEnu = (E-ET(ii))/(ET(ii+1)-ET(ii))

         if (ipar.eq.16) then
            if (e.gt.et(125)) then ! larger than maximum energy
               ii = 125
               CRSLOW   = CRSRFGLT(jj  ,ii)
               CRSHIG   = CRSRFGLT(jj+1,ii)
               FNELSPAU_INTPOL = CRSLOW + (CRSHIG-CRSLOW)*WeightPF
            else
               CRSLOW   = CRSRFGLT(jj  ,ii) 
     $              + (CRSRFGLT(jj  ,ii+1)-CRSRFGLT(jj  ,ii))*WeightEnu
               CRSHIG   = CRSRFGLT(jj+1,ii) 
     $              + (CRSRFGLT(jj+1,ii+1)-CRSRFGLT(jj+1,ii))*WeightEnu
               FNELSPAU_INTPOL = CRSLOW + (CRSHIG-CRSLOW)*WeightPF
            endif
         else
            if (e.gt.et(125)) then ! larger than maximum energy
               ii = 125
               CRSLOW   = CRSRFGRT(jj  ,ii)
               CRSHIG   = CRSRFGRT(jj+1,ii)
               FNELSPAU_INTPOL = CRSLOW + (CRSHIG-CRSLOW)*WeightPF
            else
               CRSLOW   = CRSRFGRT(jj  ,ii) 
     $              + (CRSRFGRT(jj  ,ii+1)-CRSRFGRT(jj  ,ii))*WeightEnu
               CRSHIG   = CRSRFGRT(jj+1,ii) 
     $              + (CRSRFGRT(jj+1,ii+1)-CRSRFGRT(jj+1,ii))*WeightEnu
               FNELSPAU_INTPOL = CRSLOW + (CRSHIG-CRSLOW)*WeightPF
            endif
         endif
         return
      endif

      RETURN
      END
