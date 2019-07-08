************************************************************************
*     -------------------------
      REAL FUNCTION FNELSPAU_EFFSF(E,IPAR)
*     -------------------------
*     
*     (Purpose)
*     TOTAL CROSS SECTION OF ELASTIC SCATTERING
*     FERMI MOMENTUM AND PAULI PRINCIPLE IS CONSIDERED.
*     
*     (Input)
*     E      : NEUTRINO ENERGY ( GEV )
*     IPAR   : NEUTRINO TYPE
*     12 : NEU-E
*     -12 : NEU-E-BAR            
*     14 : NEU-MU
*     -14 : NEU-MU-BAR            
*     16 : NEU-TAU
*     -16 : NEU-TAU-BAR
*     
*     (Output)
*     FNELSPAU_INTPOL : CROSS SECTION ( 10**-38cm^2 )
*     
*     (Creation Date and Author)
*     19??.??.?? ; M.Nakahata
*     1987.06.22 ; N.SATO     FOR TAU
*     1988.06.25 ; T.KAJITA   MA 0.95 -> 1.01
*     1988.06.28 ; M.TAKITA   UNIFIED
*     1995.02.16 ; K. KANEYUKI  FOR S.K.
*     ELSPAU -> FNELSPAU
*     TEMPORARY VERSION ( FNNUCL = 0. )
*     CROSS SECTIONS OF NEU_TAU_BAR ARE LARGE AT
*     HIGH ENERGY REGION, BUT WE ADD NO CORRECTION.
*     2000.09.05 ; K. KANEYUKI FIX FNNUCL, UPDATE CF
*     
C     
C     2002.09.17 ; Y.Hayato ; 
C     Cross-section was calculated with programs in ver 99b.
C     modifications : dnels.F  
C     line 65 was commented out
C     line 57 do not subtract V
C     fixed potential
C     dnelsq2.F
C     MA=1.1 
C     
C     2007.11,05 ; G.Mitsuka - cross section for each energy and
C     target is set after loading cross section table
C     And rename to FNELSPAU_INTPOL
C     2009.01.19 ; Y.Hayato  - PF table change from XX0 to XX5
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

      INTEGER I,II,J,JJ,jj1,IDUM,ILOADEFFSFQECRS
      REAL RDUM
      CHARACTER*1024 DUMSTR
      CHARACTER*4    QEMASTR
      integer*4      LENSTR
      external       lenstr
      integer*4      lenpath
      SAVE ILOADEFFSFQECRS
      DATA ILOADEFFSFQECRS/0/

      INTEGER*4 IDTAG
      SAVE      IDTAG
      DATA      IDTAG/0/

      INTEGER*4 IDVECFORM,IDX_TARG

      REAL FERMISURF(9)
      DATA FERMISURF/0.195, 0.205, 0.215, 0.225, 0.235, 
     &     0.245, 0.255, 0.265, 0.275/

C     Enu for e/mu-nu
      REAL EE(210) ,EE1(50) ,EE2(50) ,EE3(50) ,EE4(50), EE5(10)
C     Enu for tau-nu
      REAL ET(125) ,ET1(50) ,ET2(50) ,ET3(25)

      REAL WeightEnu,WeightPF
      REAL CRSSFLOW,CRSSFHIG

      INTEGER*4 IDDUM
      REAL*4    PFDUM(9)

      real*4    fnelspau_intpol
      external  fnelspau_intpol


      EQUIVALENCE (EE(1),EE1(1))
      EQUIVALENCE (EE(51),EE2(1))
      EQUIVALENCE (EE(101),EE3(1))
      EQUIVALENCE (EE(151),EE4(1))
      EQUIVALENCE (EE(201),EE5(1))

      EQUIVALENCE (ET(1),ET1(1))
      EQUIVALENCE (ET(51),ET2(1))
      EQUIVALENCE (ET(101),ET3(1))

C     
C     --  ENERGY OF ELECTRON/MUON NEUTRINO
C     
      DATA EE1/0.025,0.075,0.125,0.175,0.225, 
     $     0.275,0.325,0.375,0.425,0.475, 
     $     0.525,0.575,0.625,0.675,0.725, 
     $     0.775,0.825,0.875,0.925,0.975, 
     $     1.025,1.075,1.125,1.175,1.225, 
     $     1.275,1.325,1.375,1.425,1.475, 
     $     1.525,1.575,1.625,1.675,1.725, 
     $     1.775,1.825,1.875,1.925,1.975, 
     $     2.025,2.075,2.125,2.175,2.225, 
     $     2.275,2.325,2.375,2.425,2.475/

      DATA EE2/ 2.525,2.575,2.625,2.675,2.725, 
     $     2.775,2.825,2.875,2.925,2.975, 
     $     3.025,3.075,3.125,3.175,3.225, 
     $     3.275,3.325,3.375,3.425,3.475, 
     $     3.525,3.575,3.625,3.675,3.725, 
     $     3.775,3.825,3.875,3.925,3.975, 
     $     4.025,4.075,4.125,4.175,4.225, 
     $     4.275,4.325,4.375,4.425,4.475, 
     $     4.525,4.575,4.625,4.675,4.725, 
     $     4.775,4.825,4.875,4.925,4.975/

      DATA EE3/ 5.025,5.075,5.125,5.175,5.225, 
     $     5.275,5.325,5.375,5.425,5.475, 
     $     5.525,5.575,5.625,5.675,5.725, 
     $     5.775,5.825,5.875,5.925,5.975, 
     $     6.025,6.075,6.125,6.175,6.225, 
     $     6.275,6.325,6.375,6.425,6.475, 
     $     6.525,6.575,6.625,6.675,6.725, 
     $     6.775,6.825,6.875,6.925,6.975, 
     $     7.025,7.075,7.125,7.175,7.225, 
     $     7.275,7.325,7.375,7.425,7.475/
      
      DATA EE4/ 7.525,7.575,7.625,7.675,7.725, 
     $     7.775,7.825,7.875,7.925,7.975, 
     $     8.025,8.075,8.125,8.175,8.225, 
     $     8.275,8.325,8.375,8.425,8.475, 
     $     8.525,8.575,8.625,8.675,8.725, 
     $     8.775,8.825,8.875,8.925,8.975, 
     $     9.025,9.075,9.125,9.175,9.225, 
     $     9.275,9.325,9.375,9.425,9.475, 
     $     9.525,9.575,9.625,9.675,9.725, 
     $     9.775,9.825,9.875,9.925,9.975/

      DATA EE5/ 11.00,12.00,13.00,14.00,15.00,
     $     16.00,17.00,18.00,19.00,20.00/

C     
C     --  ENERGY OF TAU NEUTRINO
C     
      DATA ET1 /.270E+01,.280E+01,.290E+01,.300E+01,.310E+01,
     &     .320E+01,.330E+01,.340E+01,.350E+01,.360E+01,
     &     .370E+01,.380E+01,.390E+01,.400E+01,.410E+01,
     &     .420E+01,.430E+01,.440E+01,.450E+01,.460E+01,
     &     .470E+01,.480E+01,.490E+01,.500E+01,.510E+01,
     &     .520E+01,.530E+01,.540E+01,.550E+01,.560E+01,
     &     .570E+01,.580E+01,.590E+01,.600E+01,.610E+01,
     &     .620E+01,.630E+01,.640E+01,.650E+01,.660E+01,
     &     .670E+01,.680E+01,.690E+01,.700E+01,.710E+01,
     &     .720E+01,.730E+01,.740E+01,.750E+01,.760E+01/
      DATA ET2 /.770E+01,.780E+01,.790E+01,.800E+01,.810E+01,
     &     .820E+01,.830E+01,.840E+01,.850E+01,.860E+01,
     &     .870E+01,.880E+01,.890E+01,.900E+01,.910E+01,
     &     .920E+01,.930E+01,.940E+01,.950E+01,.960E+01,
     &     .970E+01,.980E+01,.990E+01,.100E+02,.100E+02,
     &     .110E+02,.120E+02,.130E+02,.140E+02,.150E+02,
     &     .160E+02,.170E+02,.180E+02,.190E+02,.200E+02,
     &     .210E+02,.220E+02,.230E+02,.240E+02,.250E+02,
     &     .260E+02,.270E+02,.280E+02,.290E+02,.300E+02,
     &     .310E+02,.320E+02,.330E+02,.340E+02,.350E+02/
      DATA ET3 /.360E+02,.370E+02,.380E+02,.390E+02,.400E+02,
     &     .410E+02,.420E+02,.430E+02,.440E+02,.450E+02,
     &     .460E+02,.470E+02,.480E+02,.490E+02,.500E+02,
     &     .510E+02,.520E+02,.530E+02,.540E+02,.550E+02,
     &     .560E+02,.570E+02,.580E+02,.590E+02,.600E+02/

      lenpath = lenstr(crstblpath)

      if (mod(MDLQE, 1000)/100 .EQ. 6 
     $     .OR. mod(MDLQE, 1000)/100 .EQ. 7) then
C
C - Spectral function mode
C
         if (iloadeffsfqecrs.ne.0) then
            goto 300
         endif
      else
C
C- Relativistic Fermi Gas mode
C         
         write(*,*) 'fnelspau_effsf: mode (',mdlqe,
     $        ' is not spectral function.'
         stop
      endif
C
C------- load Effective SF cross-section ( only needed )
C******** Effective SF based on existing SF code - C. Wilkinson 2014
C     
C
      if (mod(MDLQE, 1000)/100 .EQ. 7) then
C         DUMSTR = 'qelSfData/temSFTotXsec/'
           write(QEMASTR,'(F4.2)') XMAQE
           DUMSTR = 'qelSfData.'
           DUMSTR = DUMSTR(1:10) // QEMASTR
           DUMSTR = DUMSTR(1:14) // '/temSFTotXsec/'
      else
C         DUMSTR = 'qelSfData/effSFTotXsec/'
           write(QEMASTR,'(F4.2)') XMAQE
           DUMSTR = 'qelSfData.'
           DUMSTR = DUMSTR(1:10) // QEMASTR
           DUMSTR = DUMSTR(1:14) // '/effSFTotXsec/'
      endif
      
      DUMSTR=DUMSTR(1:28) // 'ALL'
      DUMSTR=DUMSTR(1:31) // '_1000000000'
      DUMSTR=DUMSTR(1:42) // '_cc.csv'
      DUMSTR=CRSTBLPATH(1:lenpath)//DUMSTR
      open(LUN, file=DUMSTR, status="old")
      
c     ----- IDUM - neutrino pdg code
c     ----- RDUM - neutrino energy
c     ----- note - values of EE can be overwritten here if the binning of tables changes in the future
c     ----- CRSSFLE1 is re-used here for storing cross-section values for SF implementation
      do i = 1,1090
         if(i.le.210) then ! nue
            read(LUN,*) IDUM,RDUM,
     &           CRSSFLE(1,i),CRSSFLE(2,i),CRSSFLE(3,i),CRSSFLE(4,i),
     &           CRSSFLE(5,i),CRSSFLE(6,i),CRSSFLE(7,i),CRSSFLE(8,i),
     &           CRSSFLE(9,i),CRSSFLE(10,i),CRSSFLE(11,i),CRSSFLE(12,i)
            
            if ( abs(RDUM-EE(I)).gt.0.001 ) then
               write(*,*) "nue ",
     $              'Cross-section table seems to be corrupted.'
               stop
            endif
            
         elseif(i.gt.210 .and. i.le.420) then ! nue-bar
            read(LUN,*) IDUM,RDUM,
     &           CRSSFRE(1,i-210),CRSSFRE(2,i-210),CRSSFRE(3,i-210),
     &           CRSSFRE(4,i-210),CRSSFRE(5,i-210),CRSSFRE(6,i-210),
     &           CRSSFRE(7,i-210),CRSSFRE(8,i-210),CRSSFRE(9,i-210),
     &           CRSSFRE(10,i-210),CRSSFRE(11,i-210),CRSSFRE(12,i-210)
            
            if ( abs(RDUM-EE(I-210)).gt.0.001 ) then
               write(*,*) "nuebar ",
     $              'Cross-section table seems to be corrupted.'
               stop
            endif
            
         elseif(i.gt.420 .and. i.le.630) then ! numu
            read(LUN,*) IDUM,RDUM,
     &           CRSSFLM(1,i-420),CRSSFLM(2,i-420),CRSSFLM(3,i-420),
     &           CRSSFLM(4,i-420),CRSSFLM(5,i-420),CRSSFLM(6,i-420),
     &           CRSSFLM(7,i-420),CRSSFLM(8,i-420),CRSSFLM(9,i-420),
     &           CRSSFLM(10,i-420),CRSSFLM(11,i-420),CRSSFLM(12,i-420)

            if ( abs(RDUM-EE(I-420)).gt.0.001 ) then
               write(*,*) "numu ", 
     $              'Cross-section table seems to be corrupted.'
               stop
            endif
            
         elseif(i.gt.630 .and. i.le.840) then ! numu-bar
            read(LUN,*) IDUM,RDUM,
     &           CRSSFRM(1,i-630),CRSSFRM(2,i-630),CRSSFRM(3,i-630),
     &           CRSSFRM(4,i-630),CRSSFRM(5,i-630),CRSSFRM(6,i-630),
     &           CRSSFRM(7,i-630),CRSSFRM(8,i-630),CRSSFRM(9,i-630),
     &           CRSSFRM(10,i-630),CRSSFRM(11,i-630),CRSSFRM(12,i-630)
            
            if ( abs(RDUM-EE(I-630)).gt.0.001 ) then
               write(*,*) "numubar ", 
     $              'Cross-section table seems to be corrupted.'
               stop
            endif
            
         elseif(i.gt.840 .and. i.le.965) then ! nutau
            read(LUN,*) IDUM,RDUM,
     &           CRSSFLT(1,i-840),CRSSFLT(2,i-840),CRSSFLT(3,i-840),
     &           CRSSFLT(4,i-840),CRSSFLT(5,i-840),CRSSFLT(6,i-840),
     &           CRSSFLT(7,i-840),CRSSFLT(8,i-840),CRSSFLT(9,i-840),
     &           CRSSFLT(10,i-840),CRSSFLT(11,i-840),CRSSFLT(12,i-840)
            
            if ( abs(RDUM-ET(I-840)).gt.0.001 ) then
               write(*,*) "nutau ",
     $              'Cross-section table seems to be corrupted.'
                write (*,*) RDUM, ET(i-840)
               stop
            endif
            
         elseif(i.gt.965 .and. i.le.1090) then ! nutau-bar
            read(LUN,*) IDUM,RDUM,
     &           CRSSFRT(1,i-965),CRSSFRT(2,i-965),CRSSFRT(3,i-965),
     &           CRSSFRT(4,i-965),CRSSFRT(5,i-965),CRSSFRT(6,i-965),
     &           CRSSFRT(7,i-965),CRSSFRT(8,i-965),CRSSFRT(9,i-965),
     &           CRSSFRT(10,i-965),CRSSFRT(11,i-965),CRSSFRT(12,i-965)
            
            if ( abs(RDUM-ET(I-965)).gt.0.001 ) then
               write(*,*) "nutaubar ",
     $              'Cross-section table seems to be corrupted.'
               stop
            endif
            
         endif
      enddo
      
      write (*,*) "CCQE cross section table for the Effective SF loaded"
      iloadeffsfqecrs = 1

 300  continue        
C
C-- Try to return spectral function if available
C
      FNELSPAU_EFFSF=0.
      WeightEnu=0.
      WeightPF=0.
      ii=0
      jj=0
      
      IDX_TARG = 0

      if (NUMATOM .EQ. 2) then
         IDX_TARG = 1
      else if (NUMATOM .EQ. 3) then
         IDX_TARG = 2
      else if (NUMATOM .EQ. 4) then
         IDX_TARG = 3
      else if (NUMATOM .EQ. 12) then
         IDX_TARG = 4
      else if (NUMATOM .EQ. 16) then
         IDX_TARG = 5
      else if (NUMATOM .EQ. 20) then
         IDX_TARG = 6
      else if (NUMATOM .EQ. 27) then
         IDX_TARG = 7
      else if (NUMATOM .EQ. 40) then
         IDX_TARG = 8
      else if (NUMATOM .EQ. 56) then
         IDX_TARG = 9
      else if (NUMATOM .EQ. 63) then
         IDX_TARG = 10
      else if (NUMATOM .EQ. 64) then
         IDX_TARG = 11
      else if (NUMATOM .EQ. 208) then
         IDX_TARG = 12
      endif
      
      if (IDX_TARG.eq.0) THEN
C
C-- Table is not available : use Relativistic Fermi Gas
C
         write(*,*) 'This target (A,Z)=(',NUMBNDP,',',NUMATOM,')'
         goto 100
      endif

      jj        = idx_targ
      if (idx_targ.lt.12) then
         jj1    = idx_targ+1
      else
         jj1    = idx_targ
      endif
      weightpf  = 0.0

C     
C     ++ FOR E-NEUTRINO
C     
*     set weight for interpolation of neutrino energy
      if (iabs(ipar).eq.12) then
         if (e.le.ee(1)) return ! smaller than energy threshold
         ii = 1
         do 11 i=1,209
            ii = i
            if (e.gt.ee(i) .and. e.le.ee(i+1)) then
               goto 101
            endif
 11      continue
101      WeightEnu = (E-EE(ii))/(EE(ii+1)-EE(ii))
         
         if (CRSSFLE(jj, ii) .eq. 0)then
           FNELSPAU_EFFSF = 0
           return
         endif

C     - uwe weightpf == 0 to select primary target
         
         if (ipar.eq.12) then
           if (CRSSFLE(jj, ii) .eq. 0)then
             FNELSPAU_EFFSF = 0
             return
           endif
            if (e.gt.ee(210)) then ! larger than maximum energy
               ii = 210
               CRSSFLOW   = CRSSFLE(jj  ,ii)
               if (idx_targ.ne.9) then
                  CRSSFHIG   = CRSSFLE(jj1,ii)
               else
                  CRSSFHIG   = CRSSFLE(jj,ii)
               endif
               FNELSPAU_EFFSF 
     $              = CRSSFLOW + (CRSSFHIG-CRSSFLOW)*WeightPF
            else
               CRSSFLOW   = CRSSFLE(jj  ,ii) 
     $              + (CRSSFLE(jj  ,ii+1)-CRSSFLE(jj  ,ii))
     $              *WeightEnu
               CRSSFHIG   = CRSSFLE(jj1,ii) 
     $              + (CRSSFLE(jj1,ii+1)-CRSSFLE(jj1,ii))
     $              *WeightEnu
               FNELSPAU_EFFSF = CRSSFLOW 
     $              + (CRSSFHIG-CRSSFLOW)*WeightPF
            endif
         else
           if (CRSSFRE(jj, ii) .eq. 0)then
             FNELSPAU_EFFSF = 0
             return
           endif
            if (e.gt.ee(210)) then ! larger than maximum energy
               ii = 210
               CRSSFLOW   = CRSSFRE(jj  ,ii)
               CRSSFHIG   = CRSSFRE(jj1,ii)
               FNELSPAU_EFFSF = CRSSFLOW 
     $              + (CRSSFHIG-CRSSFLOW)*WeightPF
            else
               CRSSFLOW   = CRSSFRE(jj  ,ii) 
     $              + (CRSSFRE(jj  ,ii+1)-CRSSFRE(jj  ,ii))
     $              *WeightEnu
               CRSSFHIG   = CRSSFRE(jj1,ii) 
     $              + (CRSSFRE(jj1,ii+1)-CRSSFRE(jj1,ii))
     $              *WeightEnu
               FNELSPAU_EFFSF = CRSSFLOW 
     $              + (CRSSFHIG-CRSSFLOW)*WeightPF
            endif
         endif
         return
      endif
C     
C     ++ FOR MU-NEUTRINO
C     
*     set weight for interpolation of neutrino energy
      if (iabs(ipar).eq.14) then
         if (e.le.ee(1)) return ! smaller than energy threshold
         ii = 1
         do 12 i=1,209
            ii = i
            if (e.gt.ee(i) .and. e.le.ee(i+1)) then
               goto 102
            endif

12       continue
102      WeightEnu = (E-EE(ii))/(EE(ii+1)-EE(ii))
      
         
         if (ipar.eq.14) then
c        force a sharp cut off at the turn-on, as the binning isn't fine enough
           if (CRSSFLM(jj, ii) .eq. 0)then
             FNELSPAU_EFFSF = 0
             return
           endif
            if (e.gt.ee(210)) then ! larger than maximum energy
               ii = 210
               CRSSFLOW   = CRSSFLM(jj  ,ii)
               CRSSFHIG   = CRSSFLM(jj1,ii)
               FNELSPAU_EFFSF = CRSSFLOW 
     $              + (CRSSFHIG-CRSSFLOW)*WeightPF
            else
               CRSSFLOW   = CRSSFLM(jj  ,ii) 
     $              + (CRSSFLM(jj  ,ii+1)-CRSSFLM(jj  ,ii))
     $              *WeightEnu
               CRSSFHIG   = CRSSFLM(jj1,ii) 
     $              + (CRSSFLM(jj1,ii+1)-CRSSFLM(jj1,ii))
     $              *WeightEnu
               FNELSPAU_EFFSF = CRSSFLOW 
     $              + (CRSSFHIG-CRSSFLOW)*WeightPF
            endif
         else
c        force a sharp cut off at the turn-on, as the binning isn't fine enough
           if (CRSSFRM(jj, ii) .eq. 0)then
             FNELSPAU_EFFSF = 0
             return
           endif
            if (e.gt.ee(210)) then ! larger than maximum energy
               ii = 210
               CRSSFLOW   = CRSSFRM(jj  ,ii)
               CRSSFHIG   = CRSSFRM(jj1,ii)
               FNELSPAU_EFFSF = CRSSFLOW
     $              + (CRSSFHIG-CRSSFLOW)*WeightPF
            else
               CRSSFLOW   = CRSSFRM(jj  ,ii) 
     $              + (CRSSFRM(jj  ,ii+1)-CRSSFRM(jj  ,ii))
     $              *WeightEnu
               CRSSFHIG   = CRSSFRM(jj1,ii) 
     $              + (CRSSFRM(jj1,ii+1)-CRSSFRM(jj1,ii))
     $              *WeightEnu
               FNELSPAU_EFFSF = CRSSFLOW +
     $              (CRSSFHIG-CRSSFLOW)*WeightPF
            endif
         endif
         return
      endif
C     
C     ++ FOR TAU-NEUTRINO
C     
*     set weight for interpolation of neutrino energy
      if (iabs(ipar).eq.16) then
         if (e.le.et(1)) return ! smaller than energy threshold
         ii = 1
         do 13 i=1,124
            ii = i
            if (e.gt.et(i) .and. e.le.et(i+1)) then
               goto 103
            endif
 13      continue
103      WeightEnu = (E-ET(ii))/(ET(ii+1)-ET(ii))
         
         if (ipar.eq.16) then
           if (CRSSFLT(jj, ii) .eq. 0)then
             FNELSPAU_EFFSF = 0
             return
           endif
            if (e.gt.et(125)) then ! larger than maximum energy
               ii = 125
               CRSSFLOW   = CRSSFLT(jj  ,ii)
               CRSSFHIG   = CRSSFLT(jj1,ii)
               FNELSPAU_EFFSF = CRSSFLOW 
     $              + (CRSSFHIG-CRSSFLOW)*WeightPF
            else
               CRSSFLOW   = CRSSFLT(jj  ,ii) 
     $              + (CRSSFLT(jj  ,ii+1)-CRSSFLT(jj  ,ii))
     $              *WeightEnu
               CRSSFHIG   = CRSSFLT(jj1,ii) 
     $              + (CRSSFLT(jj1,ii+1)-CRSSFLT(jj1,ii))
     $              *WeightEnu
               FNELSPAU_EFFSF = CRSSFLOW 
     $              + (CRSSFHIG-CRSSFLOW)*WeightPF
            endif
         else
           if (CRSSFRT(jj, ii) .eq. 0)then
             FNELSPAU_EFFSF = 0
             return
           endif
            if (e.gt.et(125)) then ! larger than maximum energy
               ii = 125
               CRSSFLOW   = CRSSFRT(jj  ,ii)
               CRSSFHIG   = CRSSFRT(jj1,ii)
               FNELSPAU_EFFSF = CRSSFLOW 
     $              + (CRSSFHIG-CRSSFLOW)*WeightPF
            else
               CRSSFLOW   = CRSSFRT(jj  ,ii) 
     $              + (CRSSFRT(jj  ,ii+1)-CRSSFRT(jj  ,ii))
     $              *WeightEnu
               CRSSFHIG   = CRSSFRT(jj1,ii) 
     $              + (CRSSFRT(jj1,ii+1)-CRSSFRT(jj1,ii))
     $              *WeightEnu
               FNELSPAU_EFFSF = CRSSFLOW 
     $              + (CRSSFHIG-CRSSFLOW)*WeightPF
            endif
         endif
         return
      endif

      write(*,*) 'fnelspau_intpol:Invalid Neutrino particle code?',
     $     ipar
      stop

**************************************
C     call old code if the target is not supported.
C     
C----------Based on Smith-Moniz formula / Code by H.Nakamura
C     
      
 100  continue
      write (*,*) "using RFG"

      FNELSPAU_effsf=fnelspau_intpol(e,ipar)
      return

      END
