************************************************************************
*     -------------------------
      REAL FUNCTION FNMEC_NIEVES(E,IPAR)
*     -------------------------
*     
*     (Purpose)
*     Returns TOTAL CROSS SECTION OF npnh per nucleon (using Neives model).
*     
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
*     FNMEC_NIEVES : CROSS SECTION  ( 10**-41cm^2 )
*     
*     (Creation Date and Author)
C     
C     
C------------------------------------------------------------------
C     2013.01.25  A.Redij - if Atomic mass number ANUM is given,mec
C     X-section is interpolated for desired value of ANUM,
C     from available data for Carbon and Oxygen, Calcium. 
C     Cross-section tables available for 
C     nue, nuebar, numu, numubar 
C     C12 O16,Ca40 
C     
C     
C     
************************************************************************
      IMPLICIT NONE
#include "necard.h"      
#include "neutparams.h"      
#include "meccrs.h"
#include "neutmodel.h"
#include "neutfilepath.h"

      REAL E
      INTEGER IPAR
      
      INTEGER LUN
      DATA LUN/80/

      INTEGER I,II,J,JJ,IDUM,ILOADMECCRS
      REAL RDUM
      CHARACTER*1024 DUMSTR
      integer*4      LENSTR
      external       lenstr
      integer*4      lenpath

      SAVE ILOADMECCRS
      DATA ILOADMECCRS/0/

      INTEGER*4 IDTAG
      SAVE      IDTAG
      DATA      IDTAG/0/
      
      
      INTEGER*4 NUMTAR

      INTEGER*4 MASSNUM(3)
      DATA MASSNUM/12, 16, 40/

C     Enu for e/mu-nu
      REAL EE1(92) 

      REAL WeightEnu,WeightA
      REAL CRSLOW,CRSHIG
      REAL CRS11,CRS12,CRS21,CRS22

      INTEGER*4 IDDUM
      REAL*4    ADUM(3)

      INTEGER*4 NEN
      NEN=92
C     
C     --  ENERGY OF ELECTRON/MUON NEUTRIN
C     

      DATA  EE1     /0.10  , 0.1025, 0.1050, 0.1075, 0.11  ,
     $     0.1125, 0.1150, 0.1175, 0.12  , 0.1225,
     $     0.1250, 0.1275, 0.13  , 0.1325, 0.1350,
     $     0.1375, 0.14  , 0.1425, 0.145 , 0.1475,
     $     0.15  , 0.1525, 0.155 , 0.1575, 0.16  ,
     $     0.1625, 0.165 , 0.17  , 0.175 , 0.18  ,
     $     0.185 , 0.19  , 0.2   , 0.205 , 0.210 ,
     $     0.22  , 0.23  , 0.24  , 0.25  , 0.26  ,
     $     0.275 , 0.3   , 0.325 , 0.35  , 0.375 ,
     $     0.4   , 0.425 , 0.45  , 0.475 , 0.5   ,
     $     0.525 , 0.55  , 0.6   , 0.65  , 0.7   ,
     $     0.75  , 0.8   , 0.85  , 0.9   , 0.95  ,
     $     1     , 1.05  , 1.1   , 1.15  , 1.2   , 
     $     1.25  , 1.3   , 1.4   , 1.49  , 1.4999999,
     $     1.5000001, 1.502 , 1.5075, 1.513 , 1.515 ,
     $     1.52  , 1.525 , 1.55  , 1.575 , 1.6   ,
     $     1.65  , 1.8   , 1.9   , 2     , 2.5   ,
     $     3     , 4     , 6     , 7.5   , 10    ,
     $     20    , 30/

C-------------------------------------------------------------
      FNMEC_NIEVES=0.

C
C-- at this moment nutau is not supported
C      
      if ((abs(ipar).ne.12).and.(abs(ipar).ne.14)) then
         return
      endif

      lenpath = lenstr(crstblpath)      
C     
C     check is energy is within permissible range.
C     

      if( E .lt. 0.1 .OR. E .gt. 30.0) then
         return
      endif

C--------------------------------------------------------------
C     Cross-section files are read only once
      if (iloadmeccrs.ne.0) goto 100
      DUMSTR='cc_xsec_2p2h_RPA1.dat'
      IDTAG=1

C-----------------------------------------------------------------
C Read data first time the function is called.
      DUMSTR=CRSTBLPATH(1:lenpath)//DUMSTR
      print*, 'MEC cross-section read from file: ', 
     $     DUMSTR(1:LENSTR(DUMSTR))
      open(LUN,file=DUMSTR,form='formatted',status='old')
      read(LUN,*) IDDUM
     $     ,RDUM,ADUM(1),ADUM(2),ADUM(3)
      
      if (IDDUM.ne.IDTAG) then
         write(*,*) 'MEC cross-section table not compatible.'
         stop
      endif
      
      do 200 i = 1, 3
         if (abs(ADUM(i)-MASSNUM(i)).gt.0.001) then
            write(*,*) 'MEC cross-section table is not compatible.'
            stop
         endif
 200  continue
      
*     nue
      do i = 1,NEN

         read(LUN,*) IDUM,RDUM,CRSMLE(1,i),CRSMLE(2,i),CRSMLE(3,i)
         
         if ( abs(RDUM-EE1(I)).gt.0.001 ) then
            write(*,*) 'npnh Cross-section table seems to be corrupted.',RDUM , EE1(I)
            stop
         endif
         
      enddo

*     nue bar
      do i = 1,NEN

         read(LUN,*) IDUM,RDUM,CRSMRE(1,i),CRSMRE(2,i),CRSMRE(3,i)
         
         if ( abs(RDUM-EE1(I)).gt.0.001 ) then
            write(*,*) 'npnh Cross-section table seems to be corrupted.'
            stop
         endif
         
      enddo
*     numu 
      do i = 1,NEN

         read(LUN,*) IDUM,RDUM,CRSMLM(1,i),CRSMLM(2,i),CRSMLM(3,i)
         
         if ( abs(RDUM-EE1(I)).gt.0.001 ) then
            write(*,*) 'npnh Cross-section table seems to be corrupted.'
            stop
         endif
         
      enddo

*     numu bar
      do i = 1,NEN
         read(LUN,*) IDUM,RDUM,CRSMRM(1,i),CRSMRM(2,i),CRSMRM(3,i)
         
         if ( abs(RDUM-EE1(I)).gt.0.001 ) then
            write(*,*) 'npnh Cross-section table seems to be corrupted.'
            stop
         endif
         
      enddo
      

      close(LUN)
      iloadmeccrs = 1
      write(*,*) 'Completed npnh cross-section loading.'
 100  continue
C--------------------------------------------------------------------
C initialise variables.
      WeightEnu=0.
      WeightA=0.
      ii=0
      jj=0
      CRSLOW=0.
      CRSHIG=0.
C---------------------------------------------------------------------
C return if Mass number less than 10 and heavier than 70
      NUMTAR = NUMBNDN+NUMBNDP
C     set weight for interpolation over A
      if((NUMTAR.lt.10).or.(NUMTAR.gt.208)) then
C         write(*,*) 'Cross-section table not available for this mode.'
         return
      endif


      do 10 j=1,2
         if( (NUMTAR.ge.MASSNUM(j))
     $        .and.  
     $        (NUMTAR.lt.MASSNUM(j+1))) then
            jj = j
            goto 10
         else if(NUMTAR.ge.MASSNUM(3))then
            jj=2
         else if(NUMTAR.lt.MASSNUM(1))then
            jj=1
         endif
 10   continue
      
       WeightA  = (REAL(NUMTAR)-REAL(MASSNUM(jj)))
     $     /(REAL(MASSNUM(jj+1))-REAL(MASSNUM(jj)))
       
        

C-----------------------------------------------------------------------------------
      ii = 1
      do 11 i=1,NEN-1
         if (e.gt.ee1(i) .and. e.le.ee1(i+1)) then
            ii = i
            goto 11
         endif
 11   continue
      WeightEnu = (E-EE1(ii))/(EE1(ii+1)-EE1(ii))
      

C     
C     ++ FOR E-NEUTRINO
C     
*     set weight for interpolation of neutrino energy

      if (iabs(ipar).eq.12) then
         if (ipar.eq.12) then
            if (e.gt.ee1(NEN)) then ! larger than maximum energy
               ii = NEN
               CRS11   = CRSMLE(jj  ,ii)
               CRS21   = CRSMLE(jj+1,ii)
            else
               CRS11   = CRSMLE(jj,  ii)
               CRS12   = CRSMLE(jj,  ii+1)
               CRS21   = CRSMLE(jj+1,ii)
               CRS22   = CRSMLE(jj+1,ii+1)
            endif
         else
            if (e.gt.ee1(NEN)) then ! larger than maximum energy
               ii = NEN
               CRS11   = CRSMRE(jj  ,ii)
               CRS21   = CRSMRE(jj+1,ii)
            else
               CRS11   = CRSMRE(jj  ,ii)
               CRS12   = CRSMRE(jj  ,ii+1)
               CRS21   = CRSMRE(jj+1,ii)
               CRS22   = CRSMRE(jj+1,ii+1)
               
            endif
         endif
      endif

C     
C     ++ FOR MU-NEUTRINO
C     
*     set weight for interpolation of neutrino energy
      if (iabs(ipar).eq.14) then
         if (ipar.eq.14) then
            if (e.gt.ee1(NEN)) then ! larger than maximum energy
               ii = NEN
               CRS11   = CRSMLM(jj  ,ii)
               CRS21   = CRSMLM(jj+1,ii)
            else
               CRS11   = CRSMLM(jj,  ii)
               CRS12   = CRSMLM(jj,  ii+1)
               CRS21   = CRSMLM(jj+1,ii)
               CRS22   = CRSMLM(jj+1,ii+1)

            endif
         else
            if (e.gt.ee1(NEN)) then ! larger than maximum energy
               ii = NEN
               CRS11   = CRSMRM(jj  ,ii)
               CRS21   = CRSMRM(jj+1,ii)
            else
               CRS11   = CRSMRM(jj  ,ii)
               CRS12   = CRSMRM(jj  ,ii+1)
               CRS21   = CRSMRM(jj+1,ii)
               CRS22   = CRSMRM(jj+1 ,ii+1)

            endif
         endif
      endif

C     Bilinear Interpolation

      if (e.gt.ee1(NEN)) then ! larger than maximum energy
         CRSLOW   = CRS11
         CRSHIG   = CRS21
      else
         CRSLOW   = CRS11
     $        + (CRS12-CRS11)*WeightEnu
         CRSHIG   = CRS21
     $        + (CRS22-CRS21)*WeightEnu
      endif


      FNMEC_NIEVES = CRSLOW
     $     + (CRSHIG-CRSLOW)*WeightA

C     Interpolated total Cross-section is divided by total number of neucleon.

      FNMEC_NIEVES = FNMEC_NIEVES/(NUMTAR+NUMFREP)
      

      RETURN
      END
