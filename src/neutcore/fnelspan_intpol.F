************************************************************************
*     -------------------------
      REAL FUNCTION FNELSPAN_INTPOL(E,IPAR,INUC)
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
*       FNELSPAN_INTPOL : CROSS SECTION ( 10**-38cm^2 )
*
*     (Creation Date and Author)
C       2011.03.03 ; Y.Hayato  - For neutral current
************************************************************************
      IMPLICIT NONE

#include "neutparams.h"      
#include "elcrs.h"
#include "neutmodel.h"
#include "necard.h"
#include "neutfilepath.h"

      REAL E
      INTEGER IPAR,INUC

      INTEGER LUN(4)
      DATA LUN/80,81,82,83/

      INTEGER I,II,J,JJ,IDUM,ILOADELCRS,ILOADSFELCRS,IQEMA
      REAL RDUM
      CHARACTER*80 DUMSTR
      integer*4      LENSTR
      external       lenstr
      integer*4      lenpath
      SAVE ILOADELCRS
      DATA ILOADELCRS/0/
      SAVE ILOADSFELCRS
      DATA ILOADSFELCRS/0/

      INTEGER*4 IDTAG
      SAVE      IDTAG
      DATA      IDTAG/0/

      REAL FERMISURF(9)
      DATA FERMISURF/0.205, 0.215, 0.225, 0.235, 0.245,
     &               0.255, 0.265, 0.275, 0.285/


C     Enu
      REAL EE(200) ,EE1(50) ,EE2(50) ,EE3(50) ,EE4(50)

C     Cross-section
      REAL CRSSEC(200)

      REAL WeightEnu,WeightPF
      REAL CRSLOW,CRSHIG

      INTEGER*4 IDDUM
      REAL*4    PFDUM(9)

      EQUIVALENCE (EE(1),EE1(1))
      EQUIVALENCE (EE(51),EE2(1))
      EQUIVALENCE (EE(101),EE3(1))
      EQUIVALENCE (EE(151),EE4(1))

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

      character*15 fname_prefix_nu(2)
      character*16 fname_prefix_nub(2)
      character*7  fname_vecform(3)
      character*9  fname_ma(4)
      
      data fname_prefix_nu   /"ncel_nu_p_xsec_", "ncel_nu_n_xsec_"/
      data fname_prefix_nub  /"ncel_nub_p_xsec_","ncel_nub_n_xsec_"/
      data fname_vecform  /"dipole_","bbba05_","bbba07_"/
      data fname_ma       /"ma1.0.dat","ma1.1.dat",
     $                     "ma1.2.dat","ma1.3.dat"/

      character*40 fname_crsdat(4)

      INTEGER*4 IDX_VECFORM,IDX_MA,IDX_NUC,NCELMODE

      lenpath = lenstr(CRSTBLPATH)

       if (MDLQE/100 .EQ. 4 .AND. (numbndp .eq. 6 .or. 
     &    numbndp .eq. 8 .or. numbndp .eq. 26)) then
c      if (MDLQE/100 .EQ. 4) then
        if (iloadsfelcrs.ne.0) goto 300
      else
        if (iloadelcrs.ne.0) goto 100
      endif

*******  Spectral function implementation *****
* A.Furmanski - 2013
*
      if (MDLQE/100 .EQ. 4 .AND. (numbndp .eq. 6 .or. 
     &    numbndp .eq. 8 .or. numbndp .eq. 26)) then
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
        
        DUMSTR=DUMSTR(1:32) // '_nc.csv'
        DUMSTR=CRSTBLPATH(1:lenpath)//DUMSTR
        open(LUN(1), file=DUMSTR, status="old")
c        write (*,*) "reading cross section table for sf"
c ----- RDUM here represents the number of entries in the cross-section table
c ----- note - values of EE can be overwritten here if the binning of tables changes in the future
c ----- CRSSEC is re-used here for storing cross-section values for SF implementation
        read(LUN(1),*) RDUM
        do i = 1, RDUM
          read(LUN(1),*) EE(i), CRSSEC(i)
        enddo
        close(LUN(1))
        iloadsfelcrs = 1
c        write (*,*) "cross section table for sf loaded"

300     continue

c ----- find energy wanted and interpolate to get value of total xsec
        ii = 0
        if (E.le.EE(1)) return ! smaller than energy threshold
        ii = 1
        do 21 i=1,199
          if (E.gt.EE(i) .and. E.le.EE(i+1)) then
            ii = i
            goto 21
          endif
21      continue
        WeightEnu = (E-EE(ii))/(EE(ii+1)-EE(ii))       

        if (e.gt.ee(200)) then ! larger than maximum energy - assume linear extrapolation
          ii = 200
          FNELSPAN_INTPOL = CRSSEC(ii-1)
     $              + (CRSSEC(ii)-CRSSEC(ii-1))*WeightEnu/10
        else
          FNELSPAN_INTPOL = CRSSEC(ii) 
     $              + (CRSSEC(ii+1)-CRSSEC(ii))*WeightEnu/10
        endif
        return
      endif

********************************************************************

C ------ old code - used on nuclei where new SF implementation is not possible

C
C----------Based on Spectral Func. / Code by T.Mori ( Okayama Univ. )
C

      if (iloadelcrs.ne.0) goto 100
      NCELMODE = mod(MDLQE,100)/10
      if (MDLQE/100 .eq. 4) then
      NCELMODE = mod(MDLQE,400)
      endif


      print*, 'Loading Cross section table for NCEL'
      if (NCELMODE.eq.1) then
         IDX_VECFORM=1
         write (*,*) ' --- NC Vector form factor : Dipole'
      else if (NCELMODE.eq.2) then
         IDX_VECFORM=2
         write (*,*) ' --- NC Vector form factor : bbba05'
      else if (NCELMODE.eq.3) then
         IDX_VECFORM=3
         write (*,*) ' --- NC Vector form factor : bbba07'
      else
         write(*,*) 'Invalid MDLQE was specified ',MDLQE
         stop
      endif

      IQEMA = INT(XMAQE*100+0.5)
c      IDTAG = IQEMA
      IDTAG = INT(XMAQE*100-0.5)

      if (IQEMA.eq.101) then
         IDX_MA=1
         write (*,*) ' --- NC MA : 1.01'
      else if (IQEMA.eq.111) then
         IDX_MA=2
         write (*,*) ' --- NC MA : 1.11'
      else if (IQEMA.eq.121) then
         IDX_MA=3
         write (*,*) ' --- NC MA : 1.21'
      else if (IQEMA.eq.131) then
         IDX_MA=4
         write (*,*) ' --- NC MA : 1.31'
      else
         write(*,*) 'Invalid MA for QE was specified ',XMAQE
         stop
      endif

      DO 1000 I=1,4
         if (I.le.2) THEN
            fname_crsdat(I)=
     $           CRSTBLPATH(1:lenpath)//
     $           fname_prefix_nu(I)//
     $           fname_vecform(IDX_VECFORM)//
     $           fname_ma(IDX_MA)
         else
            fname_crsdat(I)=
     $           CRSTBLPATH(1:lenpath)//
     $           fname_prefix_nub(I-2)//
     $           fname_vecform(IDX_VECFORM)//
     $           fname_ma(IDX_MA)
         endif            
         open(LUN(I),file=fname_crsdat(I),form='formatted',
     $        status='old')
C-- skip 2 lines
         read(LUN(I),*) DUMSTR
         read(LUN(I),*) DUMSTR

         read(LUN(I),*) IDDUM,RDUM,
     $        PFDUM(1),PFDUM(2),PFDUM(3),PFDUM(4),PFDUM(5),
     $        PFDUM(6),PFDUM(7),PFDUM(8),PFDUM(9)

         if (IDDUM.ne.IDTAG) then
            write(*,*) 'This NCEL file is not compatible.(Model)',
     $           fname_crsdat(I)
          write (*,*)"file: ",IDDUM,", code: ",IDTAG,", xmaqe= ",XMAQE
            stop
         endif

         do 200 j = 1, 9
            if (abs(PFDUM(j)-FERMISURF(j)).gt.0.001) then
               write(*,*) 'This NCEL file is not compatible.(PF)',
     $              fname_crsdat(I)
               stop
            endif
 200     continue

 1000 continue

      do i = 1,200
         do J=1,2
            read(LUN(J),*) IDUM,RDUM,
     &           CRSRFGELNU(J,1,i),CRSRFGELNU(J,2,i),CRSRFGELNU(J,3,i),
     &           CRSRFGELNU(J,4,i),CRSRFGELNU(J,5,i),CRSRFGELNU(J,6,i),
     $           CRSRFGELNU(J,7,i),CRSRFGELNU(J,8,i),CRSRFGELNU(J,9,i)
            if ( abs(RDUM-EE(I)).gt.0.001 ) then
               write(*,*) 'Cross-section table seems to be corrupted.'
               stop
            endif

            read(LUN(2+J),*) IDUM,RDUM,
     &           CRSRFGELAN(J,1,i),CRSRFGELAN(J,2,i),CRSRFGELAN(J,3,i),
     &           CRSRFGELAN(J,4,i),CRSRFGELAN(J,5,i),CRSRFGELAN(J,6,i),
     $           CRSRFGELAN(J,7,i),CRSRFGELAN(J,8,i),CRSRFGELAN(J,9,i)
            
            if ( abs(RDUM-EE(I)).gt.0.001 ) then
               write(*,*) 'Cross-section table seems to be corrupted.'
               stop
            endif
         enddo
      enddo

      DO 1010 I=1,4
         close(LUN(I))
 1010 continue

      iloadelcrs = 1
      write(*,*) 'Completed NCEL cross-section loading.'

 100  continue

      IF (INUC.eq.2212) then
         IDX_NUC=1
      else if (INUC.eq.2112) then
         IDX_NUC=2
      else
         write(*,*) 'Invalid TARGET ID was specified',INUC
         stop
      endif

      FNELSPAN_INTPOL=0.
      WeightEnu=0.
      WeightPF=0.
      ii=0
      jj=0
      CRSLOW=0.
      CRSHIG=0.

C     set weight for interpolation of fermi surface momentum
      if( (PFSURF.le.FERMISURF(1)) .or. 
     $     (PFSURF.gt.FERMISURF(9))) then
         write(*,*) 'fnelspan_intpol : Invalid Fermi surface mom.',
     $       PFSURF, FERMISURF(1), FERMISURF(9)
         stop
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
C ++ FOR NEUTRINO
C
*     set weight for interpolation of neutrino energy
      if (e.le.ee(1)) return ! smaller than energy threshold
      ii = 1
      do 11 i=1,199
         if (e.gt.ee(i) .and. e.le.ee(i+1)) then
            ii = i
            goto 11
         endif
 11   continue
      WeightEnu = (E-EE(ii))/(EE(ii+1)-EE(ii))
      
      if (ipar.gt.0) then
         if (e.ge.ee(200)) then ! larger than maximum energy
            ii = 200
            CRSLOW   = CRSRFGELNU(IDX_NUC,jj  ,ii)
            CRSHIG   = CRSRFGELNU(IDX_NUC,jj+1,ii)
            FNELSPAN_INTPOL = CRSLOW + (CRSHIG-CRSLOW)*WeightPF
         else
            CRSLOW   = CRSRFGELNU(IDX_NUC,jj  ,ii) 
     $           + ( CRSRFGELNU(IDX_NUC,jj  ,ii+1)
     $              -CRSRFGELNU(IDX_NUC,jj  ,ii))*WeightEnu
            CRSHIG   = CRSRFGELNU(IDX_NUC,jj+1,ii) 
     $           + ( CRSRFGELNU(IDX_NUC,jj+1,ii+1)
     $              -CRSRFGELNU(IDX_NUC,jj+1,ii))*WeightEnu
            FNELSPAN_INTPOL = CRSLOW + (CRSHIG-CRSLOW)*WeightPF
         endif
      else
         if (e.ge.ee(200)) then ! larger than maximum energy
            ii = 200
            CRSLOW   = CRSRFGELAN(IDX_NUC,jj  ,ii)
            CRSHIG   = CRSRFGELAN(IDX_NUC,jj+1,ii)
            FNELSPAN_INTPOL = CRSLOW + (CRSHIG-CRSLOW)*WeightPF
         else
            CRSLOW   = CRSRFGELAN(IDX_NUC,jj  ,ii) 
     $              + ( CRSRFGELAN(IDX_NUC,jj  ,ii+1)
     $                 -CRSRFGELAN(IDX_NUC,jj  ,ii))*WeightEnu
            CRSHIG   = CRSRFGELAN(IDX_NUC,jj+1,ii) 
     $              + ( CRSRFGELAN(IDX_NUC,jj+1,ii+1)
     $                 -CRSRFGELAN(IDX_NUC,jj+1,ii))*WeightEnu
            FNELSPAN_INTPOL = CRSLOW + (CRSHIG-CRSLOW)*WeightPF
         endif
      endif

      RETURN
      END
