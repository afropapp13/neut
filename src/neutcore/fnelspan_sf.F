************************************************************************
*     -------------------------
      REAL FUNCTION FNELSPAN_SF(E,IPAR,INUC)
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
*       FNELSPAN_SF : CROSS SECTION ( 10**-38cm^2 )
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

      INTEGER LUN
      DATA LUN/80/

      INTEGER I,II,J,JJ,IDUM,ILOADSFELCRS,IQEMA
      REAL RDUM
      CHARACTER*1024 DUMSTR
      CHARACTER*32   QEMASTR

      integer*4      LENSTR
      external       lenstr
      logical HaveSpecFuncForNuclei
      external HaveSpecFuncForNuclei

      integer*4      lenpath

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

      REAL WeightEnu,WeightPF
      REAL CRSLOW,CRSHIG

      INTEGER*4 IDDUM
      REAL*4    PFDUM(9)

      real*4   fnelspan_intpol
      external fnelspan_intpol

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

      character*40 fname_crsdat(4)

      INTEGER*4 IDX_VECFORM,IDX_MA,IDX_NUC,IDX_TARG,NCELMODE

      LENPATH = LENSTR(CRSTBLPATH)

      IF (INUC.eq.2212) then
         IDX_NUC=1
      else if (INUC.eq.2112) then
         IDX_NUC=2
      else
         write(*,*) 'Invalid TARGET ID was specified',INUC
         stop
      endif

      if (mod(MDLQE, 1000)/100 .EQ. 4) then
         if (iloadsfelcrs.ne.0) goto 300
      else
         write(*,*) 'fnelspan_sf: mode (',mdlqe,
     $        ' is not spectral function.'
         stop
      endif

*******  Spectral function implementation *****
* A.Furmanski - 2013
*
      if ( (mod(MDLQE, 1000)/100 .EQ. 4) .AND.
     $     HaveSpecFuncForNuclei() ) then
c ----- load xsec table required
C     DUMSTR = 'qelSfData/totXsec/'
         write(QEMASTR,'(F4.2)') XMAQE
         DUMSTR = 'qelSfData.'
         DUMSTR = DUMSTR(1:10) // QEMASTR
         DUMSTR = DUMSTR(1:14) // '/totXsec/'

        DUMSTR=DUMSTR(1:23) // 'ALL'
        DUMSTR=DUMSTR(1:26) // '_1000000000'
        DUMSTR=DUMSTR(1:37) // '_nc.csv'
        DUMSTR=CRSTBLPATH(1:lenpath)//DUMSTR
        open(LUN, file=DUMSTR, status="old")
c        write (*,*) "reading cross section table for sf"
c ----- RDUM here represents the number of entries in the cross-section table
c ----- note - values of EE can be overwritten here if the binning of tables changes in the future
c ----- CRSSEC is re-used here for storing cross-section values for SF implementation
        do j=1,4
           do i=1,200
              if (j .LE. 2) then
                 read(LUN,*) IDUM,RDUM,
     &             CRSSFELNU(j,1,i),CRSSFELNU(j,2,i),CRSSFELNU(j,3,i),
     &             CRSSFELNU(j,4,i),CRSSFELNU(j,5,i),CRSSFELNU(j,6,i),
     $             CRSSFELNU(j,7,i),CRSSFELNU(j,8,i),CRSSFELNU(j,9,i)
C     print*, 'debug:CRSSFELNU1',i,CRSSFELNU(j,1,i),
C     $             CRSSFELNU(j,2,i),CRSSFELNU(j,3,i)
              else if (j .GT. 2 .AND. j .LE. 4) then
                 read(LUN,*) IDUM,RDUM,
     &         CRSSFELAN(j-2,1,i),CRSSFELAN(j-2,2,i),CRSSFELAN(j-2,3,i),
     &         CRSSFELAN(j-2,4,i),CRSSFELAN(j-2,5,i),CRSSFELAN(j-2,6,i),
     $         CRSSFELAN(j-2,7,i),CRSSFELAN(j-2,8,i),CRSSFELAN(j-2,9,i)
              endif
              if ( abs(RDUM-EE(I)).gt.0.001 ) then
                 write(*,*)
     $                'Cross-section table seems to be corrupted.'
                 stop
              endif
           enddo
        enddo

        write (*,*) "NCEL cross section table for sf loaded"
        iloadsfelcrs = 1
      endif

 300  continue

c ----- find energy wanted and interpolate to get value of total xsec
      FNELSPAN_SF=0.
      WeightEnu=0.
      WeightPF=0.
      ii=0
      jj=0

      IDX_TARG = 0

      if (NUMBNDP .EQ. 6) then
         IDX_TARG = 1
      else if (NUMBNDP .EQ. 8) then
         IDX_TARG = 2
      else if (NUMBNDP .EQ.26)then
         IDX_TARG = 3
      endif

      if (IDX_TARG.ne.0) THEN

         ii = 0
         if (E.le.EE(1)) return ! smaller than energy threshold
         ii = 1
         do 21 i=1,199
            ii = i
            if (E.gt.EE(i) .and. E.le.EE(i+1)) then
               goto 201
            endif
 21      continue
201      WeightEnu = (E-EE(ii))/(EE(ii+1)-EE(ii))


        if (IPAR .GT. 0) then
c ----- force a sharp turn on, to avoid falsely claiming non-zero cross-section
             if (CRSSFELNU(IDX_NUC,IDX_TARG,ii) .eq. 0) then
               FNELSPAN_sf = 0
               return
             endif

             if (e.gt.ee(200)) then ! larger than maximum energy - assume linear extrapolation
                ii = 200
                FNELSPAN_SF = CRSSFELNU(IDX_NUC,IDX_TARG,ii)
             else
                FNELSPAN_SF = CRSSFELNU(IDX_NUC,IDX_TARG,ii)
     $           + (CRSSFELNU(IDX_NUC,IDX_TARG,ii+1)
     $           -CRSSFELNU(IDX_NUC,IDX_TARG,ii))
     $           * WeightEnu
             endif
             return
         else if (IPAR .LT. 0) then
c ----- force a sharp turn on, to avoid falsely claiming non-zero cross-section
             if (CRSSFELAN(IDX_NUC,IDX_TARG,ii) .eq. 0) then
               FNELSPAN_sf = 0
               return
             endif

             if (e.gt.ee(200)) then ! larger than maximum energy - assume linear extrapolation
                ii = 200
                FNELSPAN_SF = CRSSFELAN(IDX_NUC,IDX_TARG,ii)
             else
                FNELSPAN_SF = CRSSFELAN(IDX_NUC,IDX_TARG,ii)
     $           + (CRSSFELAN(IDX_NUC,IDX_TARG,ii+1)
     $           -CRSSFELAN(IDX_NUC,IDX_TARG,ii))
     $           * WeightEnu
             endif
             return
         endif
      else
         write (*,*) "RFG xsec for nc, A =", NUMBNDP
         FNELSPAN_sf=fnelspan_intpol(e,ipar,inuc)
      endif

      END
