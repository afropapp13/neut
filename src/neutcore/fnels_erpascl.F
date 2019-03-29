************************************************************************
*     -------------------------
      REAL FUNCTION FNELS_ERPASCL(dE,IPAR)
*     -------------------------
*
*     (Purpose)
*       Returns RPA scaling for TOTAL CROSS SECTION.
*
*
*     (Input)
*       dE   : E ( GEV )
*       IPAR   : NEUTRINO TYPE
*                  12 : NEU-E
*                 -12 : NEU-E-BAR
*                  14 : NEU-MU
*                 -14 : NEU-MU-BAR
*
*     (Output)
*       FNELS_ERPASCL : [sigma(CCQE+RPA)] / [sigma(CCQE)]
*
*     (Creation Date and Author)
*       2013-Dec Asmita Redij
C
C
************************************************************************
      IMPLICIT NONE
#include "necard.h"
#include "neutparams.h"
#include "qeerpa.h"
#include "neutmodel.h"
#include "neutfilepath.h"


      REAL dE
      INTEGER IPAR

      INTEGER LUN
      DATA LUN/80/

      INTEGER I,II,J,JJ,IDUM,ILOADQERPAE
      REAL RDUM
      CHARACTER*1024 DUMSTR
      SAVE ILOADQERPAE
      DATA ILOADQERPAE/0/

      integer*4      LENSTR
      external       lenstr
      integer*4      lenpath

      INTEGER*4 IDTAG
      SAVE      IDTAG
      DATA      IDTAG/0/


      INTEGER*4 NUMTAR

      INTEGER*4 MASSNUM(2)
      DATA MASSNUM/12, 16/

C     Enu for e/mu-nu
      REAL E(210)

      REAL WeightE,WeightA
      REAL RPALOW,RPAHIG
      REAL RPA11,RPA12,RPA21,RPA22

      INTEGER*4 IDDUM
      REAL*4    ADUM(2)

      INTEGER*4 NEN

C
C --  E table
C

      DATA  E  /0.025,0.075,0.125,0.175,0.225,
     $          0.275,0.325,0.375,0.425,0.475,
     $          0.525,0.575,0.625,0.675,0.725,
     $          0.775,0.825,0.875,0.925,0.975,
     $          1.025,1.075,1.125,1.175,1.225,
     $          1.275,1.325,1.375,1.425,1.475,
     $          1.525,1.575,1.625,1.675,1.725,
     $          1.775,1.825,1.875,1.925,1.975,
     $          2.025,2.075,2.125,2.175,2.225,
     $          2.275,2.325,2.375,2.425,2.475,
     $          2.525,2.575,2.625,2.675,2.725,
     $          2.775,2.825,2.875,2.925,2.975,
     $          3.025,3.075,3.125,3.175,3.225,
     $          3.275,3.325,3.375,3.425,3.475,
     $          3.525,3.575,3.625,3.675,3.725,
     $          3.775,3.825,3.875,3.925,3.975,
     $          4.025,4.075,4.125,4.175,4.225,
     $          4.275,4.325,4.375,4.425,4.475,
     $          4.525,4.575,4.625,4.675,4.725,
     $          4.775,4.825,4.875,4.925,4.975,
     $          5.025,5.075,5.125,5.175,5.225,
     $          5.275,5.325,5.375,5.425,5.475,
     $          5.525,5.575,5.625,5.675,5.725,
     $          5.775,5.825,5.875,5.925,5.975,
     $          6.025,6.075,6.125,6.175,6.225,
     $          6.275,6.325,6.375,6.425,6.475,
     $          6.525,6.575,6.625,6.675,6.725,
     $          6.775,6.825,6.875,6.925,6.975,
     $          7.025,7.075,7.125,7.175,7.225,
     $          7.275,7.325,7.375,7.425,7.475,
     $          7.525,7.575,7.625,7.675,7.725,
     $          7.775,7.825,7.875,7.925,7.975,
     $          8.025,8.075,8.125,8.175,8.225,
     $          8.275,8.325,8.375,8.425,8.475,
     $          8.525,8.575,8.625,8.675,8.725,
     $          8.775,8.825,8.875,8.925,8.975,
     $          9.025,9.075,9.125,9.175,9.225,
     $          9.275,9.325,9.375,9.425,9.475,
     $          9.525,9.575,9.625,9.675,9.725,
     $          9.775,9.825,9.875,9.925,9.975,
     $          11.00,12.00,13.00,14.00,15.00,
     $          16.00,17.00,18.00,19.00,20.00/

      NEN=210

      lenpath = lenstr(crstblpath)

C-------------------------------------------------------------
C
C check is energy is within permissible range.
C

        FNELS_ERPASCL=1.

        if(mod(MDLQE, 10000)/1000 .NE. 1) then
            return
         endif

        if( dE.lt. 0.0 .OR. dE.gt. 30.0) then
        return
        endif

C--------------------------------------------------------------
       if (iloadqerpae.ne.0) goto 100
       DUMSTR = 'RPA_corr_nu_E.dat'
       DUMSTR=CRSTBLPATH(1:lenpath)//DUMSTR
        open(LUN,file=DUMSTR,form='formatted',status='old')
        do i = 1,NEN
            read(LUN,*)IDUM,RDUM,RPALE_E(1,i),RPALE_E(2,i)
            if ( abs(RDUM-E(I)).gt.0.001 ) then
               write(*,*) 'RPA table seems to be corrupted.'
               stop
            endif
        enddo
        do i = 1,NEN
              read(LUN,*)IDUM,RDUM,RPARE_E(1,i),RPARE_E(2,i)

            if ( abs(RDUM-E(I)).gt.0.001 ) then
               write(*,*) 'RPA table seems to be corrupted.'
               stop
            endif
         enddo
         do i = 1,NEN
            read(LUN,*)IDUM,RDUM,RPALM_E(1,i),RPALM_E(2,i)

            if ( abs(RDUM-E(I)).gt.0.001 ) then
               write(*,*) 'RPA table seems to be corrupted.'
               stop
            endif
         enddo
         do i = 1,NEN
            read(LUN,*)IDUM,RDUM,RPARM_E(1,i),RPARM_E(2,i)

            if ( abs(RDUM-E(I)).gt.0.001 ) then
               write(*,*) 'RPA table seems to be corrupted.'
               stop
            endif
         enddo

       close(LUN)
       iloadqerpae= 1
       write(*,*) 'Completed loading RPA correction vs Energy table.'
 100   continue
C--------------------------------------------------------------------
       NUMTAR = NUMBNDN+NUMBNDP

        WeightE=0.
        WeightA=0.
        ii=0
        jj=1
        RPALOW=0.
        RPAHIG=0.
C---------------------------------------------------------------------
C       this condition will be removed on getting Ca table, once A independence is studied
        if(iabs(ipar).eq.12 .AND.(NUMTAR.ne.12) .and.(NUMTAR.ne.16) ) then
         write(*,*) 'RPA correction tables not available for this mode.'
         return
        endif
        if(iabs(ipar).eq.14 .AND.(NUMTAR.ne.12) .and. (NUMTAR.ne.16)) then
         write(*,*) 'RPA correction tables not available for this mode.'
         return
        endif

           if(dE.lt.E(1)) return

           do 11 i=1,NEN-1

           if (dE.ge.E(i) .and. dE.lt.E(i+1))then
               ii = i
               goto 11
            endif

 11        continue


          if (dE.ge.E(NEN)) then
           ii = NEN-1
          endif

 12       WeightE = (dE-E(ii))/(E(ii+1)-E(ii))

          WeightA  = (REAL(NUMTAR)-REAL(MASSNUM(jj)))
     $     /(REAL(MASSNUM(jj+1))-REAL(MASSNUM(jj)))

C
C ++ FOR E-NEUTRINO
C

          if (iabs(ipar).eq.12) then
          if (ipar.eq.12) then

               RPA11   = RPALE_E(jj,  ii)
               RPA12   = RPALE_E(jj,  ii+1)
               RPA21   = RPALE_E(jj+1,ii)
               RPA22   = RPALE_E(jj+1,ii+1)
          else
               RPA11   = RPARE_E(jj  ,ii)
               RPA12   = RPARE_E(jj  ,ii+1)
               RPA21   = RPARE_E(jj+1,ii)
              RPA22   = RPARE_E(jj+1,ii+1)
          endif
          endif

C
C ++ FOR MU-NEUTRINO
C

         if (iabs(ipar).eq.14) then
         if (ipar.eq.14) then
               RPA11   = RPALM_E(jj,  ii)
               RPA12   = RPALM_E(jj,  ii+1)
               RPA21   = RPALM_E(jj+1,ii)
               RPA22   = RPALM_E(jj+1,ii+1)
         else
               RPA11   = RPARM_E(jj  ,ii)
               RPA12   = RPARM_E(jj  ,ii+1)
               RPA21   = RPARM_E(jj+1,ii)
               RPA22   = RPARM_E(jj+1,ii+1)
         endif
         endif

C    Bilinear Interpolation

         RPALOW   = RPA11
     $        + (RPA12-RPA11)*WeightE
         RPAHIG   = RPA21
     $        + (RPA22-RPA21)*WeightE

         FNELS_ERPASCL = RPALOW
     $     + (RPAHIG-RPALOW)*WeightA

         RETURN
         END
