************************************************************************
*     -------------------------
      REAL FUNCTION FNELS_RPASCL(dQQ,IPAR)
*     -------------------------
*
*     (Purpose)
*       Returns RPA scaling of differential cross-section.
*       
*
*     (Input)
*       dQQ   : Q**2 ( GEV**2 )
*       IPAR   : NEUTRINO TYPE
*                  12 : NEU-E
*                 -12 : NEU-E-BAR            
*                  14 : NEU-MU
*                 -14 : NEU-MU-BAR            
*     
*     (Output)
*       FNELS_RPASCL : [dsigma/dQ2(CCQE+RPA)] / [dsigma/dQ2(CCQE)]
*
*     (Creation Date and Author)
*     
*       2013.12.30 A. Redij             
C                         
C------------------------------------------------------------------
C       2013.11.19  A.Redij - if Atomic mass number ANUM is given,
C		     X-section can be interpolated for desired value of ANUM,
C 		     from available data for Carbon and oxygen. 
C 		     At present this feature is disabled and only X-section 
C		     values for C12 and O16 is returned; For rest of the elements it returns zero
C		     
C		      Scaling tables available for 
C			nue C12 O16, 
C			nuebar C12 O16,
C		        numu C12,O16 
C			numubar C12 O16
C		      
C			
C
************************************************************************
      IMPLICIT NONE
#include "necard.h"      
#include "neutparams.h"      
#include "qerpa.h"
#include "neutmodel.h"

   

      REAL dQQ
      INTEGER IPAR
      
      INTEGER LUN
      DATA LUN/80/

      INTEGER I,II,J,JJ,IDUM,ILOADQERPA
      REAL RDUM
      CHARACTER*80 DUMSTR
      SAVE ILOADQERPA
      DATA ILOADQERPA/0/

      INTEGER*4 IDTAG
      SAVE      IDTAG
      DATA      IDTAG/0/
         
      
      INTEGER*4 NUMTAR

      INTEGER*4 MASSNUM(2)
      DATA MASSNUM/12, 16/

C     Enu for e/mu-nu
      REAL QQ(300) 

      REAL WeightQQ,WeightA
      REAL RPALOW,RPAHIG
      REAL RPA11,RPA12,RPA21,RPA22

      INTEGER*4 IDDUM
      REAL*4    ADUM(2)

      INTEGER*4 NEN
      NEN=300
     
C     
C --  Q**2 table
C

      DATA  QQ  /0.002, 0.01, 0.02, 0.03, 0.04, 0.05, 
     $		     0.06, 0.07, 0.08, 0.09, 0.10,
     $               0.11, 0.12, 0.13, 0.14, 0.15, 
     $		0.16, 0.17, 0.18, 0.19, 0.20  ,
     $               0.21, 0.22, 0.23, 0.24, 0.25, 
     $		0.26, 0.27, 0.28, 0.29, 0.30  ,
     $               0.31, 0.32, 0.33, 0.34, 0.35, 
     $		0.36, 0.37, 0.38, 0.39, 0.40  ,
     $               0.41, 0.42, 0.43, 0.44, 0.45, 
     $		0.46, 0.47, 0.48, 0.49, 0.50  ,
     $               0.51, 0.52, 0.53, 0.54, 0.55, 
     $		0.56, 0.57, 0.58, 0.59, 0.60  ,
     $               0.61, 0.62, 0.63, 0.64, 0.65, 
     $		0.66, 0.67, 0.68, 0.69, 0.70  ,
     $               0.71, 0.72, 0.73, 0.74, 0.75, 
     $		0.76, 0.77, 0.78, 0.79, 0.80  ,
     $               0.81, 0.82, 0.83, 0.84, 0.85, 
     $		0.86, 0.87, 0.88, 0.89, 0.90  ,
     $               0.91, 0.92, 0.93, 0.94, 0.95, 
     $		0.96, 0.97, 0.98, 0.99, 1.00  ,
     $               1.01, 1.02, 1.03, 1.04, 1.05, 
     $		1.06, 1.07, 1.08, 1.09, 1.10  ,
     $               1.11, 1.12, 1.13, 1.14, 1.15, 
     $		1.16, 1.17, 1.18, 1.19, 1.20  ,
     $               1.21, 1.22, 1.23, 1.24, 1.25, 
     $		1.26, 1.27, 1.28, 1.29, 1.30  ,
     $               1.31, 1.32, 1.33, 1.34, 1.35, 
     $		1.36, 1.37, 1.38, 1.39, 1.40  ,
     $               1.41, 1.42, 1.43, 1.44, 1.45, 
     $		1.46, 1.47, 1.48, 1.49, 1.50  ,
     $               1.51, 1.52, 1.53, 1.54, 1.55, 
     $		1.56, 1.57, 1.58, 1.59, 1.60  ,
     $               1.61, 1.62, 1.63, 1.64, 1.65, 
     $		1.66, 1.67, 1.68, 1.69, 1.70  ,
     $               1.71, 1.72, 1.73, 1.74, 1.75, 
     $		1.76, 1.77, 1.78, 1.79, 1.80  ,
     $               1.81, 1.82, 1.83, 1.84, 1.85, 
     $		1.86, 1.87, 1.88, 1.89, 1.90  ,
     $               1.91, 1.92, 1.93, 1.94, 1.95, 
     $		1.96, 1.97, 1.98, 1.99, 2.00  ,
     $		     2.01, 2.02, 2.03, 2.04, 2.05, 
     $		2.06, 2.07, 2.08, 2.09, 2.10  ,
     $               2.11, 2.12, 2.13, 2.14, 2.15, 
     $		2.16, 2.17, 2.18, 2.19, 2.20  ,
     $               2.21, 2.22, 2.23, 2.24, 2.25, 
     $		2.26, 2.27, 2.28, 2.29, 2.30  ,
     $               2.31, 2.32, 2.33, 2.34, 2.35, 
     $		2.36, 2.37, 2.38, 2.39, 2.40  ,
     $               2.41, 2.42, 2.43, 2.44, 2.45, 
     $ 		2.46, 2.47, 2.48, 2.49, 2.50  ,
     $               2.51, 2.52, 2.53, 2.54, 2.55, 
     $		2.56, 2.57, 2.58, 2.59, 2.60  ,
     $               2.61, 2.62, 2.63, 2.64, 2.65, 
     $		2.66, 2.67, 2.68, 2.69, 2.70  ,
     $               2.71, 2.72, 2.73, 2.74, 2.75, 
     $		2.76, 2.77, 2.78, 2.79, 2.80  ,
     $               2.81, 2.82, 2.83, 2.84, 2.85, 
     $		2.86, 2.87, 2.88, 2.89, 2.90  ,
     $               2.91, 2.92, 2.93, 2.94, 2.95, 
     $		2.96, 2.97, 2.98, 2.99 /
CC, 3.00  ,
C     $		     3.01, 3.02, 3.03, 3.04, 3.05,
C     $		 3.06, 3.07, 3.08, 3.09, 3.10  ,
C     $               3.11, 3.12, 3.13, 3.14, 3.15, 
C     $		3.16, 3.17, 3.18, 3.19, 3.20  ,
C     $               3.21, 3.22, 3.23, 3.24, 3.25, 
C     $		3.26, 3.27, 3.28, 3.29, 3.30  ,
C     $               3.31, 3.32, 3.33, 3.34, 3.35, 
C     $		3.36, 3.37, 3.38, 3.39, 3.40  ,
C     $               3.41, 3.42, 3.43, 3.44, 3.45, 
C     $		3.46, 3.47, 3.48, 3.49, 3.50  ,
C     $               3.51, 3.52, 3.53, 3.54, 3.55, 
C     $		3.56, 3.57, 3.58, 3.59, 3.60  ,
C     $               3.61, 3.62, 3.63, 3.64, 3.65, 
C     $		3.66, 3.67, 3.68, 3.69, 3.70  ,
C     $               3.71, 3.72, 3.73, 3.74, 3.75, 
C     $		3.76, 3.77, 3.78, 3.79, 3.80  ,
C     $               3.81, 3.82, 3.83, 3.84, 3.85, 
C     $		3.86, 3.87, 3.88, 3.89, 3.90  ,
C     $               3.91, 3.92, 3.93, 3.94, 3.95, 
C     $		3.96, 3.97, 3.98, 3.99/

         

C-------------------------------------------------------------
        
C
C check is energy is within permissible range.
C
        FNELS_RPASCL=1.

        if(mod(MDLQE, 10000)/1000 .NE. 1) then 
            return  
         endif

       if( dQQ.lt. 0.0 .OR. dQQ.gt. 2.99) then
       return 
       endif

C--------------------------------------------------------------
C   Cross-section files are read only once when ILOADQERPA = 0
        if (iloadqerpa.ne.0) goto 100
               
	DUMSTR='RPA_corr_nu_Q2.dat' 
        IDTAG=1

         print*, 'RPA cross-section read from file: ', DUMSTR
         open(LUN,file=DUMSTR,form='formatted',status='old')
         read(LUN,*) IDDUM
     $           ,RDUM,ADUM(1),ADUM(2)

         if (IDDUM.ne.IDTAG) then
          write(*,*) 'RPA correction table is not compatible.'
          stop
         endif

         do 200 i = 1, 2
          if (abs(ADUM(i)-MASSNUM(i)).gt.0.001) then
             write(*,*) 'RPA correction table is not compatible.'
             stop
          endif
 200     continue

*nue
         do i = 1,NEN

           read(LUN,*) IDUM,RDUM,RPALE_QQ(1,i),RPALE_QQ(2,i)
		
            if ( abs(RDUM-QQ(I)).gt.0.001 ) then
               write(*,*) 'RPA correction table seems to be corrupted 0.',RDUM , QQ(I),I
               stop
            endif
     
C            print*, 'debug:RPALM_QQ',i,RPALM_QQ(1,i),RPALM_QQ(2,i)
          enddo

*nue bar
       do i = 1,NEN

           read(LUN,*) IDUM,RDUM,RPARE_QQ(1,i),RPARE_QQ(2,i)
		
            if ( abs(RDUM-QQ(I)).gt.0.001 ) then
        write(*,*) 'RPA correction table seems to be corrupted 0.',RDUM , QQ(I),I
               stop
            endif
C            print*, 'debug:RPARE_QQ',i,RPARE_QQ(1,i),RPARE_QQ(2,i)
     
        enddo
*numu 
         do i = 1,NEN

           read(LUN,*) IDUM,RDUM,RPALM_QQ(1,i),RPALM_QQ(2,i)
		
            if ( abs(RDUM-QQ(I)).gt.0.001 ) then
               write(*,*) 'RPA correction table seems to be corrupted 0 '
               stop
            endif
C            print*, 'debug:RPALM_QQ',i,RPALM_QQ(1,i),RPALM_QQ(2,i)
     
         enddo

*numu bar
       do i = 1,NEN
          read(LUN,*) IDUM,RDUM,RPARM_QQ(1,i),RPARM_QQ(2,i)
		
            if ( abs(RDUM-QQ(I)).gt.0.001 ) then
               write(*,*) 'RPA correction table seems to be corrupted 0'
               stop
            endif
C            print*, 'debug:RPALM_QQ1',i,RPARM_QQ(1,i),RPARM_QQ(2,i)
     
        enddo
	

       close(LUN)
       iloadqerpa = 1
       write(*,*) 'Completed loading RPA diff. table.'
 100   continue
C--------------------------------------------------------------------
       NUMTAR = NUMBNDN+NUMBNDP

        WeightQQ=0.
        WeightA=0.
        ii=0
        jj=0
        RPALOW=0.
        RPAHIG=0.
C---------------------------------------------------------------------

C       check if the nucleus is isoscaler, if not return 0.      
C        if(NUMBNDN.ne.NUMBNDP) then
C        return
C        endif

        
C       Following condition will change depending on what is decided after meeting nieves
C       when these conditions are removed this function will return interpolated value for given A.
        if(iabs(ipar).eq.12 .AND.(NUMTAR.ne.12) .and.(NUMTAR.ne.16) ) then
         write(*,*) 'RPA correction tables not available for this mode.'
         return 
        endif
        if( iabs(ipar).eq.14 .AND. (NUMTAR.ne.12) .and. (NUMTAR.ne.16)) then
         write(*,*) 'RPA correction tables not available for this mode.'
         return 
        endif

C---------------------------------------------------------------------
C      set weight for interpolation over A
        if((NUMTAR.lt.MASSNUM(1)).or.(NUMTAR.gt.MASSNUM(2))) then
         return
        endif

        do 10 j=1,1
         if( (NUMTAR.ge.MASSNUM(j) )
     $        .and.  
     $        (NUMTAR.lt.MASSNUM(j+1))) then
            jj = j
            goto 10
        else if(NUMTAR.eq.MASSNUM(2))then
 	    jj=1
         endif
 10     continue
      
         WeightA
     $     = (NUMTAR-MASSNUM(jj))
     $     /(MASSNUM(jj+1)-MASSNUM(jj))

C-----------------------------------------------------------------------------------
           if (dQQ.le.QQ(1))then
               ii= 1
               goto 12
           endif
           do 11 i=1,NEN-1
           
            if (dQQ.gt.QQ(i) .and. dQQ.le.QQ(i+1)) then
               ii = i
               goto 11
            endif
 11        continue
                     
 12           WeightQQ = (dQQ-QQ(ii))/(QQ(ii+1)-QQ(ii))
	 
C
C ++ FOR E-NEUTRINO
C
*     set weight for interpolation of Q**2

        if (iabs(ipar).eq.12) then
          if (ipar.eq.12) then
              RPA11   = RPALE_QQ(jj,  ii)
              RPA12   = RPALE_QQ(jj,  ii+1)
              RPA21   = RPALE_QQ(jj+1,ii)
              RPA22   = RPALE_QQ(jj+1,ii+1)
          else
               RPA11   = RPARE_QQ(jj  ,ii)
               RPA12   = RPARE_QQ(jj  ,ii+1)
               RPA21   = RPARE_QQ(jj+1,ii)
               RPA22   = RPARE_QQ(jj+1,ii+1)
               
          endif
        endif

C
C ++ FOR MU-NEUTRINO
C
*      set weight for interpolation of Q**2
        if (iabs(ipar).eq.14) then
         if (ipar.eq.14) then
               RPA11   = RPALM_QQ(jj,  ii)
               RPA12   = RPALM_QQ(jj,  ii+1)
               RPA21   = RPALM_QQ(jj+1,ii)
               RPA22   = RPALM_QQ(jj+1,ii+1)

         else
               RPA11   = RPARM_QQ(jj  ,ii)
               RPA12   = RPARM_QQ(jj  ,ii+1)
               RPA21   = RPARM_QQ(jj+1,ii)
               RPA22   = RPARM_QQ(jj+1 ,ii+1)

          endif
        endif

C    Bilinear Interpolation

           RPALOW   = RPA11
     $              + (RPA12-RPA11)*WeightQQ
           RPAHIG   = RPA21
     $              + (RPA22-RPA21)*WeightQQ


         FNELS_RPASCL = RPALOW
     $        + (RPAHIG-RPALOW)*WeightA
 
                  
        RETURN
        END
