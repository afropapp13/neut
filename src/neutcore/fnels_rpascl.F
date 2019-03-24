************************************************************************
*     -------------------------
      REAL FUNCTION FNELS_RPASCL(dEE,dQQ,IPAR)
*     -------------------------
*
*     (Purpose)
*       Returns RPA scaling of differential cross-section.
*       
*
*     (Input)
*       dE    : Energy (GeV) in rest frame.
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
C		     values for C12 and O16 is returned; For rest of the elements it returns 1
C		     
C		      Scaling table is available for
C			nue C12 O16, O16, for energies values given by array E1
C			nue C12 for energy values give by array E2
C		      Studies showed that the RPA scaling does not differ for numu and nue.
C		        
C					      
C			
C
************************************************************************
      IMPLICIT NONE
#include "necard.h"      
#include "neutparams.h"      
#include "qerpa.h"
#include "neutmodel.h"
#include "neutfilepath.h"
   

      REAL dQQ
      REAL dEE,dE
      INTEGER IPAR
      
      INTEGER LUN1,LUN2,LUN3,LUN4
      DATA LUN1/80/
      DATA LUN2/81/
      DATA LUN3/82/
      DATA LUN4/83/

      INTEGER I,II,J,JJ,K,KK,L,LL,IDUM,ILOADQERPA
      REAL RDUM
      CHARACTER*1024 DUMSTR1,DUMSTR2,DUMSTR3,DUMSTR4
      SAVE ILOADQERPA
      DATA ILOADQERPA/0/

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
      REAL QQ(300) 
      REAL E1(27)
      
      REAL WeightQQ,WeightE,WeightA
      REAL FNELS_RPASCL_C,FNELS_RPASCL_O
      REAL RPALOW,RPAHIG
      REAL RPA11,RPA12,RPA21,RPA22

      INTEGER*4 IDDUM
      REAL*4    ADUM(2)

      INTEGER*4 NEN
      INTEGER*4 NQQ
      NQQ=300
      NEN=27
      
     
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


C  Energy bins          


         DATA E1/ 0.15, 0.25, 0.35, 0.45, 0.55,
     $         	  0.65, 0.75, 0.85, 0.95, 1.05,
     $         	  1.5,  1.25, 1.35, 1.5,  1.6,  
     $            1.7,  1.8,  1.9,  2.0,  2.2,
     $            2.4,  2.6,  3.0,  4.0,  5.0,
     $         	  7.5,  10.0/



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
C   All cross-section files are read only once when ILOADQERPA = 0
 
       if (iloadqerpa.ne.0) goto 100
C Open four files                
         lenpath = lenstr(crstblpath)
       
    	 DUMSTR1='RPA_corr_E-Q2_numu_c.dat' 
         IDTAG=1
         DUMSTR1=CRSTBLPATH(1:lenpath)//DUMSTR1
         print*, 'RPA cross-section read from file: ', DUMSTR1
         open(LUN1,file=DUMSTR1,form='formatted',status='old')

         DUMSTR2='RPA_corr_E-Q2_numub_c.dat' 
         IDTAG=1
         DUMSTR2=CRSTBLPATH(1:lenpath)//DUMSTR2
         print*, 'RPA cross-section read from file: ', DUMSTR2
         open(LUN2,file=DUMSTR2,form='formatted',status='old')

         DUMSTR3='RPA_corr_E-Q2_numu_o.dat' 
         IDTAG=1
         DUMSTR3=CRSTBLPATH(1:lenpath)//DUMSTR3
         print*, 'RPA cross-section read from file: ', DUMSTR3
         open(LUN3,file=DUMSTR3,form='formatted',status='old')

         DUMSTR4='RPA_corr_E-Q2_numub_o.dat' 
         IDTAG=1
         DUMSTR4=CRSTBLPATH(1:lenpath)//DUMSTR4
         print*, 'RPA cross-section read from file: ', DUMSTR4
         open(LUN4,file=DUMSTR4,form='formatted',status='old')


C Reading table for,
C numu on Carbon
         do i = 1,NQQ

         read(LUN1,*) IDUM,RDUM,RPAL_C_QQ(1,i),RPAL_C_QQ(2,i),
     $  	RPAL_C_QQ(3,i),RPAL_C_QQ(4,i),RPAL_C_QQ(5,i),
     $          RPAL_C_QQ(6,i),RPAL_C_QQ(7,i),RPAL_C_QQ(8,i),
     $          RPAL_C_QQ(9,i),RPAL_C_QQ(10,i),RPAL_C_QQ(11,i),
     $          RPAL_C_QQ(12,i),RPAL_C_QQ(13,i),RPAL_C_QQ(14,i),
     $          RPAL_C_QQ(15,i),RPAL_C_QQ(16,i),RPAL_C_QQ(17,i),
     $          RPAL_C_QQ(18,i),RPAL_C_QQ(19,i),RPAL_C_QQ(20,i),
     $          RPAL_C_QQ(21,i),RPAL_C_QQ(22,i),RPAL_C_QQ(23,i),
     $          RPAL_C_QQ(24,i),RPAL_C_QQ(25,i), RPAL_C_QQ(26,i),
     $  	RPAL_C_QQ(27,i)

         if ( abs(RDUM-QQ(I)).gt.0.001 ) then
         write(*,*)'RPA correction table corrupted 0',RDUM,QQ(I),I
         stop
         endif
         enddo
        
C nue bar on Carbon
         do i = 1,NQQ

         read(LUN2,*) IDUM,RDUM,RPAR_C_QQ(1,i),RPAR_C_QQ(2,i),
     $  	RPAR_C_QQ(3,i),RPAR_C_QQ(4,i),RPAR_C_QQ(5,i),
     $          RPAR_C_QQ(6,i),RPAR_C_QQ(7,i),RPAR_C_QQ(8,i),
     $          RPAR_C_QQ(9,i),RPAR_C_QQ(10,i),RPAR_C_QQ(11,i),
     $          RPAR_C_QQ(12,i),RPAR_C_QQ(13,i),RPAR_C_QQ(14,i),
     $          RPAR_C_QQ(15,i),RPAR_C_QQ(16,i),RPAR_C_QQ(17,i),
     $          RPAR_C_QQ(18,i),RPAR_C_QQ(19,i),RPAR_C_QQ(20,i),
     $          RPAR_C_QQ(21,i),RPAR_C_QQ(22,i),RPAR_C_QQ(23,i),
     $          RPAR_C_QQ(24,i),RPAR_C_QQ(25,i),RPAR_C_QQ(26,i),
     $          RPAR_C_QQ(27,i)

         if ( abs(RDUM-QQ(I)).gt.0.001 ) then
            write(*,*)'RPA correction table corrupted 1.',RDUM ,QQ(I),I
         stop
         endif
C            print*, 'debug:RPAR_C_QQ',i,RPAR_C_QQ(1,i),RPAR_C_QQ(2,i)
     
         enddo


C numu on Oxygen
         do i = 1,NQQ

         read(LUN3,*) IDUM,RDUM,RPAL_O_QQ(1,i),RPAL_O_QQ(2,i),
     $  	RPAL_O_QQ(3,i),RPAL_O_QQ(4,i),RPAL_O_QQ(5,i),
     $          RPAL_O_QQ(6,i),RPAL_O_QQ(7,i),RPAL_O_QQ(8,i),
     $          RPAL_O_QQ(9,i),RPAL_O_QQ(10,i),RPAL_O_QQ(11,i),
     $          RPAL_O_QQ(12,i),RPAL_O_QQ(13,i),RPAL_O_QQ(14,i),
     $          RPAL_O_QQ(15,i),RPAL_O_QQ(16,i) ,RPAL_O_QQ(17,i),
     $          RPAL_O_QQ(18,i),RPAL_O_QQ(19,i),RPAL_O_QQ(20,i),
     $          RPAL_O_QQ(21,i),RPAL_O_QQ(22,i),RPAL_O_QQ(23,i),
     $          RPAL_O_QQ(24,i),RPAL_O_QQ(25,i), RPAL_O_QQ(26,i),
     $          RPAL_O_QQ(27,i)

         if ( abs(RDUM-QQ(I)).gt.0.001 ) then
            write(*,*)'RPA correction table corrupted 2.',RDUM ,QQ(I),I
         stop
         endif
         enddo
C nue bar on Oxygen
         do i = 1,NQQ

         read(LUN4,*) IDUM,RDUM,RPAR_O_QQ(1,i),RPAR_O_QQ(2,i),
     $  	RPAR_O_QQ(3,i),RPAR_O_QQ(4,i),RPAR_O_QQ(5,i),
     $          RPAR_O_QQ(6,i),RPAR_O_QQ(7,i),RPAR_O_QQ(8,i),
     $          RPAR_O_QQ(9,i),RPAR_O_QQ(10,i),RPAR_O_QQ(11,i),
     $          RPAR_O_QQ(12,i),RPAR_O_QQ(13,i),RPAR_O_QQ(14,i),
     $          RPAR_O_QQ(15,i),RPAR_O_QQ(16,i),RPAR_O_QQ(17,i),
     $          RPAR_O_QQ(18,i),RPAR_O_QQ(19,i),RPAR_O_QQ(20,i),
     $          RPAR_O_QQ(21,i),RPAR_O_QQ(22,i),RPAR_O_QQ(23,i),
     $          RPAR_O_QQ(24,i),RPAR_O_QQ(25,i),RPAR_O_QQ(26,i),
     $          RPAR_O_QQ(27,i)


         if ( abs(RDUM-QQ(I)).gt.0.001 ) then
           write(*,*)'RPA correction table corrupted 3.',RDUM ,QQ(I),I
         stop
         endif
C            print*, 'debug:RPAR_C_QQ',i,RPAR_C_QQ(1,i),RPAR_C_QQ(2,i)
     
         enddo

        close(LUN1)
        close(LUN2)
        close(LUN3)
        close(LUN4)

        iloadqerpa = 1
        write(*,*) 'Completed loading RPA diff. table.'
 100    continue
C--------------------------------------------------------------------
        NUMTAR = NUMBNDN+NUMBNDP

C Initialize variables
   
        WeightQQ=0.
        WeightE=0.
	WeightA=0.
        ii=0
        jj=0
        RPALOW=0.
        RPAHIG=0.
	FNELS_RPASCL_C=0.
	FNELS_RPASCL_O=0.
        dE=dEE
C---------------------------------------------------------------------

C       check if the nucleus is isoscaler, if not return 0.      
C        if(NUMBNDN.ne.NUMBNDP) then
C        return
C        endif

        
C       Following condition will change depending on what is decided after meeting nieves
C       when these conditions are removed this function will return interpolated value for given A.
        if(iabs(ipar).eq.12.AND.(NUMTAR.ne.12).and.(NUMTAR.ne.16))then
        write(*,*) 'RPA correction tables not available for this mode.'
        return 
        endif
        if(iabs(ipar).eq.14.AND.(NUMTAR.ne.12).and. (NUMTAR.ne.16))then
        write(*,*) 'RPA correction tables not available for this mode.'
        return 
        endif

C---------------------------------------------------------------------
C      set weight for interpolation over A
        if((NUMTAR.lt.MASSNUM(1)).or.(NUMTAR.gt.MASSNUM(2)))then
        return
        endif

        do 5 l=1,1
        if( (NUMTAR.ge.MASSNUM(l) )
     $        .and.  
     $        (NUMTAR.lt.MASSNUM(l+1)))then
            ll = l
        goto 5
        else if(NUMTAR.eq.MASSNUM(2))then
 	    ll=1
        endif
 5      continue
      
        WeightA
     $     = (REAL(NUMTAR)-REAL(MASSNUM(ll)))
     $     /(REAL(MASSNUM(ll+1))-REAL(MASSNUM(ll)))



C------------------------------------------------------------------------------------
       	        
	if (dE.le.E1(1))then
               jj= 1
        dE=E1(1)
        goto 16
        endif
       
        do 15 j=1,NEN-1
        if( (dE.gt.E1(j) )
     $        .and.  
     $        (dE.le.E1(j+1)))then
            jj = j
        goto 15
        else if(dE.gt.E1(NEN))then
 	    jj=NEN-1
            dE=E1(NEN)
        endif
 15     continue
      
 16    continue 
        WeightE    = (dE-E1(jj))
     $         /(E1(jj+1)-E1(jj))
 

C-----------------------------------------------------------------------------------
        if (dQQ.le.QQ(1))then
               ii= 1
        goto 12
        endif

        do 11 i=1,NQQ-1
        if (dQQ.gt.QQ(i) .and. dQQ.le.QQ(i+1))then
               ii = i
        goto 11
        else if(dQQ.gt.QQ(NQQ)) then
               ii = NQQ-1 
        goto 11
        endif
 11     continue
                     
 12     WeightQQ = (dQQ-QQ(ii))/(QQ(ii+1)-QQ(ii))
	 
C
C ++ FOR Carbon target
C
*       set weight for interpolation of Q**2

        if (iabs(ipar).eq.12 .OR. iabs(ipar).eq.14)then
        if (ipar.eq.12 .OR. ipar.eq.14) then
              RPA11   = RPAL_C_QQ(jj,  ii)
              RPA12   = RPAL_C_QQ(jj,  ii+1)
              RPA21   = RPAL_C_QQ(jj+1,ii)
              RPA22   = RPAL_C_QQ(jj+1,ii+1)
            
        else 
               RPA11   = RPAR_C_QQ(jj  ,ii)
               RPA12   = RPAR_C_QQ(jj  ,ii+1)
               RPA21   = RPAR_C_QQ(jj+1,ii)
               RPA22   = RPAR_C_QQ(jj+1,ii+1)
 
        endif
        endif
C    Bilinear Interpolation

               RPALOW   = RPA11
     $              + (RPA12-RPA11)*WeightQQ
               RPAHIG   = RPA21
     $              + (RPA22-RPA21)*WeightQQ


               FNELS_RPASCL_C = RPALOW
     $        + (RPAHIG-RPALOW)*WeightE

C
C ++ FOR Oxygen target
C
*      set weight for interpolation of Q**2
        if (iabs(ipar).eq.12 .OR. iabs(ipar).eq.14)then
        if (ipar.eq.12 .OR. ipar.eq.14) then
               RPA11   = RPAL_O_QQ(jj,  ii)
               RPA12   = RPAL_O_QQ(jj,  ii+1)
               RPA21   = RPAL_O_QQ(jj+1,ii)
               RPA22   = RPAL_O_QQ(jj+1,ii+1)
        else
               RPA11   = RPAR_O_QQ(jj  ,ii)
               RPA12   = RPAR_O_QQ(jj  ,ii+1)
               RPA21   = RPAR_O_QQ(jj+1,ii)
               RPA22   = RPAR_O_QQ(jj+1 ,ii+1)

        endif
        endif
  
C    Bilinear Interpolation

               RPALOW   = RPA11
     $              + (RPA12-RPA11)*WeightQQ
               RPAHIG   = RPA21
     $              + (RPA22-RPA21)*WeightQQ

C-------------------
               FNELS_RPASCL_O = RPALOW
     $         + (RPAHIG-RPALOW)*WeightE
C-------------------
                FNELS_RPASCL = FNELS_RPASCL_C
     $         + (FNELS_RPASCL_O-FNELS_RPASCL_C)*WeightA
  
                  
        RETURN
        END
