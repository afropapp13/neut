**********************************************************************
*     ------------------------
      FUNCTION FNPKNO2(N,XMEAN,CKNO,ALP)
*     ------------------------
*     
*     ( purpose )
*        PROBABILITY OF PION MULTIPLICITY    M.NAKAHATA
*        (EXECUTED FROM KNO SCALING)
*
*     ( input ) 
*        N     ; hadron MULTIPLICITY
*        XMEAN ; MEAN PION MULTIPLICITY
*        CKNO, ALP : model parameters
*
*     ( output )
*        FNPKNO : PROBABILIRY OF PION MULTIPLICITY
*
*     ( creation date and author )
*       198?.??.?? ; First version by M.Nakahata
*       1988.06.27 ; T.KAJITA      PKNO(N,W) ---> PKNO(N,XMEAN)
*       1995.03.03 ; K.Kaneyuki for S.K.
*                    PKNO -> FNPKNO
*       2016.07.25 ; C. Bronner
*                    Different multiplicity models based on card file
*                    Create FNPKNO2 as FNPKNO used elsewhere
*
**********************************************************************
#include "necard.h"
      DIMENSION PROB(16)
      DATA PROB/1.E-1 , 1.7E-1 , 5.5E-1 , 1.3 , 2. , 2.1 , 1.8 , 1.2 ,
     &     6.5E-1 , 2.6E-1 ,1.E-1 ,3.9E-2 , 1.5E-2 , 6.1E-3 ,
     &     2.5E-3 , 1.E-3/
      real *8 X, CX 
     
   
      FNPKNO2 = 0.
*     Reduced variable
      if ((NEMULT.eq.0).or.(NEMULT.eq.2))then
         Z=FLOAT(N)/XMEAN
         RATIO=1./XMEAN
      else if (NEMULT.eq.1) then
         Z=(FLOAT(N)-ALP)/(XMEAN-ALP)
         RATIO=1./(XMEAN-ALP)
      endif

      if (NEMULT.eq.0) then
*  NEUT old method using interpolation of prob array
         IBIN=IFIX(Z/0.2)+1
         IF(IBIN.GE.16) RETURN        
         XLOG=ALOG10(PROB(IBIN+1))
         XLOG1=ALOG10(PROB(IBIN))
         PLOG=(XLOG-XLOG1)*(Z-0.2*FLOAT(IBIN-1))/0.2+XLOG1
         FNPKNO2=10**PLOG
         RETURN
      else if ((NEMULT.eq.1).or.(NEMULT.eq.2))then
* Other models use Levy function
         X=CKNO*Z+1
         CX=EXP(X*ALOG(CKNO))
         FNPKNO2=2*EXP(-CKNO)*CX/GAMMA(X)
         FNPKNO2=FNPKNO2*RATIO
      endif
      
      RETURN
      END
