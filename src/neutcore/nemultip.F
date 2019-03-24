**********************************************************************
*     --------------------
      FUNCTION NEMULTIP(W,IPAR,ITARG)
*     --------------------
*     
*     ( purpose )
*        GIVE PION MULTIPLICITY     M.NAKAHATA
*
*     ( input ) 
*        W      : Invariant mass of hadronic system( GEV )
*        IPAR   : NEUTRINO SPECIES
*        ITARG  : TARGET NUCLEON
*
*     ( output )
*        NEMULTIP : PION MULTIPLICITY
*
*     ( creation date and author )
*       198?.??.?? ; First version by M.Nakahata
*       1988.06.27 ; T.KAJITA      PKNO(N,W) ---> PKNO(N,XMEAN)
*       1995.03.03 ; K.Kaneyuki for S.K.
*                    MULTIP -> NEMULTIP
*       2016.07.25 ; C.Bronner
*                    Add possibility to use different multiplicity models,
*                    with choice in card file
*
**********************************************************************
#include "necard.h"
      DIMENSION BUF(10)
      real *4 ALPHAKNO
      real *4 CKNO, XMEAN
      real *4 MAXPROB

*     Get parameter values
      CALL NEKNOPAR(W,IPAR,ITARG,XMEAN,CKNO,ALPHAKNO)
     
*     Fill array with the cumulative probabilities of multiplicities
      BUF(1)=FNPKNO2(2,XMEAN,CKNO,ALPHAKNO)
      DO 100 I=2,9
         BUF(I)=BUF(I-1)+FNPKNO2(I+1,XMEAN,CKNO,ALPHAKNO)
  100 CONTINUE

*     Randomly generate probability using this array
      IF (NEMULT.eq.0) then  ! --- old neut model has an offset of 1, so need to start lower
         BUF(10)=BUF(9)+FNPKNO2(1,XMEAN,CKNO,ALPHAKNO)
         MAXPROB=BUF(10)
      else
         BUF(10)=0
         MAXPROB=BUF(9)
      endif

  110 RAN=RLU(DUM)*MAXPROB
      I=1

  150 IF(BUF(I).GE.RAN)GO TO 200
      I=I+1
      IF(I.EQ.10.and.NEMULT.ne.0) GO TO 110
      if (I.EQ.10.and.BUF(10).GE.RAN) then ! --- old neut model has an offset of 1, so need to start lower
         NEMULTIP=2
         RETURN
      endif
           
      GO TO 150
  200 NEMULTIP=I+1
*     In old NEUT model, probability was for number of pions and not hadrons
      IF (NEMULT.eq.0) then   
         NEMULTIP=I+2      
      endif
      RETURN
      END
