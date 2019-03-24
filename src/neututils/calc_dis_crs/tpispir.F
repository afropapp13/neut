C
      function tpispir(whadron,IPAR,ITARG)
C
C     return cross section ratio of
C     more than 1 pi production and 1 pi production.
C
C     INPUT:
C     whadron  invariant mass of hadron system
C     ipar     neutrino type
C     itarg    target nucleon
C     OUTPUT:
C     tpispir  [ >=2pi production] / [ >=1pi production]
C
      Implicit none
      Real tpispir
      Real whadron
      Real fnpkno, fnpkno2

      Real meanm, total, totalp, prob, probs, probt, prob2
      Integer i, ipar, itarg

      real *8 ALPHAKNO
      real *4 CKNO, XMEAN

*     Get parameter values
      CALL NEKNOPAR(whadron,IPAR,ITARG,XMEAN,CKNO,ALPHAKNO)
  
      
      probs = 0.   ! --- probability of 1 pion
      probt = 0.   ! --- probability of more than 1 pion

*     Get probability for the different multiplicities 

      Do i=1,20
C CB  Old NEUT style, KNO for number of pions
C CB  Should be changed for new models
         prob = fnpkno2(i,XMEAN,CKNO,ALPHAKNO)
         If (i.ge.1) probs=probs+prob
         If (i.ge.2) probt=probt+prob   


      Enddo




*     Compute ratio
      tpispir = probt/probs

      Return
      End
