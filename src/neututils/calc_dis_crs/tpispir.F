C
      function tpispir(whadron)
C
C     return cross section ratio of
C     more than 1 pi production and 1 pi production.
C
C     INPUT:
C     whadron  invariant mass of hadron system
C     OUTPUT:
C     tpispir  [ >=2pi production] / [ >=1pi production]
C
      Implicit none
      Real tpispir
      Real whadron
      Real fnpkno

      Real meanm, total, totalp, prob, probs, probt
      Integer i

      logical verbose/.false./
      meanm = 0.09 + 1.83*alog(whadron**2.)
      if(verbose) print *, 'mean multiplicity=',meanm,whadron
      total = 0.
      totalp = 0.
      probs = 0.
      probt = 0.

      Do i=1,20
         prob = fnpkno(i,meanm)
         total = total+prob
      Enddo

      Do i=1,20
         prob = fnpkno(i,meanm)
         totalp=totalp+prob/total
         If (i.ge.1) probs=probs+prob
         If (i.ge.2) probt=probt+prob
C         print *, 'P(',i,',',meanm,') =',prob/total
      Enddo
      if(verbose) print *,probt,probs
      tpispir = probt/probs
C      print *, 'total P, tpi/spi:',totalp,tpispir
      Return
      End
