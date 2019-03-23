      SUBROUTINE structg(ipar,it,en,X,Y,f2,xf3,au,ad)
*
      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER (xmp    = 0.93827231D0)
      PARAMETER (xmn    = 0.93956563D0)

      f2 = 0.D0
      xf3 = 0.D0

      if(it.eq.2212) then
        AM    = xmp
      elseif(it.eq.2112) then
        AM    = xmn
      else
        return
      endif

      Q2    = 2*AM*X*Y*EN

*     GRV94
      CALL QGDISG(X,Q2,G,U,D,AU,AD,S,C,B,T)

      AS=S
      AC=C
      AB=B
      AT=T

      if(ipar.gt.0) then
        if(it.eq.2212) then
C nu-proton
          F2 = 2.D0*(D+S+B+AU+AC+AT)
          xF3 = 2.D0*(D+S+B-AU-AC-AT)
        else
C nu-neutron
          F2 = 2.D0*(U+S+B+AD+AC+AT)
          xF3 = 2.D0*(U+S+B-AD-AC-AT)
        endif
      else
        if(it.eq.2212) then
C nubar-proton
	  F2 = 2.D0*(U+C+T+AD+AS+AB)
	  xF3 = 2.D0*(U+C+T-AD-AS-AB)
        else
C nubar-neutron
	  F2 = 2.D0*(D+C+T+AU+AS+AB)
          xF3 = 2.D0*(D+C+T-AU-AS-AB)
        endif
      endif

      RETURN
      END
