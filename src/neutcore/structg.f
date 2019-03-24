************************************************************************
*     --------------------------------------
      SUBROUTINE structg(ipar,it,itype,en,X,Y,f2,xf3)
*
*     --------------------------------------
*
*     (Purpose)
*     ++ Computes structure functions F2 and xF3
*
*     (Input)
*       IPAR   : NEUTRINO SPECIES
*       IT     : Nucleon SPECIES
*       ITYPE  : INTERACTION TYPE
*               =1 CHARGED CURRENT
*               =0 NEUTRAL CURRENT
*       EN     : NEUTRINO ENERGY ( GEV )
*       X      : Q**2/(2Mv) BJORKEN x
*       Y      : v/E        BJORKEN y
*     
*     (Output)
*       F2     : structure function f2
*       xF3    : structure function xf3
*
*     (Creation Date and Author)
*       ????.??.?? ; ?? first version
*       2016.07.21 ; C. Bronner
*                    Separate CC and NC structure functions
*                    Add CKM mixing for CC case, with check for charm quarks
*
************************************************************************



      IMPLICIT REAL*8 (A-H,O-Z)
*     Nucleon masses
      PARAMETER (xmp    = 0.93827231D0)
      PARAMETER (xmn    = 0.93956563D0)

*     D0 meson mass (PDG 2015)
      PARAMETER (xmd = 1.86484D0)

*     CKM parameters. 
*     PDG 2015 values
      PARAMETER (Vud    = 0.97425D0)
      PARAMETER (Vus    = 0.2253D0)
      PARAMETER (Vcd    = 0.225D0)
      PARAMETER (Vcs    = 0.986D0)

*     GENIE 2.10.0 values (used for comparisons)
C      PARAMETER (Vud    = 0.97377D0)
C      PARAMETER (Vus    = 0.2257D0)
C      PARAMETER (Vcd    = 0.230D0)
C      PARAMETER (Vcs    = 0.957D0)


*     Square of the CKM matrix elements used to compute CC structure functions    
      Vud2=Vud*Vud
      Vus2=Vus*Vus
      
*     Check if a hadron containing a charm quark can be produced
*     If not, set the CKM elements corresponding to charm production to 0
*     Criteria: W large enough to produce a proton and a D0
      Wref=xmp+xmd
      
      Q2=X*(2.*xmp*Y*en)
      W2=-Q2+xmp**2+2.*xmp*Y*en
      W=SQRT(W2)
     
      if (W.gt.Wref) then
         Vcd2=Vcd*Vcd
         Vcs2=Vcs*Vcs
      else
        Vcd2=0
        Vcs2=0 
      endif
      
      
      f2 = 0.D0
      xf3 = 0.D0

*     Set target mass 
      if(it.eq.2212) then
        AM    = xmp
      elseif(it.eq.2112) then
        AM    = xmn
      else
        return
      endif

   
*     Get PDF, a priori from GRV98 LO + B-Y corrections, except if specified
*     otherwise in card file
      Q2    = 2*AM*X*Y*EN
      CALL QGDISG(X,Q2,G,U,D,AU,AD,S,C,B,T)

*     For sea quarks, take anti-quark probabilities to be the same as quark one
      AS=S
      AC=C
      AB=B
      AT=T

*     Compute the structure functions from PDF
*     PDFs are for proton, so for neutron target exchange U<->D
*     Assume no actual contribution from sea charm, so don't put
*     CKM coefficients for C/AC
      if (itype.eq.1) then      ! Charged current structure functions
         if(ipar.gt.0) then     ! neutrino
            if(it.eq.2212) then 
C     nu-proton
               F2 = 2.D0*(D*(Vud2+Vcd2)+S*(Vus2+Vcs2)+B+AU*(Vud2+Vus2)
     &              +AC+AT)
               xF3 = 2.D0*(D*(Vud2+Vcd2)+S*(Vus2+Vcs2)+B-AU*(Vud2+Vus2)
     &              -AC-AT)
            else
C     nu-neutron
               F2 = 2.D0*(U*(Vud2+Vcd2)+S*(Vus2+Vcs2)+B+AD*(Vud2+Vus2)
     &              +AC+AT)
               xF3 = 2.D0*(U*Vud2+S*Vus2+B-AD*(Vud2+Vus2)-AC-AT)
            endif
         else                   ! anti-neutrino
            if(it.eq.2212) then
C     nubar-proton
               F2 = 2.D0*(U*(Vud2+Vus2)+C+T+AD*(Vud2+Vcd2)
     &              +AS*(Vus2+Vcs2)+AB)
               xF3 = 2.D0*(U*(Vud2+Vus2)+C+T-AD*Vud2-AS*Vus2-AB)
            else
C     nubar-neutron
               F2 = 2.D0*(D*(Vud2+Vus2)+C+T+AU*(Vud2+Vcd2)
     &              +AS*(Vus2+Vcs2)+AB)
               xF3 = 2.D0*(D*(Vud2+Vus2)+C+T-AU*(Vud2+Vcd2)
     &               -AS*(Vus2+Vcs2)-AB)
            endif
         endif
      else                      ! Neutral current case
         if(ipar.gt.0) then     ! neutrino
            if(it.eq.2212) then 
C     nu-proton
               F2 = 2.D0*(D+S+B+AU+AC+AT)
               xF3 = 2.D0*(D+S+B-AU-AC-AT)
            else
C     nu-neutron
               F2 = 2.D0*(U+S+B+AD+AC+AT)
               xF3 = 2.D0*(U+S+B-AD-AC-AT)
            endif
         else                   ! anti-neutrino
            if(it.eq.2212) then
C     nubar-proton
               F2 = 2.D0*(U+C+T+AD+AS+AB)
               xF3 = 2.D0*(U+C+T-AD-AS-AB)
            else
C     nubar-neutron
               F2 = 2.D0*(D+C+T+AU+AS+AB)
               xF3 = 2.D0*(D+C+T-AU-AS-AB)
            endif
         endif
      endif

      RETURN
      END
