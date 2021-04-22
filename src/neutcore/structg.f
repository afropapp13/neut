************************************************************************
*     --------------------------------------
      SUBROUTINE structg(ipar,it,itype,en,X,Y,f2v,f2a,xf3)
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
*     ????.??.?? ; ?? first version
*     2016.07.21 ; C. Bronner
*                  Separate CC and NC structure functions
*                  Add CKM mixing for CC case, with check for charm quarks
*      
*     2017.05.15 ; C. Bronner
*                  Put real expression of NC structure functions
*
*     2019.10.11 ; J. Xia
*                  Moved most of the Bodek-Yang correction from grv98_lo.f
*
*     2020.11.19 ; C. Bronner
*                  Rearrange and use subroutines
************************************************************************



      IMPLICIT REAL*8 (A-H,O-Z)

C      REAL QS
      REAL K_val_u,K_val_d,K_sea_u,K_sea_d
      REAL K_ax_val,K_ax_sea
      REAL X_H, H_Factor
      REAL*8 K_LW
      
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
C     PARAMETER (Vcs    = 0.957D0)

      INCLUDE 'necard.h'   

*     Electroweak couplings
      PARAMETER (sin2w  = 0.2223D0)   !--CODATA2014 value

      REAL*8 gau, gvu, gad, gvd
      gau=0.5D0
      gad=-0.5D0
      gvu=0.5D0-4.D0/3.D0*sin2w
      gvd=-0.5D0+2.D0/3.D0*sin2w
      

      
*     Square of the CKM matrix elements used to compute CC structure functions    
      Vud2=Vud*Vud
      Vus2=Vus*Vus
      
*     Check if a hadron containing a charm quark can be produced
*     If not, set the CKM elements corresponding to charm production to 0
*     Criteria: W large enough to produce a proton and a D0
      Wref=xmp+xmd

*     Set target mass 
      if(it.eq.2212) then
        AM    = xmp
      elseif(it.eq.2112) then
        AM    = xmn
      else
        return
      endif
      
      
      Q2=X*(2.*AM*Y*en)

C     Save the bjoken X for calculation of H factor
      X_H   = X
      
      W2=-Q2+AM**2+2.*AM*Y*en
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



   
*     Get PDF, a priori from GRV98 LO + B-Y corrections, except if specified
*     otherwise in card file
      Q2    = 2*AM*X*Y*EN
      CALL QGDISG(X,Q2,G,U,D,AU,AD,S,C,B,T)

*     Vector variables to have the original PDFs still available for axial terms
      Uvect=U
      Dvect=D
      AUvect=AU
      ADvect=AD
      Svect=S
      Cvect=C
      Bvect=B
      Tvect=T

*     BY correction for vector terms
      if(NEBODEK.ge.1) then
         CALL BYKVEC(Q2,K_val_u,K_val_d,K_sea_u,K_sea_d)

*     BY low W factor (for resonances) - not used by default
         if(NEBODEK.eq.2.and.NEBYLW.eq.1) then
            W_LW_MIN = 1.4D0
            if(W.gt.W_LW_MIN) then
               BJNU2 = DBLE(EN)**2*DBLE(Y)**2 !--- nu
               C_LW = 0.218D0
               K_LW = (BJNU2 + C_LW) / BJNU2
               K_val_u = K_val_u*K_LW
               K_val_d = K_val_d*K_LW
            endif
         endif
         
*     Scaling factor for valence quarks
         Uvect = Uvect*K_val_u
         Dvect = Dvect*K_val_d
*     Scaling factor for sea quarks
         AUvect = AUvect*K_sea_u
         ADvect = ADvect*K_sea_d
         Svect = Svect*K_sea_d
      endif

*     PDFs for axial F2
      if(NEBODEK.eq.2.and.NEBYFORCET1.eq.0) then     !--- specific corrections in new BY
         CALL BYKAX(Q2, K_ax_val,K_ax_sea)
         Uax = U*K_ax_val
         Dax = D*K_ax_val
         AUax = AU*K_ax_sea
         ADax = AD*K_ax_sea
         Sax = S*K_ax_sea
      else     !---- same PDFs for vector and axial
         Uax = Uvect
         Dax = Dvect
         AUax = AUvect
         ADax = ADvect
         Sax = Svect        
      endif
      Cax=C
      Bax=B
      Tax=T
      

*     Add valence and sea components for up and down 
      Uvect = Uvect + AUvect
      Dvect = Dvect + ADvect
      Uax = Uax + AUax
      Dax = Dax + ADax
      
*     For sea quarks, take anti-quark probabilities to be the same as quark one
      ASvect=Svect
      ACvect=Cvect
      ABvect=Bvect
      ATvect=Tvect
      ASax=Sax
      ACax=Cax
      ABax=Bax
      ATax=Tax



*     Compute the structure functions from PDF
*     PDFs are for proton, so for neutron target exchange U<->D
*     Assume no actual contribution from sea charm, so don't put
*     CKM coefficients for C/AC
      if (itype.eq.1) then      ! Charged current structure functions
         if(ipar.gt.0) then     ! neutrino
            if(it.eq.2212) then 
C     nu-proton
               F2V = Dvect*(Vud2+Vcd2)+Svect*(Vus2+Vcs2)+Bvect
     &              +AUvect*(Vud2+Vus2)+ACvect+ATvect
               F2A = Dax*(Vud2+Vcd2)+Sax*(Vus2+Vcs2)+Bax
     &              +AUax*(Vud2+Vus2)+ACax+ATax
               xF3 = 2.D0*(Dvect*(Vud2+Vcd2)+Svect*(Vus2+Vcs2)+Bvect
     &              -AUvect*(Vud2+Vus2)-ACvect-ATvect)
            else
C     nu-neutron
               F2V = Uvect*(Vud2+Vcd2)+Svect*(Vus2+Vcs2)+Bvect
     &              +ADvect*(Vud2+Vus2)+ACvect+ATvect
               F2A = Uax*(Vud2+Vcd2)+Sax*(Vus2+Vcs2)+Bax
     &              +ADax*(Vud2+Vus2)+ACax+ATax
               xF3 = 2.D0*(Uvect*Vud2+Svect*Vus2+Bvect
     &              -ADvect*(Vud2+Vus2)-ACvect-ATvect)
            endif
         else                   ! anti-neutrino
            if(it.eq.2212) then
C     nubar-proton
               F2V = Uvect*(Vud2+Vus2)+Cvect+Tvect+ADvect*(Vud2+Vcd2)
     &              +ASvect*(Vus2+Vcs2)+ABvect
               F2A = Uax*(Vud2+Vus2)+Cax+Tax+ADax*(Vud2+Vcd2)
     &              +ASax*(Vus2+Vcs2)+ABax
               xF3 = 2.D0*(Uvect*(Vud2+Vus2)+Cvect+Tvect
     &              -ADvect*Vud2-ASvect*Vus2-ABvect)
            else
C     nubar-neutron
               F2V = Dvect*(Vud2+Vus2)+Cvect+Tvect+AUvect*(Vud2+Vcd2)
     &              +ASvect*(Vus2+Vcs2)+ABvect
               F2A = Dax*(Vud2+Vus2)+Cax+Tax+AUax*(Vud2+Vcd2)
     &              +ASax*(Vus2+Vcs2)+ABax
               xF3 = 2.D0*(Dvect*(Vud2+Vus2)+Cvect+Tvect
     &               -AUvect*(Vud2+Vcd2)-ASvect*(Vus2+Vcs2)-ABvect)
            endif
         endif
      else                      ! Neutral current case
C     Not clear if should also separate axial and vector for NC...
C     Implement the separation for now, but should reconsider in the future    
         if(it.eq.2212) then 
C     proton target
            F2V = 0.5D0*((gvu*gvu+gau*gau)*(Uvect+Cvect+Tvect+AUvect
     &           +ACvect+ATvect)+ (gvd*gvd+gad*gad)*(Dvect+Svect+Bvect
     &           +ADvect+ASvect+ABvect))
            F2A = 0.5D0*((gvu*gvu+gau*gau)*(Uax+Cax+Tax+AUax+ACax+ATax)
     &           + (gvd*gvd+gad*gad)*(Dax+Sax+Bax+ADax+ASax+ABax))
            xF3 = 2.D0*gvu*gau*(Uvect+Cvect+Tvect-AUvect-ACvect-ATvect)
     &          + 2.D0*gvd*gad*(Dvect+Svect+Bvect-ADvect-ASvect-ABvect)
         else
C     neutron target
            F2V = 0.5D0*((gvu*gvu+gau*gau)*(Dvect+Cvect+Tvect+ADvect           
     &           +ACvect+ATvect)+ (gvd*gvd+gad*gad)*(Uvect+Svect+Bvect
     &           +AUvect+ASvect+ABvect))
            F2A = 0.5D0*((gvu*gvu+gau*gau)*(Dax+Cax+Tax+ADax+ACax+ATax)
     &           + (gvd*gvd+gad*gad)*(Uax+Sax+Bax+AUax+ASax+ABax))
            xF3 = 2.D0*gvu*gau*(Dvect+Cvect+Tvect-ADvect-ACvect-ATvect)
     &           + 2.D0*gvd*gad*(Uvect+Svect+Bvect-AUvect-ASvect-ABvect)
         endif      
      endif

      if(NEBODEK.eq.2) then
         H_Factor = 0.914 + 0.296*X_H - 0.374*X_H*X_H
     &        + 0.165*X_H*X_H*X_H
         xF3 = H_Factor*xF3
      endif

      RETURN
      END
