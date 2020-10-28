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
*       F2v    : structure function f2 vector part
*       F2a    : structure function f2 axial part     
*       xF3    : structure function xf3
*
*     (Creation Date and Author)
*       ????.??.?? ; ?? first version
*       2016.07.21 ; C. Bronner
*                    Separate CC and NC structure functions
*                    Add CKM mixing for CC case, with check for charm quarks
*       2019.10.11 ; J. Xia
*                    Moved most of the Bodek-Yang correction to from grv98_lo.f
*                    Kept the xi scaling and d/u assymetry in the old place  
*                    Separated vector and axial part F2
************************************************************************
      IMPLICIT REAL*8 (A-H,K,L,O-Z)
C      IMPLICIT NONE
C#include "necard.h"
*     Nucleon masses
C      REAL*8 xmp, xmn, xmd, xmc, am, Q2
C      REAL*8 AM
C      REAL*8 Vud, Vus, Vcd, Vcs
C      REAL*8 f2, xf3
C      REAL BJNU2, X_H
C      REAL UTOT, DTOT,
C      REAL AD_F3, AU_F3, U_F3, D_F3, S_F3, C_F3, AS_F3, AC_F3
C      REAL W_LW_MIN, Wref, W2, W
            
      PARAMETER (xmp    = 0.93827231D0)
      PARAMETER (xmn    = 0.93956563D0)

*     D0 meson mass (PDG 2015)
      PARAMETER (xmd = 1.86484D0)
      
*     Charm quark mass (PDG 2019)      
      PARAMETER (xmc = 1.2705D0)

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

      INCLUDE 'necard.h'
      
*     Square of the CKM matrix elements used to compute CC structure functions    
      Vud2=Vud*Vud
      Vus2=Vus*Vus
      
*     Check if a hadron containing a charm quark can be produced
*     If not, set the CKM elements corresponding to charm production to 0
*     Criteria: W large enough to produce a proton and a D0
      W_LW_MIN = 1.4D0
      Wref     = xmp+xmd
*     Switch to change target mass
      AM       = xmp
      if(it.eq.2112) then
        AM     = xmn
      endif
   
*     Get PDF, a priori from GRV98 LO + B-Y corrections, except if specified
*     otherwise in card file
      Q2    = 2.D0*AM*DBLE(X)*DBLE(Y)*DBLE(EN)
*     Save the bjoken X for calculation of H factor
      X_H   = DBLE(X)

      W2=-Q2+AM**2+2.D0*AM*DBLE(Y)*DBLE(en)
      W=SQRT(W2)
     
      if (W.gt.Wref) then
         Vcd2=Vcd*Vcd
         Vcs2=Vcs*Vcs
      else
         Vcd2=0
         Vcs2=0
      endif
      
C     Initialize parameters
      f2v         = 0.D0
      f2a         = 0.D0
      xf3         = 0.D0
      K_LW        = 1.D0
      K_vec_sea_c = 1.D0
      H_factor    = 1.D0
      K_vec_val_u = 1.D0
      K_vec_val_d = 1.D0
      K_vec_sea_u = 1.D0
      K_vec_sea_d = 1.D0
      K_vec_sea_c = 1.D0
      K_axi_val   = 1.D0
      K_axi_sea   = 1.D0
      
      C_vec_val1_d = 0.D0
      C_vec_val2_d = 0.D0
      C_vec_val1_u = 0.D0
      C_vec_val2_u = 0.D0
      C_vec_sea_d  = 0.D0
      C_vec_sea_u  = 0.D0
      C_axi_sea    = 0.D0
      P_axi_sea    = 0.D0
      P_axi_val    = 0.D0

*     Setting bjoken nu for K_LW calculation
      BJNU2 = DBLE(EN)**2*DBLE(Y)**2
      
*     Start to apply Bode-Yang corrections: K_vec, K_axi, H(x,Q2), K_LW

*     Bodek-Yang vector C constants

      if(nebodek.eq.1) then ! ---> Using parameters from arXiv:hep-ex/0508007
         C_vec_val1_d = 0.202D0
         C_vec_val2_d = 0.255D0
         C_vec_val1_u = 0.291D0
         C_vec_val2_u = 0.189D0
         C_vec_sea_d  = 0.621D0
         C_vec_sea_u  = 0.363D0
      else if(nebodek.ge.2) then ! ---> Using new parameters from BY2019 paper in prep.
         C_vec_val1_d = 0.341D0
         C_Vec_val2_d = 0.323D0
         C_vec_val1_u = 0.417D0
         C_vec_val2_u = 0.264D0
         C_vec_sea_d  = 0.561D0
         C_vec_sea_u  = 0.369D0
      endif

*     Updated LW scaling factors
      
      if(nebodek.ge.5.and.W.gt.W_LW_MIN) then

         C_LW = 0.218D0
         K_LW = (BJNU2 + C_LW) / BJNU2
         
      endif

      if(nebodek.eq.6) then
         K_vec_sea_c = Q2 / (Q2 + xmc*xmc)
      endif

*     Compute correction factors
         K_vec_val_u = (1.D0-1.D0/(1.D0+Q2/0.71D0)**4)*
     &    ((Q2+C_vec_val2_u)/(Q2+C_vec_val1_u))

         K_vec_val_d = (1.D0-1.D0/(1.D0+Q2/0.71D0)**4)*
     &    ((Q2+C_vec_val2_d)/(Q2+C_vec_val1_d))

         K_vec_sea_u = Q2/(Q2 + C_vec_sea_u)
         K_vec_sea_d = Q2/(Q2 + C_vec_sea_d)
*     Updated parameters separating for F2 axial parts, Bodek-Yang 2019 Paper in prep
      if(nebodek.ne.4) then   
         C_axi_sea=0.75D0
         P_axi_sea=0.55D0
         P_axi_val=0.018D0

         K_axi_sea = (Q2 + P_axi_sea * C_axi_sea)/(Q2 + C_axi_sea)
         K_axi_val = (Q2 + p_axi_val)/(Q2 + 0.18D0)
      endif

*     Updated scaling factor for xF3
      
      if(nebodek.ne.3.and.nebodek.ne.4) then
         H_factor = 0.914D0 + 0.296D0*X_H - 0.374D0*X_H*X_H
     &        + 0.165D0*X_H*X_H*X_H
      endif

      
      K_vec_val_u = K_vec_val_u*K_LW
      K_vec_val_d = K_vec_val_d*K_LW
      
      CALL QGDISG(X,Q2,G,U,D,AU,AD,S,C,B,T)

CX    TESTING SEA QUARK FACTION ON SYS ERR
C      QUARK_TOT = U + D + AU + AD + S + C + B + T
C      AU = AU*1.05
C      AD = AD*1.05
C      S = S*1.05
c      C = C*1.05
C      B = B*1.05
C      T = T*1.05
c      FRAC_VAL = U/(U+D)
C      U = (QUARK_TOT-AU-AD-S-C-B-T)*FRAC_VAL
C      D = (QUARK_TOT-AU-AD-S-C-B-T)*(1.-FRAC_VAL)

*     Contracting the expressions for quarks

      AD_F3 = 2.D0*AD
      AU_F3 = 2.D0*AU
      D_F3  = 2.D0*D + AD_F3
      U_F3  = 2.D0*U + AU_F3
      S_F3  = 2.D0*S
      C_F3  = 2.D0*C
      AS_F3 = S_F3
      AC_F3 = C_F3


      AD_v = AD
      AD_a = AD
      AU_v = AU
      AU_a = AU
      D_F2_v = D + AD
      D_F2_a = D + AD
      U_F2_v = U + AU
      U_F2_a = U + AU
      S_F2_v = S
      S_F2_a = S
      C_F2_v = C
      C_F2_a = C
      AS_v = S
      AS_a = S
      AC_v = C
      AC_a = C

      if(nebodek.ge.1) then
         AD_F3 = K_vec_sea_d*AD_F3
         AU_F3 = K_vec_sea_u*AU_F3
         D_F3  = 2.*K_vec_val_d*D + AD_F3
         U_F3  = 2.*K_Vec_val_u*U + AU_F3
         S_F3  = K_vec_sea_d*S_F3
         C_F3  = K_vec_sea_c*C_F3
         AS_F3 = S_F3
         AC_F3 = C_F3
         
         if(nebodek.eq.4) then 
            
            AD_v = 0.5D0*AD_F3
            AD_a = 0.5D0*AD_F3
            AU_v = 0.5D0*AU_F3
            AU_a = 0.5D0*AU_F3
            AS_v = 0.5D0*AS_F3
            AS_a = 0.5D0*AS_F3
            AC_v = 0.5D0*AC_F3
            AC_a = 0.5D0*AC_F3
            D_F2_v = 0.5D0*D_F3
            D_F2_a = 0.5D0*D_F3
            U_F2_v = 0.5D0*U_F3
            U_F2_a = 0.5D0*U_F3
            S_F2_v = 0.5D0*S_F3
            S_F2_a = 0.5D0*S_F3
            C_F2_v = 0.5D0*C_F3
            C_F2_a = 0.5D0*C_F3
            
         else if(nebodek.ge.2) then 
            AD_v = K_vec_sea_d*AD
            AD_a = K_axi_sea*AD
            AU_v = K_vec_sea_u*AU
            AU_a = K_axi_sea*AU
            D_F2_v = K_vec_val_d*D + AD_v
            D_F2_a = K_axi_val*D + AD_a
            U_F2_v = K_vec_val_u*U + AU_v
            U_F2_a = K_axi_val*U + AU_a
            S_F2_v = K_vec_sea_d*S
            S_F2_a = K_axi_sea*S
            C_F2_v = K_vec_sea_c*C
            C_F2_a = K_axi_sea*C
         endif
      endif
*     For sea quarks, take anti-quark probabilities to be the same as quark one
      AS_v=S_F2_v
      AS_a=S_F2_a 
      AC_v=C_F2_v
      AC_a=C_F2_a
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
Cx                  F2 = 2.D0*(D*(Vud2+Vcd2)+S*(Vus2+Vcs2)+B
Cx     &                 +AU*(Vud2+Vus2)+AC+AT)
Cx                  xF3 = 2.D0*(D*(Vud2+Vcd2)+S*(Vus2+Vcs2)+B
Cx     &                 -AU*(Vud2+Vus2)-AC-AT)

               F2v = D_F2_v*(Vud2+Vcd2)+S_F2_v*(Vus2+Vcs2)+B
     &              +AU_v*(Vud2+Vus2)+AC_v+AT
               F2a = D_F2_a*(Vud2+Vcd2)+S_F2_a*(Vus2+Vcs2)+B
     &              +AU_a*(Vud2+Vus2)+AC_a+AT
               xF3 = D_F3*(Vud2+Vcd2)+S_F3*(Vus2+Vcs2)
     &               +2.D0*B-AU_F3*(Vud2+Vus2)-AC_F3-2.D0*AT
            else
C     nu-neutron
Cx               F2 = (UTOT*(Vud2+Vcd2)+S*(Vus2+Vcs2)+B
Cx     &              +AD*(Vud2+Vus2)+AC+AT)
Cx               xF3 = (U_F3*Vud2+S_F3*Vus2+B
Cx     &               -AD_F3*(Vud2+Vus2)-AC_F3-AT)

               F2v = (U_F2_v*(Vud2+Vcd2)+S_F2_v*(Vus2+Vcs2)+B
     &              +AD_v*(Vud2+Vus2)+AC_v+AT)
               F2a = (U_F2_a*(Vud2+Vcd2)+S_F2_a*(Vus2+Vcs2)+B
     &              +AD_a*(Vud2+Vus2)+AC_a+AT)
               xF3 = (U_F3*Vud2+S_F3*Vus2+2.D0*B
     &               -AD_F3*(Vud2+Vus2)-AC_F3-2.D0*AT)

            endif
         else                   ! anti-neutrino
            if(it.eq.2212) then
C     nubar-proton
Cx               F2 = (UTOT*(Vud2+Vus2)+C+T+AD*(Vud2+Vcd2)
Cx     &              +AS*(Vus2+Vcs2)+AB)
Cx               xF3 = (U_F3*(Vud2+Vus2)+C_F3+T-AD_F3*Vud2
Cx     &              -AS_F3*Vus2-AB)
               F2v = (U_F2_v*(Vud2+Vus2)+C_F2_v+T+AD_v*(Vud2+Vcd2)
     &              +AS_v*(Vus2+Vcs2)+AB)
               F2a = (U_F2_a*(Vud2+Vus2)+C_F2_a+T+AD_a*(Vud2+Vcd2)
     &              +AS_a*(Vus2+Vcs2)+AB)
               xF3 = (U_F3*(Vud2+Vus2)+C_F3+2.D0*T-AD_F3*Vud2
     &              -AS_F3*Vus2-2.D0*AB)

            else
C     nubar-neutron
Cx               F2 = (DTOT*(Vud2+Vus2)+C+T+AU*(Vud2+Vcd2)
Cx     &              +AS*(Vus2+Vcs2)+AB)
Cx               xF3 = (D_F3*(Vud2+Vus2)+C_F3+T-AU_F3*(Vud2+Vcd2)
Cx     &               -AS_F3*(Vus2+Vcs2)-AB)
               F2v = (D_F2_v*(Vud2+Vus2)+C_F2_v+T+AU_v*(Vud2+Vcd2)
     &              +AS_v*(Vus2+Vcs2)+AB)
               F2a = (D_F2_a*(Vud2+Vus2)+C_F2_a+T+AU_a*(Vud2+Vcd2)
     &              +AS_a*(Vus2+Vcs2)+AB)
               xF3 = (D_F3*(Vud2+Vus2)+C_F3+2.D0*T-AU_F3*(Vud2+Vcd2)
     &               -AS_F3*(Vus2+Vcs2)-2.D0*AB)

            endif
         endif
      else                      ! Neutral current case
         if(ipar.gt.0) then     ! neutrino
            if(it.eq.2212) then 
C     nu-proton
Cx               F2 = (DTOT+S+B+AU+AC+AT)
Cx               xF3 = (D_F3+S_F3+B-AU_F3-AC_F3-AT)
               F2v = (D_F2_v+S_F2_v+B+AU_v+AC_v+AT)
               F2a = (D_F2_a+S_F2_a+B+AU_a+AC_a+AT)
               xF3 = (D_F3+S_F3+2.D0*B-AU_F3-AC_F3-2.D0*AT)
            else
C     nu-neutron
Cx               F2 = (UTOT+S+B+AD+AC+AT)
Cx               xF3 = (U_F3+S_F3+B-AD_F3-AC_F3-AT)
               F2v = (U_F2_v+S_F2_v+B+AD_v+AC_v+AT)
               F2a = (U_F2_a+S_F2_a+B+AD_a+AC_a+AT)
               xF3 = (U_F3+S_F3+2.D0*B-AD_F3-AC_F3-2.D0*AT)
            endif
         else                   ! anti-neutrino
            if(it.eq.2212) then
C     nubar-proton
Cx               F2 = (UTOT+C+T+AD+AS+AB)
Cx               xF3 = (U_F3+C_F3+T-AD_F3-AS_F3-AB)
               F2v = (U_F2_v+C_F2_v+T+AD_v+AS_v+AB)
               F2a = (U_F2_a+C_F2_a+T+AD_a+AS_a+AB)
               xF3 = (U_F3+C_F3+2.D0*T-AD_F3-AS_F3-2.D0*AB)
            else
C     nubar-neutron
Cx               F2 = (DTOT+C+T+AU+AS+AB)
Cx               xF3 = (D_F3+C_F3+T-AU_F3-AS_F3-AB)
               F2v = (D_F2_v+C_F2_v+T+AU_v+AS_v+AB)
               F2a = (D_F2_a+C_F2_a+T+AU_a+AS_a+AB)
               xF3 = (D_F3+C_F3+2.D0*T-AU_F3-AS_F3-2.D0*AB)
            endif
         endif
      endif

      xF3 = K_vec_sea_c*H_factor*xF3
            
      RETURN
      END
