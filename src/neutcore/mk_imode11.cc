/********************************************************************************************
*                                                                                           *
*     ( purpose )                                                                           *
*       return differental cross-section for \nu n -> \mu n \pi- channel (IMODE11)           *
*              for specific W, Q, neutrino energy & pion angles in the Adler frame          *
*                                                                                           *
*     ( input )                                                                             *
*       Q           : sqaured 4-momentum transfer, Q = -k*k                                 *
*       W           : Invariant mass                                                        *
*       E           : Energy of incident neutrino in the lab frame (GeV)                    *
*       t     : pion polar angle in Adler frame                                             *
*       phi   : pion Azimuthal angle in Adler frame                                         *
*                                                                                           *
*     ( output )                                                                            *
*      d^4(sig)/dwdQ dt d(phi)                                                              *
*                                                                                           *
*     ( creation date and author )                                                          *
*     Minoo Kabirnezhad, March 2019                                                         *
*                                                                                           *
*     ( comment )                                                                           *
*       Reference: M. Kabirnezhad, Phys.Rev.D 97 (2018) no.1, 013002                        *
*                  M. Kabirnezhad thesis:                                                   *
*                  https://www.ncbj.gov.pl/sites/default/files/m.kabirnezhad_thesis_0.pdf   *
*                  Rein and Sehgal,Ann. of Phys. 133(1981),79-153                           *
*                  D. Rein ,Z.Phys.C 35(1987),43-64                                         *
*                  S. L. Adler, Ann. Phys. (N.Y.) 50, 189 (1968).                           *
********************************************************************************************/
#include "mkcons.h"

extern "C"{

void xsec_im11_(float *diff_Xsec, float *Wf, float *Qf, float *Ef, float *tf, float *phif, float *lepmass) {

  double W = *Wf;
  double Q = *Qf;
  double E = *Ef;
  double t = *tf;
  double phi = *phif;
  double m_l= *lepmass;

  MAR = nemdls_.xmaspi;
  M_A = nemdls_.xmabkgm;
  M_V = nemdls_.xmvspi;
  C50[0] = nemdls_.rca5ispi;

  //kinematic
  double q_0 =(W*W - M*M + m_pi* m_pi)/(2.*W);
  double k_0 = (W*W - M*M - Q)/(2.*W) ;
  double abs_mom_q = sqrt_chk( q_0*q_0 - m_pi*m_pi); 
  double abs_mom_k = sqrt_chk(k_0*k_0 + Q) ;
  double E_2L = (M*M - W*W - Q + 2.*M*E)/(2.*M) ;
  double abs_mom_k_L = (W/M)*abs_mom_k ;
  double k_2 = (M*M + 2.*M*E - W*W -m_l*m_l)/(2.*W) ;
  double k_1 =  k_2 + (W*W - M*M - Q)/(2.*W) ;
  double p_10 = (W*W + M*M + Q)/(2.*W) ;
  double qk = q_0*k_0 - abs_mom_k*abs_mom_q*t ;
  double k_2L = sqrt_chk(E_2L*E_2L - m_l*m_l) ; //magnitude of lepton momentum in lab frame
  double k_2_iso = sqrt_chk(k_2*k_2 - m_l*m_l) ; //magnitude of lepton momentum in CM frame
  ///////////////////////////////////////////////////////////////////////////////////////////
  double cos_theta = (2.*k_1*k_2  - Q - m_l*m_l)/( 2.*k_1*k_2_iso ) ;
  double A_plus= sqrt_chk( k_1*(k_2 - k_2_iso) );
  double A_minus= sqrt_chk( k_1*(k_2 + k_2_iso) );
  double eps_1_plus = A_plus*((k_1 - k_2_iso)/abs_mom_k)*sqrt_chk(1.+ cos_theta);
  double eps_1_minus= -A_minus*((k_1 + k_2_iso)/abs_mom_k)*sqrt_chk(1.- cos_theta);
  double eps_2_plus = A_plus*sqrt_chk(1.+ cos_theta);
  double eps_2_minus= A_minus*sqrt_chk(1.- cos_theta);

  double C_L_plus =  (abs_mom_k_L/(2.*E))*(eps_1_plus  - eps_2_plus);
  double C_L_minus=  (abs_mom_k_L/(2.*E))*(eps_1_minus - eps_2_minus);
  double C_R_plus = -(abs_mom_k_L/(2.*E))*(eps_1_plus  + eps_2_plus);
  double C_R_minus= -(abs_mom_k_L/(2.*E))*(eps_1_minus + eps_2_minus);
  double C_T_minus_square = C_L_minus*C_L_minus + C_L_plus*C_L_plus;   
  double C_T_plus_square= C_R_minus *C_R_minus + C_R_plus*C_R_plus ;
  double C_LR=-(abs_mom_k_L/(2.*E))*(abs_mom_k_L/(2.*E))*(eps_1_minus*eps_1_minus + eps_1_plus*eps_1_plus - eps_2_minus*eps_2_minus - eps_2_plus*eps_2_plus);
  double C_S_plus_square = //equal to zero in massless lepton
    (abs_mom_k_L*abs_mom_k_L/ (2.*E*E))*( k_1*(k_2 - k_2_iso ))*(1. - cos_theta);
  double C_S_minus_square = //equal to sqrt_chk(2.uv) is massless lepton
    (abs_mom_k_L*abs_mom_k_L/ (2.*E*E))*( k_1*(k_2 + k_2_iso )) *(1. + cos_theta) ;

  double eps_zero_L = -1.;
  double eps_zero_R = 1.;
  double eps_z_L =  -( (k_1 - k_2_iso)/abs_mom_k );
  double eps_z_R = ( (k_1 + k_2_iso)/abs_mom_k ) ; 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /*    Bkg Contribution   */

  //Auxiliary functions
  double W_plus = W+M ;
  double W_minus = W-M ;
  double O_1_plus = sqrt_chk((W_plus*W_plus + Q)*( W_plus*W_plus - m_pi*m_pi ))/ (2.*W) ;
  double O_1_minus = sqrt_chk((W_minus*W_minus + Q)*(W_minus*W_minus - m_pi*m_pi))/ (2.*W) ;
  double O_2_plus = sqrt_chk((W_plus*W_plus + Q)/(W_plus*W_plus - m_pi*m_pi)) ;
  double O_2_minus = sqrt_chk((W_minus*W_minus + Q)/(W_minus*W_minus - m_pi*m_pi)) ;
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Vector_helicity amplituds for bkg
  double K_1_V = W_minus*O_1_plus; 
  double K_2_V = W_plus*O_1_minus; 
  double K_3_V = (q_0*q_0 - m_pi*m_pi)*W_plus*O_2_minus; 
  double K_4_V = abs_mom_q*abs_mom_q*W_minus*O_2_plus; 
  double K_5_V = 1./O_2_plus; 
  double K_6_V = 1./O_2_minus; 
  double FV_cut ; // vector cut
  if (W<1.30)
    FV_cut = 1. ;//virtual form factor to kill the bkg smothly
  else if (W>=1.30 && W<1.6)
    FV_cut =   8.08497*W*W*W  -41.6767*W*W + 66.3419*W -32.5733   ;
  else 
    FV_cut = 0.;
  // vector form-factor
  double mup= 2.792847;
  double mun= -1.913043;
  double Ln= 5.6;
  double tau= Q/(4.*M*M);
  double GEp= 1./((1.+ Q/(M_V*M_V))*(1.+ Q/(M_V*M_V))) ;
  double GEn= -(mun*tau/(1 + Ln*tau))*GEp;
  double GMn= mun*GEp;
  double GMp= mup*GEp;
  double F1n= (GEn + tau*GMn)/(1.+tau);
  double F1p= (GEp + tau*GMp)/(1.+tau);
  double muF2n= (GMn - GEn)/(1.+tau);
  double muF2p= (GMp - GEp)/(1.+tau);

  double FF_1 = 0.5*(F1p - F1n);             
  double FF_2 = 0.5*(muF2p - muF2n);
  double FF_pi = 2.*FF_1 ;

  double V_1 =  sqrt_chk(2.)*(g_A*M/f_pi)*(2.*FF_1*FV_cut)*(1./(m_pi*m_pi - 2.*q_0*p_10 - 2.*abs_mom_q*abs_mom_k*t)) +  sqrt_chk(2.)*(g_A/f_pi)*(FF_2*FV_cut/M) ;
  double V_2 =  sqrt_chk(2.)*(g_A*M/f_pi)*(2.*FF_1*FV_cut/qk)*(1./(m_pi*m_pi - 2.*q_0*p_10 - 2.*abs_mom_q*abs_mom_k*t)) ;
  double V_3 =  sqrt_chk(2.)*(g_A*M/f_pi)*(2.*FF_2*FV_cut/M)*(1./(m_pi*m_pi - 2.*q_0*p_10 - 2.*abs_mom_q*abs_mom_k*t)) ;
  double V_4 = -sqrt_chk(2.)*(g_A*M/f_pi)*(2.*FF_2*FV_cut/M)*(1./(m_pi*m_pi - 2.*q_0*p_10 - 2.*abs_mom_q*abs_mom_k*t)) ;
  double V_5 = -sqrt_chk(2.)*(g_A*M/f_pi)*(FF_pi*FV_cut/qk)*(1./(Q+2.*qk)) ;

  double V_25 = (W_plus*W_minus)*V_2 - Q*V_5; 

  double F_1 = V_1 + qk*(V_3 - V_4)/W_minus + W_minus*V_4; 
  double F_2 = -V_1 + qk*(V_3 - V_4)/W_plus  + W_plus*V_4; 
  double F_3 = (V_3 - V_4) + V_25/W_plus; 
  double F_4 = (V_3 - V_4) - V_25/W_minus; 
  double F_5 = (W_plus*W_plus + Q)*V_1/(2.*W) - qk*(W_plus*W_plus + Q + 2.*W*W_minus)*V_2/(2.*W) + (W_plus*q_0 
      - qk)*(V_3 - V_4) + (W_plus*W_plus + Q)*W_minus*V_4/(2.*W) - k_0*qk*V_5 + q_0*V_25; 
  double F_6 = -(W_minus*W_minus + Q)*V_1/(2.*W) + qk*(W_plus*W_plus + Q + 2.*W*W_minus)*V_2/(2.*W) + (W_minus*q_0 
      - qk)*(V_3 - V_4) + (W_minus*W_minus + Q)*W_plus*V_4/(2.*W) - q_0*V_25  + k_0*qk*V_5; 

  double sF_1 = K_1_V*F_1/(2.*M); 
  double sF_2 = K_2_V*F_2/(2.*M); 
  double sF_3 = K_3_V*F_3/(2.*M); 
  double sF_4 = K_4_V*F_4/(2.*M); 
  double sF_5 = K_5_V*F_5/(2.*M); 
  double sF_6 = K_6_V*F_6/(2.*M); 

  double F_plus11 = -sqrt_chk(2.)*( sqrt_chk((1-t)/2.)*(sF_1 + sF_2) + 0.5*sqrt_chk(1. - t*t)* sqrt_chk((1+t)/2.)*(sF_3 + sF_4) );
  double F_plus_11 = -sqrt_chk(2.)*( sqrt_chk((1+t)/2.)*(sF_1 - sF_2) - 0.5*sqrt_chk(1. - t*t)* sqrt_chk((1-t)/2.)*(sF_3 - sF_4) );
  double F_plus1_1 = (1./sqrt_chk(2.))*sqrt_chk(1. - t*t)* sqrt_chk((1-t)/2.)*(sF_3 - sF_4);
  double F_plus_1_1 = (1./sqrt_chk(2.))*sqrt_chk(1. - t*t)* sqrt_chk((1+t)/2.)*(sF_3 + sF_4);

  double F_minus11 = (1./sqrt_chk(2.))*sqrt_chk(1. - t*t)* sqrt_chk((1+t)/2.)*(sF_3 + sF_4);
  double F_minus_11 = -(1./sqrt_chk(2.))*sqrt_chk(1. - t*t)* sqrt_chk((1-t)/2.)*(sF_3 - sF_4);
  double F_minus1_1 = sqrt_chk(2.)*( sqrt_chk((1+t)/2.)*(sF_1 - sF_2) - 0.5*sqrt_chk(1. - t*t)* sqrt_chk((1-t)/2.)*(sF_3 - sF_4) );
  double F_minus_1_1 = -sqrt_chk(2.)*( sqrt_chk((1-t)/2.)*(sF_1 + sF_2) + 0.5*sqrt_chk(1. - t*t)* sqrt_chk((1+t)/2.)*(sF_3 + sF_4) );

  double F_zero_minus11 = sqrt_chk((1+t)/2.)*(k_0*eps_z_L - abs_mom_k*eps_zero_L)*(sF_5 + sF_6); 
  double F_zero_minus_11 = - sqrt_chk((1-t)/2.)*(k_0*eps_z_L - abs_mom_k*eps_zero_L)*(sF_5 - sF_6); 
  double F_zero_minus1_1 = - sqrt_chk((1-t)/2.)*(k_0*eps_z_L - abs_mom_k*eps_zero_L)*(sF_5 - sF_6);
  double F_zero_minus_1_1 = - sqrt_chk((1+t)/2.)*(k_0*eps_z_L - abs_mom_k*eps_zero_L)*(sF_5 + sF_6);

  double F_zero_plus11 = sqrt_chk((1+t)/2.)*(k_0*eps_z_R - abs_mom_k*eps_zero_R)*(sF_5 + sF_6); 
  double F_zero_plus_11 = - sqrt_chk((1-t)/2.)*(k_0*eps_z_R - abs_mom_k*eps_zero_R)*(sF_5 - sF_6);
  double F_zero_plus1_1 = - sqrt_chk((1-t)/2.)*(k_0*eps_z_R - abs_mom_k*eps_zero_R)*(sF_5 - sF_6);
  double F_zero_plus_1_1 = - sqrt_chk((1+t)/2.)*(k_0*eps_z_R - abs_mom_k*eps_zero_R)*(sF_5 + sF_6);
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Axial helicity amp for bkg
  double K_1_A = abs_mom_q*O_2_plus ;
  double K_2_A = abs_mom_q* O_2_minus ;
  double K_3_A = abs_mom_q*O_1_minus ;
  double K_4_A = abs_mom_q* O_1_plus ;
  double K_5_A = O_1_minus ;
  double K_6_A = O_1_plus ;
  double K_7_A = O_1_minus ;
  double K_8_A = O_1_plus ;
  double Delta = k_0*(q_0*k_0 - qk)/(abs_mom_k*abs_mom_k) ;

  double FA_cut ;
  if (W<1.28)
    FA_cut = 1. ;//virtual form factor to kill the bkg smothly
  else if (W>=1.28 && W<1.5)
    FA_cut =   45.7846*W*W*W  -185.994*W*W + 246.608*W -105.945   ;
  else 
    FA_cut = 0.; 

  double FF_A =  F0_A/((1+ Q/(M_A*M_A))*(1+ Q/(M_A*M_A))) ;
  double F_rho = F0_rho/(1. + ((Q - m_pi*m_pi + 2.*qk)/(m_rho*m_rho)) ) ;

  double A_1 =  sqrt_chk(2.)*(g_A*M/f_pi)*FF_A*FA_cut*( 1./(m_pi*m_pi - 2.*q_0*p_10 - 2.*abs_mom_q*abs_mom_k*t)) ;
  double A_3 = -sqrt_chk(2.)*(g_A*M/f_pi)*FF_A*FA_cut*( 1./(m_pi*m_pi - 2.*q_0*p_10 - 2.*abs_mom_q*abs_mom_k*t)) ;
  double A_4 = (1./sqrt_chk(2.))*(1./f_pi)*(FA_cut/M)*(-g_A*FF_A + F_rho) ;
  double A_7 = -sqrt_chk(2.)*(g_A*M/f_pi)*FF_A*FA_cut*(1./(Q + m_pi*m_pi)) ;
  double A_8 = (1./sqrt_chk(2.))*(g_A/f_pi)*FF_A*FA_cut*(1./(Q + m_pi*m_pi))*(1. +  4.*M*M/(m_pi*m_pi - 2.*q_0*p_10 - 2.*abs_mom_q*abs_mom_k*t))
    -(1./sqrt_chk(2.))*(1./f_pi)*(1./(Q + m_pi*m_pi))*F_rho*FA_cut;

  double G_1 = W_plus* A_1 - M*A_4 ;
  double G_2 = (-1)*W_minus*A_1 - M*A_4 ;
  double G_3 = A_1 - A_3 ;
  double G_4 = -A_1 + A_3 ;
  double G_5 =
    (Delta + (W_plus*W_plus - m_pi*m_pi)/(2.*W) + 2.*(W*k_0*W_plus)/(W_minus*W_minus + Q))*A_1 
    + (q_0 - Delta)*A_3 - M*W_minus*A_4/(p_10-M);
  double G_6 =
    -(Delta + (W_minus*W_minus - m_pi*m_pi)/(2.*W) + 2.*W*k_0*W_minus/(W_plus*W_plus + Q))*A_1 
    - (q_0 - Delta)*A_3 - M*W_plus*A_4/(p_10 + M);
  double G_7=  ((W_plus*W_plus - m_pi*m_pi)*A_1/(2.*W)) + q_0*A_3 - M*A_4 + k_0*A_7 + W_plus* k_0*A_8 ;
  double G_8 = -((W_minus*W_minus - m_pi*m_pi)*A_1/(2.*W)) - q_0*A_3 - M*A_4 - k_0*A_7 + W_minus* k_0*A_8 ;

  //Isobaric amplitud (Scribal G)
  double sG_1 = K_1_A*G_1/(2.*M) ;
  double sG_2 = K_2_A*G_2/(2.*M) ;
  double sG_3 = K_3_A*G_3/(2.*M) ;
  double sG_4 = K_4_A*G_4/(2.*M) ;
  double sG_5 = K_5_A*G_5/(2.*M) ;
  double sG_6 = K_6_A*G_6/(2.*M) ;
  double sG_7 = K_5_A*G_7/(2.*M) ;
  double sG_8 = K_6_A*G_8/(2.*M) ;

  // helicity amplitudes
  double G_plus11 = -sqrt_chk(2.)*(  sqrt_chk((1-t)/2.)*(sG_1 + sG_2) + 0.5*sqrt_chk(1. - t*t)* sqrt_chk((1+t)/2.)*(sG_3 + sG_4) );
  double G_plus_11 = sqrt_chk(2.)*(  sqrt_chk((1+t)/2.)*(sG_1 - sG_2) - 0.5*sqrt_chk(1. - t*t)* sqrt_chk((1-t)/2.)*(sG_3 - sG_4) );
  double G_plus1_1 = (1./sqrt_chk(2.))*sqrt_chk(1. - t*t)* sqrt_chk((1-t)/2.)*(sG_3 - sG_4);
  double G_plus_1_1 = -(1./sqrt_chk(2.))*sqrt_chk(1. - t*t)* sqrt_chk((1+t)/2.)*(sG_3 + sG_4);

  double G_minus11 =  (1./sqrt_chk(2.))*sqrt_chk(1. - t*t)* sqrt_chk((1+t)/2.)*(sG_3 + sG_4);
  double G_minus_11 =  (1./sqrt_chk(2.))*sqrt_chk(1. - t*t)* sqrt_chk((1-t)/2.)*(sG_3 - sG_4);
  double G_minus1_1 =  sqrt_chk(2.)*(  sqrt_chk((1+t)/2.)*(sG_1 - sG_2) - 0.5*sqrt_chk(1. - t*t)* sqrt_chk((1-t)/2.)*(sG_3 - sG_4) );
  double G_minus_1_1 =  sqrt_chk(2.)*(  sqrt_chk((1-t)/2.)*(sG_1 + sG_2) + 0.5*sqrt_chk(1. - t*t)* sqrt_chk((1+t)/2.)*(sG_3 + sG_4) );

  double G_zero_minus11 = ( sqrt_chk((1+t)/2.)/k_0)*( (abs_mom_k*eps_z_L)*(sG_5 + sG_6) + (k_0*eps_zero_L - abs_mom_k*eps_z_L)*(sG_7 + sG_8) );
  double G_zero_minus_11 = ( sqrt_chk((1-t)/2.)/k_0)*( (abs_mom_k*eps_z_L)*(sG_5 - sG_6) + (k_0*eps_zero_L - abs_mom_k*eps_z_L)*(sG_7 - sG_8) );
  double G_zero_minus1_1 = -( sqrt_chk((1-t)/2.)/k_0)*( (abs_mom_k*eps_z_L)*(sG_5 - sG_6) + (k_0*eps_zero_L - abs_mom_k*eps_z_L)*(sG_7 - sG_8) );
  double G_zero_minus_1_1 = ( sqrt_chk((1+t)/2.)/k_0)*( (abs_mom_k*eps_z_L)*(sG_5 + sG_6) + (k_0*eps_zero_L - abs_mom_k*eps_z_L)*(sG_7 + sG_8) ); 

  double G_zero_plus11 = ( sqrt_chk((1+t)/2.)/k_0)*( (abs_mom_k*eps_z_R)*(sG_5 + sG_6) + (k_0*eps_zero_R - abs_mom_k*eps_z_R)*(sG_7 + sG_8) );
  double G_zero_plus_11 =  ( sqrt_chk((1-t)/2.)/k_0)*( (abs_mom_k*eps_z_R)*(sG_5 - sG_6) + (k_0*eps_zero_R - abs_mom_k*eps_z_R)*(sG_7 - sG_8) );
  double G_zero_plus1_1 = -( sqrt_chk((1-t)/2.)/k_0)*( (abs_mom_k*eps_z_R)*(sG_5 - sG_6) + (k_0*eps_zero_R - abs_mom_k*eps_z_R)*(sG_7 - sG_8) );
  double G_zero_plus_1_1 =  ( sqrt_chk((1+t)/2.)/k_0)*( (abs_mom_k*eps_z_R)*(sG_5 + sG_6) + (k_0*eps_zero_R - abs_mom_k*eps_z_R)*(sG_7 + sG_8) ); 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /*   Resonace contribution   */

  double N = (W+M)*(W+M) + Q ; 
  double H = M/W ;
  double X = sqrt_chk(2./Omega) ; 
  double L = X*H*abs_mom_k_L ;

  double qMR_0[nRes], abs_mom_qMR[nRes], Gamma[nRes], Vf_BW_real[nRes], Vf_BW_Im[nRes], Af_BW_real[nRes], Af_BW_Im[nRes];

  double  HV_realM11[nRes], HV_realM_11[nRes],  HV_realM1_1[nRes],  HV_realM_1_1[nRes] ;    
  double  HV_realP11[nRes], HV_realP_11[nRes],  HV_realP1_1[nRes],  HV_realP_1_1[nRes] ;   
  double  HV0_realM11[nRes],HV0_realM_11[nRes], HV0_realM1_1[nRes], HV0_realM_1_1[nRes] ;    
  double  HV0_realP11[nRes],HV0_realP_11[nRes], HV0_realP1_1[nRes], HV0_realP_1_1[nRes] ;    

  double  HA_realM11[nRes],HA_realM_11[nRes], HA_realM1_1[nRes], HA_realM_1_1[nRes] ;    
  double  HA_realP11[nRes],HA_realP_11[nRes], HA_realP1_1[nRes], HA_realP_1_1[nRes] ;   
  double  HA0_realM11[nRes],HA0_realM_11[nRes], HA0_realM1_1[nRes], HA0_realM_1_1[nRes] ;    
  double  HA0_realP11[nRes],HA0_realP_11[nRes], HA0_realP1_1[nRes], HA0_realP_1_1[nRes] ;

  double  HV_ImM11[nRes], HV_ImM_11[nRes],  HV_ImM1_1[nRes],  HV_ImM_1_1[nRes] ;    
  double  HV_ImP11[nRes], HV_ImP_11[nRes],  HV_ImP1_1[nRes],  HV_ImP_1_1[nRes] ;    
  double  HV0_ImM11[nRes],HV0_ImM_11[nRes], HV0_ImM1_1[nRes], HV0_ImM_1_1[nRes] ;    
  double  HV0_ImP11[nRes],HV0_ImP_11[nRes], HV0_ImP1_1[nRes], HV0_ImP_1_1[nRes] ;    

  double  HA_ImM11[nRes],HA_ImM_11[nRes], HA_ImM1_1[nRes], HA_ImM_1_1[nRes] ;    
  double  HA_ImP11[nRes],HA_ImP_11[nRes], HA_ImP1_1[nRes], HA_ImP_1_1[nRes] ;    
  double  HA0_ImM11[nRes],HA0_ImM_11[nRes], HA0_ImM1_1[nRes], HA0_ImM_1_1[nRes] ;    
  double  HA0_ImP11[nRes],HA0_ImP_11[nRes], HA0_ImP1_1[nRes], HA0_ImP_1_1[nRes] ; 

  double fV1, fV_1, fV3, fV_3, fV0_minus,fV0_plus, fV0R_minus,fV0R_plus;
  double fA1, fA_1, fA3, fA_3, fA0_minus,fA0_plus, fA0R_minus,fA0R_plus;
  double kapa[nRes],CA5[nRes], G_A0[nRes], G_A1[nRes], G_A2[nRes], CV4[nRes], G_V0[nRes], G_V1[nRes], G_V2[nRes]  ;
  double C0_minus[nRes], C1_minus[nRes], C2_minus[nRes], B1_minus[nRes], B2_minus[nRes] ;
  double C0_plus[nRes], C1_plus[nRes], C2_plus[nRes], B1_plus[nRes], B2_plus[nRes];
  double R0_A[nRes], R1_A[nRes], R2_A[nRes], R0_V[nRes], R1_V[nRes], R2_V[nRes];
  double T0_A[nRes], T1_A[nRes], T2_A[nRes], T0_V[nRes], T1_V[nRes], T2_V[nRes];
  double S1_KLM_minus[nRes], S2_KLM_minus[nRes], S1_KLM_plus[nRes], S2_KLM_plus[nRes];

  for (int i=0; i<nRes; i++){
    qMR_0[i] =(MR[i]*MR[i] - M*M + m_pi* m_pi)/(2.*MR[i]);
    abs_mom_qMR[i] = sqrt_chk( qMR_0[i]*qMR_0[i] - m_pi*m_pi);
    kapa[i] = (pi*W/M)*sqrt_chk(2./(JP[i]*abs_mom_q)) ;
    Gamma[i] = RWA[i]*pow((abs_mom_q/ abs_mom_qMR[i]), LP[i]);

    Vf_BW_real[i] = cos(qV[i])*sqrt_chk(BRA[i]*Gamma[i]/(2.*pi))*((W-MR[i])/((W-MR[i])*(W-MR[i]) + Gamma[i]*Gamma[i]/4.) )
      + sin(qV[i])*sqrt_chk(BRA[i]*Gamma[i]/(2.*pi))*((Gamma[i]/2.)/( (W-MR[i])*(W-MR[i]) + Gamma[i]*Gamma[i]/4.) );
    Vf_BW_Im[i] = - cos(qV[i])*sqrt_chk(BRA[i]*Gamma[i]/(2.*pi))*((Gamma[i]/2.)/( (W-MR[i])*(W-MR[i]) + Gamma[i]*Gamma[i]/4.) )
      + sin(qV[i])*sqrt_chk(BRA[i]*Gamma[i]/(2.*pi))*((W-MR[i])/( (W-MR[i])*(W-MR[i]) + Gamma[i]*Gamma[i]/4.) ) ;

    Af_BW_real[i] = cos(qA[i])*sqrt_chk(BRA[i]*Gamma[i]/(2.*pi))*((W-MR[i])/((W-MR[i])*(W-MR[i]) + Gamma[i]*Gamma[i]/4.) )
      + sin(qA[i])*sqrt_chk(BRA[i]*Gamma[i]/(2.*pi))*((Gamma[i]/2.)/( (W-MR[i])*(W-MR[i]) + Gamma[i]*Gamma[i]/4.) );
    Af_BW_Im[i] = - cos(qA[i])*sqrt_chk(BRA[i]*Gamma[i]/(2.*pi))*((Gamma[i]/2.)/( (W-MR[i])*(W-MR[i]) + Gamma[i]*Gamma[i]/4.) )
      + sin(qA[i])*sqrt_chk(BRA[i]*Gamma[i]/(2.*pi))*((W-MR[i])/( (W-MR[i])*(W-MR[i]) + Gamma[i]*Gamma[i]/4.) ) ;

    //GSK Form-factor
    CV4[i] = 1.51*CV[i]/( (1.+ Q/(4.*M_V*M_V))*(1.+ Q/(M_V*M_V))*(1.+ Q/(M_V*M_V)) );

    G_V0[i] =(1./(4.*sqrt_chk(3.)))*((M+W)/M)*((M+W)/M)*(1. + Q/((M +W)*(M+W)))* sqrt_chk(1. + Q/((M +W)*(M+W)))* CV4[i];
    G_V1[i] =(1./(4.*sqrt_chk(3.)))*((M+W)/M)*((M+W)/M)*(1. + Q/((M +W)*(M+W)))* sqrt_chk(1. + Q/((M +W)*(M+W)))* CV4[i] *sqrt_chk(1./(1.+ Q/(4.*M*M)));
    G_V2[i] =(1./(4.*sqrt_chk(3.)))*((M+W)/M)*((M+W)/M)*(1. + Q/((M +W)*(M+W)))* sqrt_chk(1. + Q/((M +W)*(M+W)))* CV4[i] *(1./(1.+ Q/(4.*M*M)));
    //(2) GR Form-factor
    CA5[i]= C50[i]/((1+ Q/(MAR*MAR))*(1+ Q/(MAR*MAR)) );

    G_A0[i] = 0.5*sqrt_chk(3.)*sqrt_chk(1.+ Q/(W_plus*W_plus))*((W*W -Q -M*M)/(2.*W*(W-M)) + W*Q/(4.*M*M*(W-M)))*CA5[i];
    G_A1[i] = 0.5*sqrt_chk(3.)*sqrt_chk(1.+ Q/(W_plus*W_plus))*(1 - (W*W -Q -M*M)/(8.*M*M))*CA5[i]*sqrt_chk(1./(1.+ Q/(4.*M*M)));
    G_A2[i] = 0.5*sqrt_chk(3.)*sqrt_chk(1.+ Q/(W_plus*W_plus))*(1 - (W*W -Q -M*M)/(8.*M*M))*CA5[i]*(1./(1.+ Q/(4.*M*M)));

    double a_aux = 1. + ((W*W + Q + M*M)/(2.*M*W)) ;
    //We multiply C,Band S to abs_mom_k
    C0_minus[i] = ( (eps_zero_L*abs_mom_k - eps_z_L*k_0)*((1./3.) + k_0/(a_aux*M)) + eps_z_L*((2./3.)*W - Q/(a_aux*M) ) + 
        (eps_zero_L*k_0 - eps_z_L*abs_mom_k)*abs_mom_k*(((2./3.)*W - Q/(a_aux*M))/(m_pi*m_pi + Q)) )*G_A0[i]/(2.*W);
    C1_minus[i] = ( (eps_zero_L*abs_mom_k - eps_z_L*k_0)*((1./3.) + k_0/(a_aux*M)) + eps_z_L*((2./3.)*W - Q/(a_aux*M) + Omega/(3.*a_aux*M) ) + 
        (eps_zero_L*k_0 - eps_z_L*abs_mom_k)*abs_mom_k*(((2./3.)*W - Q/(a_aux*M) + Omega/(3.*a_aux*M))/(m_pi*m_pi + Q)) )*G_A1[i]/(2.*W);
    C2_minus[i] = ( (eps_zero_L*abs_mom_k - eps_z_L*k_0)*((1./3.) + k_0/(a_aux*M)) + eps_z_L*((2./3.)*W - Q/(a_aux*M) + 2.*Omega/(3.*a_aux*M) ) + 
        (eps_zero_L*k_0 - eps_z_L*abs_mom_k)*abs_mom_k*(((2./3.)*W - Q/(a_aux*M) + 2.*Omega/(3.*a_aux*M))/(m_pi*m_pi + Q)) )*G_A2[i]/(2.*W);

    C0_plus[i] = ( (eps_zero_R*abs_mom_k - eps_z_R*k_0)*((1./3.) + k_0/(a_aux*M)) + eps_z_R*((2./3.)*W - Q/(a_aux*M) ) +   
        (eps_zero_R*k_0 - eps_z_R*abs_mom_k)*abs_mom_k*(((2./3.)*W - Q/(a_aux*M))/(m_pi*m_pi + Q)) )*G_A0[i]/(2.*W) ;
    C1_plus[i] = ( (eps_zero_R*abs_mom_k - eps_z_R*k_0)*((1./3.) + k_0/(a_aux*M)) + eps_z_R*((2./3.)*W - Q/(a_aux*M) + Omega/(3.*a_aux*M) ) +   
        (eps_zero_R*k_0 - eps_z_R*abs_mom_k)*abs_mom_k*(((2./3.)*W - Q/(a_aux*M) + Omega/(3.*a_aux*M))/(m_pi*m_pi + Q)) )*G_A1[i]/(2.*W) ;                     
    C2_plus[i] =  ( (eps_zero_R*abs_mom_k - eps_z_R*k_0)*((1./3.) + k_0/(a_aux*M)) + eps_z_R*((2./3.)*W - Q/(a_aux*M) + 2.*Omega/(3.*a_aux*M) ) +   
        (eps_zero_R*k_0 - eps_z_R*abs_mom_k)*abs_mom_k*(((2./3.)*W - Q/(a_aux*M) + 2.*Omega/(3.*a_aux*M))/(m_pi*m_pi + Q)) )*G_A2[i]/(2.*W) ;                     

    B1_minus[i] = sqrt_chk(Omega/2)* (eps_zero_L + eps_z_L*abs_mom_k/(a_aux*M))*G_A1[i]/(3.*W)+ 
      sqrt_chk(Omega/2)* (eps_zero_L*k_0 - eps_z_L*abs_mom_k)* ((k_0 + abs_mom_k*abs_mom_k/(2.*M*a_aux))/(m_pi*m_pi + Q)) *G_A1[i]/(3.*W) ;
    B2_minus[i] = sqrt_chk(Omega/2)* (eps_zero_L + eps_z_L*abs_mom_k/(a_aux*M))*G_A2[i]/(3.*W)+ 
      sqrt_chk(Omega/2)* (eps_zero_L*k_0 - eps_z_L*abs_mom_k)* ((k_0 + abs_mom_k*abs_mom_k/(2.*M*a_aux))/(m_pi*m_pi + Q)) *G_A2[i]/(3.*W) ;

    B1_plus[i] = sqrt_chk(Omega/2)* (eps_zero_R + eps_z_R*abs_mom_k/(a_aux*M))*G_A1[i]/(3.*W)+ 
      sqrt_chk(Omega/2)* (eps_zero_R*k_0 - eps_z_R*abs_mom_k)* ((k_0 + abs_mom_k*abs_mom_k/(2.*M*a_aux))/(m_pi*m_pi + Q)) *G_A1[i]/(3.*W) ;
    B2_plus[i] = sqrt_chk(Omega/2)* (eps_zero_R + eps_z_R*abs_mom_k/(a_aux*M))*G_A2[i]/(3.*W)+ 
      sqrt_chk(Omega/2)* (eps_zero_R*k_0 - eps_z_R*abs_mom_k)* ((k_0 + abs_mom_k*abs_mom_k/(2.*M*a_aux))/(m_pi*m_pi + Q)) *G_A2[i]/(3.*W) ;

    S1_KLM_minus[i] = abs_mom_k*(eps_z_L*k_0 - eps_zero_L*abs_mom_k)* (1. + Q/(M*M) - 3.*W/M)*( G_V1[i]/(6.*abs_mom_k_L* abs_mom_k_L)) ;
    S2_KLM_minus[i] = abs_mom_k*(eps_z_L*k_0 - eps_zero_L*abs_mom_k)* (1. + Q/(M*M) - 3.*W/M)*( G_V2[i]/(6.*abs_mom_k_L* abs_mom_k_L)) ;

    S1_KLM_plus[i] = abs_mom_k*(eps_z_R*k_0 - eps_zero_R*abs_mom_k)* (1. + Q/(M*M) - 3.*W/M)*( G_V1[i]/(6.*abs_mom_k_L* abs_mom_k_L))  ;
    S2_KLM_plus[i] = abs_mom_k*(eps_z_R*k_0 - eps_zero_R*abs_mom_k)* (1. + Q/(M*M) - 3.*W/M)*( G_V2[i]/(6.*abs_mom_k_L* abs_mom_k_L))  ;

    R0_V[i] = -sqrt_chk(2.)*H*( (W+M)*abs_mom_k_L/N )* G_V0[i] ;
    R1_V[i] = -sqrt_chk(2.)*H*( (W+M)*abs_mom_k_L/N )* G_V1[i] ;
    R2_V[i] = -sqrt_chk(2.)*H*( (W+M)*abs_mom_k_L/N )* G_V2 [i];

    R0_A[i] = (sqrt_chk(2.)/(6.*W))*(W + M)* G_A0[i] ;
    R1_A[i] =  (sqrt_chk(2.)/(6.*W))*(W + M + 2.*Omega*W/N)* G_A1[i] ;
    R2_A[i] = (sqrt_chk(2.)/(6.*W))*(W + M + 4.*Omega*W/N)* G_A2[i] ;


    T1_V[i] = -G_V1[i]/(3.*W*X) ;
    T2_V[i] = -G_V2[i]/(3.*W*X) ;

    T1_A[i] = (2./3.)*(H* abs_mom_k_L/(N*X))*G_A1[i] ;
    T2_A[i] = (2./3.)*(H* abs_mom_k_L/(N*X))*G_A2[i] ;


    if (i==0){
      //P_{33}(1232) **** IBLOCK=0  ****
      fV_3 = sqrt_chk(6.)* R0_V[0] ;
      fV_1 = sqrt_chk(2.)* R0_V[0] ;
      fV1 = -sqrt_chk(2.)*R0_V[0] ;
      fV3 = -sqrt_chk(6.)*R0_V[0] ; 
      fV0_plus= 0. ;
      fV0_minus = 0. ;
      fV0R_plus = 0. ;
      fV0R_minus = 0. ;

      fA_3 = -sqrt_chk(6.)*R0_A[0] ;
      fA_1 = -sqrt_chk(2.)*R0_A[0] ;
      fA1 = -sqrt_chk(2.)*R0_A[0] ;
      fA3 = -sqrt_chk(6.)*R0_A[0] ; 
      fA0_plus =   2.*sqrt_chk(2.)*C0_minus[0] ;
      fA0_minus =  2.*sqrt_chk(2.)*C0_minus[0] ;
      fA0R_plus =  2.*sqrt_chk(2.)*C0_plus[0] ;
      fA0R_minus = 2.*sqrt_chk(2.)*C0_plus[0] ;

    }else if (i==1){
      //p_{11}(1440) **** IBLOCK=1  **** 
      fV_3 = 0.;
      fV_1 =  -(5./6.)*sqrt_chk(3.)*L*L*R2_V[1] ;
      fV1 =   -(5./6.)*sqrt_chk(3.)*L*L*R2_V[1] ;
      fV3 = 0.;
      fV0_plus =  -sqrt_chk(3./4.)*L*L*S2_KLM_minus[1] ;
      fV0_minus = -sqrt_chk(3./4.)*L*L*S2_KLM_minus[1] ; 
      fV0R_plus = -sqrt_chk(3./4.)*L*L*S2_KLM_plus[1]  ;
      fV0R_minus= -sqrt_chk(3./4.)*L*L*S2_KLM_plus[1]  ;

      fA_3 = 0.;
      fA_1 =   (5./6.)*sqrt_chk(3.)*L*L*R2_A[1] ;
      fA1 =   -(5./6.)*sqrt_chk(3.)*L*L*R2_A[1] ;
      fA3 = 0.;
      fA0_plus =   (5./6.)*sqrt_chk(3.)*L*( L*C2_minus[1] - 2.*B2_minus[1] );
      fA0_minus = - (5./6.)*sqrt_chk(3.)*L*( L*C2_minus[1] - 2.*B2_minus[1] ); 
      fA0R_plus =  (5./6.)*sqrt_chk(3.)*L*( L*C2_plus[1] - 2.*B2_plus[1] );
      fA0R_minus= - (5./6.)*sqrt_chk(3.)*L*( L*C2_plus[1] - 2.*B2_plus[1] );

    }else if (i==2){
      //D_{13}(1520)   **** IBLOCK=2  **** 
      fV_3 = 2.*sqrt_chk(9./2.)* T1_V[2] ;
      fV_1 =    sqrt_chk(6.)* T1_V[2] - (4./sqrt_chk(3.))*L*R1_V[2] ;
      fV1  =    sqrt_chk(6.)* T1_V[2] - (4./sqrt_chk(3.))*L*R1_V[2] ;
      fV3  = 2.*sqrt_chk(9./2.)* T1_V[2] ;
      fV0_plus = -2.*sqrt_chk(3.)*L*S1_KLM_minus[2]  ;
      fV0_minus= -2.*sqrt_chk(3.)*L*S1_KLM_minus[2]  ;
      fV0R_plus= -2.*sqrt_chk(3.)*L*S1_KLM_plus[2]   ;
      fV0R_minus=-2.*sqrt_chk(3.)*L*S1_KLM_plus[2]   ;

      fA_3 = 2.*sqrt_chk(9./2.)* -T1_A[2] ;
      fA_1 =    sqrt_chk(6.)* -T1_A[2] - (4./sqrt_chk(3.))*L*-R1_A[2] ;
      fA1  =    sqrt_chk(6.)* T1_A[2] - (4./sqrt_chk(3.))*L*R1_A[2] ;
      fA3  = 2.*sqrt_chk(9./2.)* T1_A[2] ;
      fA0_plus =  (4./sqrt_chk(3.))*L*C1_minus[2] ;
      fA0_minus= - (4./sqrt_chk(3.))*L*C1_minus[2] ;
      fA0R_plus=  (4./sqrt_chk(3.))*L*C1_plus[2] ;
      fA0R_minus=-(4./sqrt_chk(3.))*L*C1_plus[2] ;

    }else if (i==3){
      //S_{11}(1535)   **** IBLOCK=3  **** 
      fV_3 =  0.;
      fV_1 =  2.*sqrt_chk(3.)*T1_V[3] + (4./sqrt_chk(6.))*L*R1_V[3] ;
      fV1  = -2.*sqrt_chk(3.)*T1_V[3] - (4./sqrt_chk(6.))*L*R1_V[3] ;
      fV3  =  0.;
      fV0_plus  =  sqrt_chk(6.)*L*S1_KLM_minus[3] ;
      fV0_minus = -sqrt_chk(6.)*L*S1_KLM_minus[3] ;
      fV0R_plus =  sqrt_chk(6.)*L*S1_KLM_plus[3]  ;
      fV0R_minus= -sqrt_chk(6.)*L*S1_KLM_plus[3]  ;

      fA_3 =  0.;
      fA_1 =  2.*sqrt_chk(3.)*-T1_A[3] + (4./sqrt_chk(6.))*L*-R1_A[3] ;
      fA1  = -2.*sqrt_chk(3.)*T1_A[3] - (4./sqrt_chk(6.))*L*R1_A[3] ;
      fA3  =  0.;
      fA0_plus  = - 2.*sqrt_chk(2./3.)*( L*C1_minus[3] - 3.*B1_minus[3] ) ;
      fA0_minus = - 2.*sqrt_chk(2./3.)*( L*C1_minus[3] - 3.*B1_minus[3] ) ;
      fA0R_plus = - 2.*sqrt_chk(2./3.)*( L*C1_plus[3] - 3.*B1_plus[3] ) ;
      fA0R_minus= - 2.*sqrt_chk(2./3.)*( L*C1_plus[3] - 3.*B1_plus[3] ) ;

    }else if (i==4){
      //S_{31}(1620)    **** IBLOCK=4  **** 
      fV_3 =  0.;
      fV_1 = -sqrt_chk(3.)* T1_V[4] + sqrt_chk(1./6.)*L*R1_V[4] ;
      fV1  =  sqrt_chk(3.)* T1_V[4] - sqrt_chk(1./6.)*L*R1_V[4] ;
      fV3  =  0.;
      fV0_plus  = -sqrt_chk(3./2.)*L*S1_KLM_minus[4] ;
      fV0_minus =  sqrt_chk(3./2.)*L*S1_KLM_minus[4] ;
      fV0R_plus = -sqrt_chk(3./2.)*L*S1_KLM_plus[4]  ;
      fV0R_minus=  sqrt_chk(3./2.)*L*S1_KLM_plus[4]  ;

      fA_3 =  0.;
      fA_1 = -sqrt_chk(3.)* -T1_A[4] + sqrt_chk(1./6.)*L*-R1_A[4] ;
      fA1  =  sqrt_chk(3.)* T1_A[4] - sqrt_chk(1./6.)*L*R1_A[4] ;
      fA3  =  0.;
      fA0_plus  = - sqrt_chk(1./6.)*( L*C1_minus[4] - 3.*B1_minus[4] ) ;
      fA0_minus = - sqrt_chk(1./6.)*( L*C1_minus[4] - 3.*B1_minus[4] ) ;
      fA0R_plus = - sqrt_chk(1./6.)*( L*C1_plus[4] - 3.*B1_plus[4]) ;
      fA0R_minus= - sqrt_chk(1./6.)*( L*C1_plus[4] - 3.*B1_plus[4] ) ;

    }else if (i==5){
      //S_{11}(1650)  **** IBLOCK=5  **** 
      fV_3 =  0.;
      fV_1 =  sqrt_chk(1./6.)*L*R1_V[5] ;
      fV1  = -sqrt_chk(1./6.)*L*R1_V[5] ;
      fV3  =  0.;
      fV0_plus  = 0. ;
      fV0_minus = 0. ;
      fV0R_plus = 0. ;
      fV0R_minus= 0. ;

      fA_3 =  0.;
      fA_1 =  sqrt_chk(1./6.)*L*-R1_A[5] ;
      fA1  = -sqrt_chk(1./6.)*L*R1_A[5] ;
      fA3  =  0.;
      fA0_plus  = sqrt_chk(2./3.)*( L*C1_minus[5]- 3.*B1_minus[5] ) ;
      fA0_minus = sqrt_chk(2./3.)*( L*C1_minus[5]- 3.*B1_minus[5] ) ;
      fA0R_plus = sqrt_chk(2./3.)*( L*C1_plus[5] - 3.*B1_plus[5] ) ;
      fA0R_minus= sqrt_chk(2./3.)*( L*C1_plus[5] - 3.*B1_plus[5] ) ;
    }else if (i==6){
      //D_{15}(1675)  **** IBLOCK=6  **** 
      fV_3 = -sqrt_chk(3./5.)*L* R1_V[6] ;
      fV_1 = -sqrt_chk(3./10.)*L* R1_V[6] ;
      fV1  =  sqrt_chk(3./10.)*L* R1_V[6] ; 
      fV3  =  sqrt_chk(3./5.)*L* R1_V[6] ; 
      fV0_plus  = 0. ;
      fV0_minus = 0. ;
      fV0R_plus = 0. ;
      fV0R_minus= 0. ;

      fA_3 =  sqrt_chk(3./5.)*L* R1_A[6] ;
      fA_1 =  sqrt_chk(3./10.)*L* R1_A[6] ;
      fA1  =  sqrt_chk(3./10.)*L* R1_A[6] ; 
      fA3  =  sqrt_chk(3./5.)*L* R1_A[6] ; 
      fA0_plus  = -sqrt_chk(6./5.)*L*C1_minus[6] ;
      fA0_minus = -sqrt_chk(6./5.)*L*C1_minus[6] ;
      fA0R_plus = -sqrt_chk(6./5.)*L*C1_plus[6] ;
      fA0R_minus= -sqrt_chk(6./5.)*L*C1_plus[6] ;

    }else if (i==7){
      //F_{15(1680)  **** IBLOCK=7  **** 
      fV_3 =  -sqrt_chk(18./5.)*L*T2_V[7] ;
      fV_1 =  -sqrt_chk(9./5.)*L*T2_V[7] + sqrt_chk(5./2.)*L*L* R2_V[7] ;
      fV1  =  -sqrt_chk(9./5.)*L*T2_V[7]  + sqrt_chk(5./2.)*L*L* R2_V[7] ; 
      fV3  =  -sqrt_chk(18./5.)*L*T2_V[7] ;
      fV0_plus  = sqrt_chk(9./10.)*L*L*S2_KLM_minus[7]  ;
      fV0_minus = sqrt_chk(9./10.)*L*L*S2_KLM_minus[7]  ;
      fV0R_plus = sqrt_chk(9./10.)*L*L*S2_KLM_plus[7]   ;
      fV0R_minus= sqrt_chk(9./10.)*L*L*S2_KLM_plus[7]   ;

      fA_3 =   sqrt_chk(18./5.)*L*T2_A[7] ;
      fA_1 =   sqrt_chk(9./5.)*L*T2_A[7]  - sqrt_chk(5./2.)*L*L* R2_A[7] ;
      fA1  =  -sqrt_chk(9./5.)*L*T2_A[7]  + sqrt_chk(5./2.)*L*L* R2_A[7] ; 
      fA3  =  -sqrt_chk(18./5.)*L*T2_A[7] ;
      fA0_plus  = - sqrt_chk(5./2.)*L*L*C2_minus[7] ;
      fA0_minus =   sqrt_chk(5./2.)*L*L*C2_minus[7] ;
      fA0R_plus = - sqrt_chk(5./2.)*L*L*C2_plus[7] ;
      fA0R_minus=   sqrt_chk(5./2.)*L*L*C2_plus[7] ;

    }else if (i==8){
      //D_{13}(1700) **** IBLOCK=8  **** 
      fV_3 = sqrt_chk(9./10.)*L* R1_V[8] ;
      fV_1 = sqrt_chk(1./30.)*L* R1_V[8] ;
      fV1 =  sqrt_chk(1./30.)*L* R1_V[8] ;
      fV3 =  sqrt_chk(9./10.)*L* R1_V[8] ;
      fV0_plus  = 0. ;
      fV0_minus = 0. ;
      fV0R_plus = 0. ;
      fV0R_minus= 0. ;

      fA_3 = -sqrt_chk(9./10.)*L* R1_A[8] ;
      fA_1 = -sqrt_chk(1./30.)*L* R1_A[8] ;
      fA1 =   sqrt_chk(1./30.)*L* R1_A[8] ;
      fA3 =   sqrt_chk(9./10.)*L* R1_A[8] ;
      fA0_plus  =  sqrt_chk(2./15.)*L*C1_minus[8] ;
      fA0_minus = -sqrt_chk(2./15.)*L*C1_minus[8] ;
      fA0R_plus =  sqrt_chk(2./15.)*L*C1_plus[8] ;
      fA0R_minus= -sqrt_chk(2./15.)*L*C1_plus[8] ;

    }else if (i==9){
      //D_{33}(1700)  **** IBLOCK=9  **** 
      fV_3 = -sqrt_chk(9./2.)* T1_V[9] ;
      fV_1 = -sqrt_chk(3./2.)* T1_V[9] - sqrt_chk(1./3.)*L* R1_V[9] ;
      fV1  = -sqrt_chk(3./2.)* T1_V [9] - sqrt_chk(1./3.)*L* R1_V[9] ;
      fV3  = -sqrt_chk(9./2.)* T1_V[9] ;
      fV0_plus  = sqrt_chk(3.)*L* S1_KLM_minus[9]  ;
      fV0_minus = sqrt_chk(3.)*L* S1_KLM_minus[9]  ;
      fV0R_plus = sqrt_chk(3.)*L* S1_KLM_plus[9]   ;
      fV0R_minus= sqrt_chk(3.)*L* S1_KLM_plus[9]   ;

      fA_3 =  sqrt_chk(9./2.)* T1_A[9] ;
      fA_1 =  sqrt_chk(3./2.)* T1_A[9] + sqrt_chk(1./3.)*L* R1_A[9] ;
      fA1  = -sqrt_chk(3./2.)* T1_A[9] - sqrt_chk(1./3.)*L* R1_A[9] ;
      fA3  = -sqrt_chk(9./2.)* T1_A[9] ;
      fA0_plus  =   sqrt_chk(1./3.)*L*C1_minus[9] ;
      fA0_minus = - sqrt_chk(1./3.)*L*C1_minus[9] ;
      fA0R_plus =   sqrt_chk(1./3.)*L*C1_plus[9] ;
      fA0R_minus= - sqrt_chk(1./3.)*L*C1_plus[9] ;
    }else if (i==10){
      //P_{11}(1710)    **** IBLOCK=10  ****
      fV_3 = 0.;
      fV_1 = sqrt_chk(2./3.)*L*L*R2_V[10] ;
      fV1 =  sqrt_chk(2./3.)*L*L*R2_V[10] ;
      fV3 = 0.;
      fV0_plus  = sqrt_chk(3./2.)*L*L*S2_KLM_minus[10]  ;
      fV0_minus = sqrt_chk(3./2.)*L*L*S2_KLM_minus[10]  ;
      fV0R_plus = sqrt_chk(3./2.)*L*L*S2_KLM_plus[10]   ;
      fV0R_minus= sqrt_chk(3./2.)*L*L*S2_KLM_plus[10]   ;

      fA_3 = 0.;
      fA_1 = -sqrt_chk(2./3.)*L*L*R2_A[10] ;
      fA1 =  -sqrt_chk(2./3.)*L*L*R2_A[10] ;
      fA3 = 0.;
      fA0_plus  = - sqrt_chk(2./3.)*L*( L*C2_minus[10] - 2.*B2_minus[10] ) ;
      fA0_minus =   sqrt_chk(2./3.)*L*( L*C2_minus[10] - 2.*B2_minus[10] ) ;
      fA0R_plus = - sqrt_chk(2./3.)*L*( L*C2_plus[10] - 2.*B2_plus[10] ) ;
      fA0R_minus=   sqrt_chk(2./3.)*L*( L*C2_plus[10] - 2.*B2_plus[10] ) ;

    }else if (i==11){
      //P_{13}(1720)   **** IBLOCK=11  ****
      fV_3 =  sqrt_chk(9./10.)*L*T2_V[11] ;
      fV_1 = -sqrt_chk(27./10.)*L*T2_V[11] - sqrt_chk(5./3.)*L*L* R2_V[11]  ;
      fV1  =  sqrt_chk(27./10.)*L*T2_V[11]  + sqrt_chk(5./3.)*L*L* R2_V[11]  ; 
      fV3  = -sqrt_chk(9./10.)*L*T2_V[11] ;
      fV0_plus  = -sqrt_chk(3./5.)*L*L*S2_KLM_minus[11]  ;
      fV0_minus =  sqrt_chk(3./5.)*L*L*S2_KLM_minus[11]   ; 
      fV0R_plus = -sqrt_chk(3./5.)*L*L*S2_KLM_plus[11]   ;
      fV0R_minus=  sqrt_chk(3./5.)*L*L*S2_KLM_plus[11]   ;

      fA_3 = -sqrt_chk(9./10.)*L*T2_A[11] ;
      fA_1 =  sqrt_chk(27./10.)*L*T2_A[11] - sqrt_chk(5./3.)*L*L* -R2_A[11]  ;
      fA1  =  sqrt_chk(27./10.)*L*T2_A[11]  + sqrt_chk(5./3.)*L*L* R2_A[11]  ; 
      fA3  = -sqrt_chk(9./10.)*L*T2_A[11] ;
      fA0_plus  =  sqrt_chk(5./3.)*L*(L*C2_minus[11] - 5.*B2_minus[11] )  ;
      fA0_minus =  sqrt_chk(5./3.)*L*(L*C2_minus[11] - 5.*B2_minus[11] )  ; 
      fA0R_plus =  sqrt_chk(5./3.)*L*(L*C2_plus [11]- 5.*B2_plus [11]) ;
      fA0R_minus=  sqrt_chk(5./3.)*L*(L*C2_plus[11] - 5.*B2_plus[11] ) ;

    }else if (i==12){
      //F_{35}(1905)    **** IBLOCK=12  ****
      fV_3 = -sqrt_chk(18./35.)*L*L* R2_V[12] ;
      fV_1 = -sqrt_chk(1./35.)*L*L* R2_V[12] ;
      fV1  = -sqrt_chk(1./35.)*L*L* R2_V[12] ;
      fV3  = -sqrt_chk(18./35.)*L*L* R2_V [12];
      fV0_plus  = 0. ;
      fV0_minus = 0. ;
      fV0R_plus = 0.;
      fV0R_minus= 0.;

      fA_3 =  sqrt_chk(18./35.)*L*L*R2_A[12] ;
      fA_1 =  sqrt_chk(1./35.)*L*L* R2_A[12] ;
      fA1  = -sqrt_chk(1./35.)*L*L* R2_A[12] ;
      fA3  = -sqrt_chk(18./35.)*L*L* R2_A [12];
      fA0_plus  = -sqrt_chk(4./35.)*L*L*C2_minus[12] ;
      fA0_minus =  sqrt_chk(4./35.)*L*L*C2_minus[12] ;
      fA0R_plus = -sqrt_chk(4./35.)*L*L*C2_plus[12] ;
      fA0R_minus=  sqrt_chk(4./35.)*L*L*C2_plus [12];

    }else if (i==13){
      //P_{31}(1910)    **** IBLOCK=13  ****
      fV_3 = 0.;
      fV_1 = sqrt_chk(1./15.)*L*L*R2_V[13] ;
      fV1  = sqrt_chk(1./15.)*L*L*R2_V[13]  ;
      fV3  = 0.;
      fV0_plus  = 0.;
      fV0_minus = 0.;
      fV0R_plus = 0.;
      fV0R_minus= 0.;

      fA_3 = 0.;
      fA_1 = -sqrt_chk(1./15.)*L*L*R2_A[13] ;
      fA1  =  sqrt_chk(1./15.)*L*L*R2_A[13]  ;
      fA3  = 0.;
      fA0_plus  =  sqrt_chk(4./15.)*L*( L*C2_minus[13]  - 5.*B2_minus[13]  );
      fA0_minus = -sqrt_chk(4./15.)*L*( L*C2_minus[13]  - 5.*B2_minus[13]  );
      fA0R_plus =  sqrt_chk(4./15.)*L*( L*C2_plus [13]  - 5.*B2_plus[13]  );
      fA0R_minus= -sqrt_chk(4./15.)*L*( L*C2_plus[13]   - 5.*B2_plus[13]  );

    }else if (i==14){
      //P_{33}(1920)   **** IBLOCK=14  ****
      fV_3 =  sqrt_chk(1./5. )*L*L* R2_V[14] ;
      fV_1 = -sqrt_chk(1./15.)*L*L* R2_V[14] ;
      fV1  =  sqrt_chk(1./15.)*L*L* R2_V[14] ;
      fV3  = -sqrt_chk(1./5. )*L*L* R2_V[14] ;
      fV0_plus  = 0.;
      fV0_minus = 0.;
      fV0R_plus = 0.;
      fV0R_minus= 0.;

      fA_3 = -sqrt_chk(1./5. )*L*L*R2_A[14] ;
      fA_1 =  sqrt_chk(1./15.)*L*L*R2_A[14] ;
      fA1  =  sqrt_chk(1./15.)*L*L* R2_A[14] ;
      fA3  = -sqrt_chk(1./5. )*L*L* R2_A[14] ;
      fA0_plus  = -sqrt_chk(4./15.)*L*( L*C2_minus[14] - 5.*B2_minus[14] );
      fA0_minus = -sqrt_chk(4./15.)*L*( L*C2_minus [14]- 5.*B2_minus[14] );
      fA0R_plus = -sqrt_chk(4./15.)*L*( L*C2_plus [14] - 5.*B2_plus[14] );
      fA0R_minus= -sqrt_chk(4./15.)*L*( L*C2_plus[14]  - 5.*B2_plus[14] );

    }else if (i==15){
      //F_{37}(1950)   **** IBLOCK=15  ****
      fV_3 =  sqrt_chk(2./7.)*L*L* R2_V[15] ;
      fV_1 =  sqrt_chk(6./35.)*L*L* R2_V[15] ;
      fV1  = -sqrt_chk(6./35.)*L*L* R2_V [15];
      fV3  = -sqrt_chk(2./7.)*L*L* R2_V [15];
      fV0_plus  = 0. ;
      fV0_minus = 0. ;
      fV0R_plus = 0. ;
      fV0R_minus= 0. ;

      fA_3 = -sqrt_chk(2./7. )*L*L* R2_A[15] ;
      fA_1 = -sqrt_chk(6./35.)*L*L* R2_A[15] ;
      fA1  = -sqrt_chk(6./35.)*L*L* R2_A [15];
      fA3  = -sqrt_chk(2./7. )*L*L* R2_A [15];
      fA0_plus  = 2.*sqrt_chk(6./35.)*L*L*C2_minus[15] ;
      fA0_minus = 2.*sqrt_chk(6./35.)*L*L*C2_minus[15] ;
      fA0R_plus = 2.*sqrt_chk(6./35.)*L*L*C2_plus[15] ;
      fA0R_minus= 2.*sqrt_chk(6./35.)*L*L*C2_plus[15] ;
    }else if (i==16){
      //P_{33}(1600)   **** IBLOCK=16  ****
      fV_3 = -sqrt_chk(1./2.)*L*L*R2_V[16] ;
      fV_1 = -sqrt_chk(1./6.)*L*L*R2_V[16] ;
      fV1  = sqrt_chk(1./6.)*L*L*R2_V[16] ;
      fV3  = sqrt_chk(1./2.)*L*L*R2_V[16] ;
      fV0_plus  = 0. ;
      fV0_minus = fV0_plus ;
      fV0R_plus = 0.;
      fV0R_minus= fV0R_plus ;

      fA_3 = sqrt_chk(1./2.)*L*L*R2_A[16] ;
      fA_1 = sqrt_chk(1./6.)*L*L*R2_A[16] ;
      fA1  = sqrt_chk(1./6.)*L*L*R2_A[16] ;
      fA3  = sqrt_chk(1./2.)*L*L*R2_A[16] ;
      fA0_plus  = -sqrt_chk(2./3.)*L*(L*C2_minus[16]-2.*B2_minus[16]) ;
      fA0_minus = fA0_plus ;
      fA0R_plus = -sqrt_chk(2./3.)*L*(L*C2_plus[16]-2.*B2_plus[16]) ;
      fA0R_minus= fA0R_plus ;

    }
    // Vector helicity amplitudes
    HV_realM11[i]  = sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_real[i]*fV3*Jsgn[i]  ;
    HV_realM_11[i] = sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_real[i]*fV3;
    HV_realM1_1[i] = sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_real[i]*fV1*Jsgn[i] ;
    HV_realM_1_1[i]= sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_real[i]*fV1 ;

    HV_realP11[i]  =  -sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_real[i]*fV_1 ;
    HV_realP_11[i] =   sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_real[i]*fV_1*Jsgn[i];
    HV_realP1_1[i] =   sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_real[i]*fV_3;
    HV_realP_1_1[i]=  -sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_real[i]*fV_3*Jsgn[i] ;

    HV0_realM11[i]  =  sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_real[i]*fV0_minus*Jsgn[i] ;
    HV0_realM_11[i] =  sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_real[i]*fV0_minus;
    HV0_realM1_1[i] = -sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_real[i]*fV0_plus*Jsgn[i] ;
    HV0_realM_1_1[i]=  sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_real[i]*fV0_plus ;

    HV0_realP11[i]  =  sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_real[i]*fV0R_minus*Jsgn[i] ;
    HV0_realP_11[i] =  sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_real[i]*fV0R_minus;
    HV0_realP1_1[i] = -sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_real[i]*fV0R_plus*Jsgn[i]  ;
    HV0_realP_1_1[i]=  sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_real[i]*fV0R_plus ;

    HV_ImM11[i]  =    sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_Im[i]*fV3*Jsgn[i] ;
    HV_ImM_11[i] =    sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_Im[i]*fV3;
    HV_ImM1_1[i] =    sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_Im[i]*fV1*Jsgn[i] ;
    HV_ImM_1_1[i]=    sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_Im[i]*fV1 ;

    HV_ImP11[i]  =  -sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_Im[i]*fV_1 ;
    HV_ImP_11[i] =   sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_Im[i]*fV_1*Jsgn[i];
    HV_ImP1_1[i] =   sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_Im[i]*fV_3 ;
    HV_ImP_1_1[i]=  -sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_Im[i]*fV_3*Jsgn[i] ;

    HV0_ImM11[i]  =  sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_Im[i]*fV0_minus*Jsgn[i] ;
    HV0_ImM_11[i] =  sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_Im[i]*fV0_minus;
    HV0_ImM1_1[i] = -sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_Im[i]*fV0_plus*Jsgn[i]  ;
    HV0_ImM_1_1[i]=  sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_Im[i]*fV0_plus ;

    HV0_ImP11[i]  =  sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_Im[i]*fV0R_minus*Jsgn[i] ;
    HV0_ImP_11[i] =  sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_Im[i]*fV0R_minus;
    HV0_ImP1_1[i] = -sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_Im[i]*fV0R_plus*Jsgn[i]  ;
    HV0_ImP_1_1[i]=  sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Vf_BW_Im[i]*fV0R_plus ; 

    //Axial helicity amplitudes
    HA_realM11[i]  = sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_real[i]*fA3*Jsgn[i]  ;
    HA_realM_11[i] = sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_real[i]*fA3;
    HA_realM1_1[i] = sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_real[i]*fA1*Jsgn[i] ;
    HA_realM_1_1[i]= sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_real[i]*fA1 ;

    HA_realP11[i]  =  -sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_real[i]*fA_1 ;
    HA_realP_11[i] =   sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_real[i]*fA_1*Jsgn[i];
    HA_realP1_1[i] =   sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_real[i]*fA_3;
    HA_realP_1_1[i]=  -sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_real[i]*fA_3*Jsgn[i] ;

    HA0_realM11[i]  =  sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_real[i]*fA0_minus*Jsgn[i] ;
    HA0_realM_11[i] =  sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_real[i]*fA0_minus;
    HA0_realM1_1[i] = -sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_real[i]*fA0_plus*Jsgn[i] ;
    HA0_realM_1_1[i]=  sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_real[i]*fA0_plus ;

    HA0_realP11[i]  =  sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_real[i]*fA0R_minus*Jsgn[i] ;
    HA0_realP_11[i] =  sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_real[i]*fA0R_minus;
    HA0_realP1_1[i] = -sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_real[i]*fA0R_plus*Jsgn[i]  ;
    HA0_realP_1_1[i]=  sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_real[i]*fA0R_plus ;

    HA_ImM11[i]  =    sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_Im[i]*fA3*Jsgn[i] ;
    HA_ImM_11[i] =    sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_Im[i]*fA3;
    HA_ImM1_1[i] =    sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_Im[i]*fA1*Jsgn[i] ;
    HA_ImM_1_1[i]=    sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_Im[i]*fA1 ;

    HA_ImP11[i]  =  -sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_Im[i]*fA_1 ;
    HA_ImP_11[i] =   sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_Im[i]*fA_1*Jsgn[i];
    HA_ImP1_1[i] =   sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_Im[i]*fA_3 ;
    HA_ImP_1_1[i]=  -sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_Im[i]*fA_3*Jsgn[i] ;

    HA0_ImM11[i]  =  sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_Im[i]*fA0_minus*Jsgn[i] ;
    HA0_ImM_11[i] =  sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_Im[i]*fA0_minus;
    HA0_ImM1_1[i] = -sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_Im[i]*fA0_plus*Jsgn[i]  ;
    HA0_ImM_1_1[i]=  sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_Im[i]*fA0_plus ;

    HA0_ImP11[i]  =  sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_Im[i]*fA0R_minus*Jsgn[i] ;
    HA0_ImP_11[i] =  sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_Im[i]*fA0R_minus;
    HA0_ImP1_1[i] = -sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_Im[i]*fA0R_plus*Jsgn[i]  ;
    HA0_ImP_1_1[i]=  sqrt_chk(2.)*JP[i]*Dsgn[i]*kapa[i]*Af_BW_Im[i]*fA0R_plus ; 
    }

    double pt_1 = 1. ;
    double pt_2 = 3.*t ;
    double pt_3 = (15./2.)*t*t -3./2. ;
    double pt_4 = (35./2.)*t*t*t - (15./2.)*t ;

    double d111 = sqrt_chk((1 + t)/2.);
    double d11_1 =  -sqrt_chk((1- t)/2.);
    double d131 = 0.;
    double d13_1 = 0.;

    double d311 = (1./2.)*sqrt_chk((1 + t)/2.)*(pt_2 - pt_1) ;
    double d31_1 =(-1./2.)*sqrt_chk((1- t)/2.)*(pt_2 + pt_1) ;
    double d331 =(-1./2.)*sqrt_chk((1- t)/2.)*(sqrt_chk(1./3.)*pt_2 + sqrt_chk(3.)*pt_1) ;
    double d33_1 =(1./2.)*sqrt_chk((1 + t)/2.)*(-sqrt_chk(1./3.)*pt_2 + sqrt_chk(3.)*pt_1) ;

    double d511 = (1./3.)* sqrt_chk((1+t)/2.)*(pt_3 - pt_2) ;
    double d51_1 = (-1./3.)* sqrt_chk((1-t)/2.)*(pt_3 + pt_2);
    double d531 = (-1./3.)* sqrt_chk((1-t)/2.)*( sqrt_chk(1./2.)*pt_3 + sqrt_chk(2.)*pt_2); 
    double d53_1 = (1./3.)* sqrt_chk((1+t)/2.)*(-sqrt_chk(1./2.)*pt_3 + sqrt_chk(2.)*pt_2); 


    double d711 = (1./4.)* sqrt_chk((1+t)/2.)*(pt_4 - pt_3);
    double d71_1 = (-1./4.)* sqrt_chk((1-t)/2.)*(pt_4 + pt_3);
    double d731 = (-1./4.)* sqrt_chk((1-t)/2.)*(sqrt_chk(3./5.)*pt_4 + sqrt_chk(5./3.)*pt_3);
    double d73_1 = (1./4.)* sqrt_chk((1+t)/2.)*(-sqrt_chk(3./5.)*pt_4 + sqrt_chk(5./3.)*pt_3) ;

    //sum of vextor HAmplitudes
    double Vsum3Re_M11  = HV_realM11[0]*  d331 + HV_realM11[4]* d131 + HV_realM11[9]*  d331 + HV_realM11[12]* d531 + HV_realM11[13]* d131 + HV_realM11[14]* d331 + HV_realM11[15]*d731 + HV_realM11[16]* d331;
    double Vsum3Re_M_11 = HV_realM_11[0]* d33_1+ HV_realM_11[4]*d13_1+ HV_realM_11[9]*d33_1 + HV_realM_11[12]*d53_1 +HV_realM_11[13]*d13_1+ HV_realM_11[14]*d33_1 + HV_realM_11[15]*d73_1 + HV_realM_11[16]*d33_1;
    double Vsum3Re_M1_1 = HV_realM1_1[0]* d311 + HV_realM1_1[4]*d111 + HV_realM1_1[9]* d311 + HV_realM1_1[12]*d511 + HV_realM1_1[13]*d111 + HV_realM1_1[14]* d311 + HV_realM1_1[15]*d711 + HV_realM1_1[16]* d311;
    double Vsum3Re_M_1_1= HV_realM_1_1[0]*d31_1+ HV_realM_1_1[4]*d11_1+HV_realM_1_1[9]*d31_1+ HV_realM_1_1[12]*d51_1+HV_realM_1_1[13]*d11_1+ HV_realM_1_1[14]*d31_1+HV_realM_1_1[15]*d71_1+ HV_realM_1_1[16]*d31_1;

    double Vsum3Re_P11  = HV_realP11[0]*d31_1 + HV_realP11[4]*d11_1 + HV_realP11[9]*d31_1 + HV_realP11[12]*d51_1 + HV_realP11[13]*d11_1 + HV_realP11[14]*d31_1 + HV_realP11[15]*d71_1 + HV_realP11[16]*d31_1;
    double Vsum3Re_P_11 = HV_realP_11[0]*d311 + HV_realP_11[4]*d111 + HV_realP_11[9]*d311 + HV_realP_11[12]*d511 + HV_realP_11[13]*d111 + HV_realP_11[14]*d311 + HV_realP_11[15]*d711 + HV_realP_11[16]*d311;
    double Vsum3Re_P1_1 = HV_realP1_1[0]*d33_1+ HV_realP1_1[4]*d13_1+ HV_realP1_1[9]*d33_1+ HV_realP1_1[12]*d53_1+ HV_realP1_1[13]*d13_1+ HV_realP1_1[14]*d33_1+ HV_realP1_1[15]*d73_1+ HV_realP1_1[16]*d33_1;
    double Vsum3Re_P_1_1= HV_realP_1_1[0]*d331+ HV_realP_1_1[4]*d131+ HV_realP_1_1[9]*d331+ HV_realP_1_1[12]*d531+ HV_realP_1_1[13]*d131+ HV_realP_1_1[14]*d331+ HV_realP_1_1[15]*d731+ HV_realP_1_1[16]*d331;

    double Vsum3Re_0M11  = HV0_realM11[0]*d311  + HV0_realM11[4]* d111 + HV0_realM11[9]* d311 + HV0_realM11[12]* d511 + HV0_realM11[13]* d111 + HV0_realM11[14]* d311 + HV0_realM11[15]*d711 + HV0_realM11[16]* d311;
    double Vsum3Re_0M_11 = HV0_realM_11[0]*d31_1+ HV0_realM_11[4]*d11_1+ HV0_realM_11[9]*d31_1+ HV0_realM_11[12]*d51_1+ HV0_realM_11[13]*d11_1+ HV0_realM_11[14]*d31_1+ HV0_realM_11[15]*d71_1+ HV0_realM_11[16]*d31_1;
    double Vsum3Re_0M1_1 = HV0_realM1_1[0]*d31_1+ HV0_realM1_1[4]*d11_1+ HV0_realM1_1[9]*d31_1+ HV0_realM1_1[12]*d51_1+ HV0_realM1_1[13]*d11_1+ HV0_realM1_1[14]*d31_1+ HV0_realM1_1[15]*d71_1+ HV0_realM1_1[16]*d31_1;
    double Vsum3Re_0M_1_1= HV0_realM_1_1[0]*d311+ HV0_realM_1_1[4]*d111+ HV0_realM_1_1[9]*d311+ HV0_realM_1_1[12]*d511+ HV0_realM_1_1[13]*d111+ HV0_realM_1_1[14]*d311+ HV0_realM_1_1[15]*d711+ HV0_realM_1_1[16]*d311;

    double Vsum3Re_0P11  = HV0_realP11[0]*d311  + HV0_realP11[4]* d111 + HV0_realP11[9]* d311 + HV0_realP11[12]* d511 + HV0_realP11[13]* d111 + HV0_realP11[14]* d311 + HV0_realP11[15]*d711 + HV0_realP11[16]* d311;
    double Vsum3Re_0P_11 = HV0_realP_11[0]*d31_1+ HV0_realP_11[4]*d11_1+ HV0_realP_11[9]*d31_1+ HV0_realP_11[12]*d51_1+ HV0_realP_11[13]*d11_1+ HV0_realP_11[14]*d31_1+ HV0_realP_11[15]*d71_1+HV0_realP_11[16]*d31_1;
    double Vsum3Re_0P1_1 = HV0_realP1_1[0]*d31_1+ HV0_realP1_1[4]*d11_1+ HV0_realP1_1[9]*d31_1+ HV0_realP1_1[12]*d51_1+ HV0_realP1_1[13]*d11_1+ HV0_realP1_1[14]*d31_1+ HV0_realP1_1[15]*d71_1+ HV0_realP1_1[16]*d31_1;
    double Vsum3Re_0P_1_1= HV0_realP_1_1[0]*d311+ HV0_realP_1_1[4]*d111+ HV0_realP_1_1[9]*d311+ HV0_realP_1_1[12]*d511+ HV0_realP_1_1[13]*d111+ HV0_realP_1_1[14]*d311+ HV0_realP_1_1[15]*d711+ HV0_realP_1_1[16]*d311;

    double Vsum3Im_M11  = HV_ImM11[0]*  d331 + HV_ImM11[4]* d131 + HV_ImM11[9]*  d331 + HV_ImM11[12]* d531 + HV_ImM11[13]* d131 + HV_ImM11[14]* d331 + HV_ImM11[15]*d731 + HV_ImM11[16]* d331 ;
    double Vsum3Im_M_11 = HV_ImM_11[0]* d33_1+ HV_ImM_11[4]*d13_1+ HV_ImM_11[9]*d33_1 + HV_ImM_11[12]*d53_1 +HV_ImM_11[13]*d13_1+ HV_ImM_11[14]*d33_1 + HV_ImM_11[15]*d73_1+ HV_ImM_11[16]*d33_1;
    double Vsum3Im_M1_1 = HV_ImM1_1[0]* d311 + HV_ImM1_1[4]*d111 + HV_ImM1_1[9]* d311 + HV_ImM1_1[12]*d511 + HV_ImM1_1[13]*d111 + HV_ImM1_1[14]* d311 + HV_ImM1_1[15]*d711+ HV_ImM1_1[16]* d311;
    double Vsum3Im_M_1_1= HV_ImM_1_1[0]*d31_1+ HV_ImM_1_1[4]*d11_1+HV_ImM_1_1[9]*d31_1+ HV_ImM_1_1[12]*d51_1+HV_ImM_1_1[13]*d11_1+ HV_ImM_1_1[14]*d31_1+HV_ImM_1_1[15]*d71_1+ HV_ImM_1_1[16]*d31_1;

    double Vsum3Im_P11  = HV_ImP11[0]*d31_1 + HV_ImP11[4]*d11_1 + HV_ImP11[9]*d31_1 + HV_ImP11[12]*d51_1 + HV_ImP11[13]*d11_1 + HV_ImP11[14]*d31_1 + HV_ImP11[15]*d71_1 + HV_ImP11[16]*d31_1;
    double Vsum3Im_P_11 = HV_ImP_11[0]*d311 + HV_ImP_11[4]*d111 + HV_ImP_11[9]*d311 + HV_ImP_11[12]*d511 + HV_ImP_11[13]*d111 + HV_ImP_11[14]*d311 + HV_ImP_11[15]*d711+ HV_ImP_11[16]*d311;
    double Vsum3Im_P1_1 = HV_ImP1_1[0]*d33_1+ HV_ImP1_1[4]*d13_1+ HV_ImP1_1[9]*d33_1+ HV_ImP1_1[12]*d53_1+ HV_ImP1_1[13]*d13_1+ HV_ImP1_1[14]*d33_1+ HV_ImP1_1[15]*d73_1+ HV_ImP1_1[16]*d33_1;
    double Vsum3Im_P_1_1= HV_ImP_1_1[0]*d331+ HV_ImP_1_1[4]*d131+ HV_ImP_1_1[9]*d331+ HV_ImP_1_1[12]*d531+ HV_ImP_1_1[13]*d131+ HV_ImP_1_1[14]*d331+ HV_ImP_1_1[15]*d731+ HV_ImP_1_1[16]*d331;

    double Vsum3Im_0M11  = HV0_ImM11[0]*d311  + HV0_ImM11[4]* d111 + HV0_ImM11[9]* d311 + HV0_ImM11[12]* d511 + HV0_ImM11[13]* d111 + HV0_ImM11[14]* d311 + HV0_ImM11[15]*d711 + HV0_ImM11[16]* d311;
    double Vsum3Im_0M_11 = HV0_ImM_11[0]*d31_1+ HV0_ImM_11[4]*d11_1+ HV0_ImM_11[9]*d31_1+ HV0_ImM_11[12]*d51_1+ HV0_ImM_11[13]*d11_1+ HV0_ImM_11[14]*d31_1+ HV0_ImM_11[15]*d71_1+ HV0_ImM_11[16]*d31_1;
    double Vsum3Im_0M1_1 = HV0_ImM1_1[0]*d31_1+ HV0_ImM1_1[4]*d11_1+ HV0_ImM1_1[9]*d31_1+ HV0_ImM1_1[12]*d51_1+ HV0_ImM1_1[13]*d11_1+ HV0_ImM1_1[14]*d31_1+ HV0_ImM1_1[15]*d71_1+ HV0_ImM1_1[16]*d31_1;
    double Vsum3Im_0M_1_1= HV0_ImM_1_1[0]*d311+ HV0_ImM_1_1[4]*d111+ HV0_ImM_1_1[9]*d311+ HV0_ImM_1_1[12]*d511+ HV0_ImM_1_1[13]*d111+ HV0_ImM_1_1[14]*d311+ HV0_ImM_1_1[15]*d711+ HV0_ImM_1_1[16]*d311;

    double Vsum3Im_0P11  = HV0_ImP11[0]*d311  + HV0_ImP11[4]* d111 + HV0_ImP11[9]* d311 + HV0_ImP11[12]* d511 + HV0_ImP11[13]* d111 + HV0_ImP11[14]* d311 + HV0_ImP11[15]*d711 + HV0_ImP11[16]* d311;
    double Vsum3Im_0P_11 = HV0_ImP_11[0]*d31_1+ HV0_ImP_11[4]*d11_1+ HV0_ImP_11[9]*d31_1+ HV0_ImP_11[12]*d51_1+ HV0_ImP_11[13]*d11_1+ HV0_ImP_11[14]*d31_1+ HV0_ImP_11[15]*d71_1+ HV0_ImP_11[16]*d31_1;
    double Vsum3Im_0P1_1 = HV0_ImP1_1[0]*d31_1+ HV0_ImP1_1[4]*d11_1+ HV0_ImP1_1[9]*d31_1+ HV0_ImP1_1[12]*d51_1+ HV0_ImP1_1[13]*d11_1+ HV0_ImP1_1[14]*d31_1+ HV0_ImP1_1[15]*d71_1+ HV0_ImP1_1[16]*d31_1;
    double Vsum3Im_0P_1_1= HV0_ImP_1_1[0]*d311+ HV0_ImP_1_1[4]*d111+ HV0_ImP_1_1[9]*d311+ HV0_ImP_1_1[12]*d511+ HV0_ImP_1_1[13]*d111+ HV0_ImP_1_1[14]*d311+ HV0_ImP_1_1[15]*d711+ HV0_ImP_1_1[16]*d311;

    //sum of Axial HAmplitudes
    double Asum3Re_M11  = HA_realM11[0]*  d331 + HA_realM11[4]* d131 + HA_realM11[9]*  d331 + HA_realM11[12]* d531 + HA_realM11[13]* d131 + HA_realM11[14]* d331 + HA_realM11[15]*d731 + HA_realM11[16]* d331;
    double Asum3Re_M_11 = HA_realM_11[0]* d33_1+ HA_realM_11[4]*d13_1+ HA_realM_11[9]*d33_1 + HA_realM_11[12]*d53_1 +HA_realM_11[13]*d13_1+ HA_realM_11[14]*d33_1 + HA_realM_11[15]*d73_1 + HA_realM_11[16]*d33_1;
    double Asum3Re_M1_1 = HA_realM1_1[0]* d311 + HA_realM1_1[4]*d111 + HA_realM1_1[9]* d311 + HA_realM1_1[12]*d511 + HA_realM1_1[13]*d111 + HA_realM1_1[14]* d311 + HA_realM1_1[15]*d711 + HA_realM1_1[16]* d311;
    double Asum3Re_M_1_1= HA_realM_1_1[0]*d31_1+ HA_realM_1_1[4]*d11_1+HA_realM_1_1[9]*d31_1+ HA_realM_1_1[12]*d51_1+HA_realM_1_1[13]*d11_1+ HA_realM_1_1[14]*d31_1+HA_realM_1_1[15]*d71_1+ HA_realM_1_1[16]*d31_1;

    double Asum3Re_P11  = HA_realP11[0]*d31_1 + HA_realP11[4]*d11_1 + HA_realP11[9]*d31_1 + HA_realP11[12]*d51_1 + HA_realP11[13]*d11_1 + HA_realP11[14]*d31_1 + HA_realP11[15]*d71_1 + HA_realP11[16]*d31_1;
    double Asum3Re_P_11 = HA_realP_11[0]*d311 + HA_realP_11[4]*d111 + HA_realP_11[9]*d311 + HA_realP_11[12]*d511 + HA_realP_11[13]*d111 + HA_realP_11[14]*d311 + HA_realP_11[15]*d711 + HA_realP_11[16]*d311;
    double Asum3Re_P1_1 = HA_realP1_1[0]*d33_1+ HA_realP1_1[4]*d13_1+ HA_realP1_1[9]*d33_1+ HA_realP1_1[12]*d53_1+ HA_realP1_1[13]*d13_1+ HA_realP1_1[14]*d33_1+ HA_realP1_1[15]*d73_1+ HA_realP1_1[16]*d33_1;
    double Asum3Re_P_1_1= HA_realP_1_1[0]*d331+ HA_realP_1_1[4]*d131+ HA_realP_1_1[9]*d331+ HA_realP_1_1[12]*d531+ HA_realP_1_1[13]*d131+ HA_realP_1_1[14]*d331+ HA_realP_1_1[15]*d731+ HA_realP_1_1[16]*d331;

    double Asum3Re_0M11  = HA0_realM11[0]*d311  + HA0_realM11[4]* d111 + HA0_realM11[9]* d311 + HA0_realM11[12]* d511 + HA0_realM11[13]* d111 + HA0_realM11[14]* d311 + HA0_realM11[15]*d711 + HA0_realM11[16]* d311;
    double Asum3Re_0M_11 = HA0_realM_11[0]*d31_1+ HA0_realM_11[4]*d11_1+ HA0_realM_11[9]*d31_1+ HA0_realM_11[12]*d51_1+ HA0_realM_11[13]*d11_1+ HA0_realM_11[14]*d31_1+ HA0_realM_11[15]*d71_1+ HA0_realM_11[16]*d31_1;
    double Asum3Re_0M1_1 = HA0_realM1_1[0]*d31_1+ HA0_realM1_1[4]*d11_1+ HA0_realM1_1[9]*d31_1+ HA0_realM1_1[12]*d51_1+ HA0_realM1_1[13]*d11_1+ HA0_realM1_1[14]*d31_1+ HA0_realM1_1[15]*d71_1+ HA0_realM1_1[16]*d31_1;
    double Asum3Re_0M_1_1= HA0_realM_1_1[0]*d311+ HA0_realM_1_1[4]*d111+ HA0_realM_1_1[9]*d311+ HA0_realM_1_1[12]*d511+ HA0_realM_1_1[13]*d111+ HA0_realM_1_1[14]*d311+ HA0_realM_1_1[15]*d711+ HA0_realM_1_1[16]*d311;

    double Asum3Re_0P11  = HA0_realP11[0]*d311  + HA0_realP11[4]* d111 + HA0_realP11[9]* d311 + HA0_realP11[12]* d511 + HA0_realP11[13]* d111 + HA0_realP11[14]* d311 + HA0_realP11[15]*d711 + HA0_realP11[16]* d311;
    double Asum3Re_0P_11 = HA0_realP_11[0]*d31_1+ HA0_realP_11[4]*d11_1+ HA0_realP_11[9]*d31_1+ HA0_realP_11[12]*d51_1+ HA0_realP_11[13]*d11_1+ HA0_realP_11[14]*d31_1+ HA0_realP_11[15]*d71_1+HA0_realP_11[16]*d31_1;
    double Asum3Re_0P1_1 = HA0_realP1_1[0]*d31_1+ HA0_realP1_1[4]*d11_1+ HA0_realP1_1[9]*d31_1+ HA0_realP1_1[12]*d51_1+ HA0_realP1_1[13]*d11_1+ HA0_realP1_1[14]*d31_1+ HA0_realP1_1[15]*d71_1+ HA0_realP1_1[16]*d31_1;
    double Asum3Re_0P_1_1= HA0_realP_1_1[0]*d311+ HA0_realP_1_1[4]*d111+ HA0_realP_1_1[9]*d311+ HA0_realP_1_1[12]*d511+ HA0_realP_1_1[13]*d111+ HA0_realP_1_1[14]*d311+ HA0_realP_1_1[15]*d711+ HA0_realP_1_1[16]*d311;

    double Asum3Im_M11  = HA_ImM11[0]*  d331 + HA_ImM11[4]* d131 + HA_ImM11[9]*  d331 + HA_ImM11[12]* d531 + HA_ImM11[13]* d131 + HA_ImM11[14]* d331 + HA_ImM11[15]*d731 + HA_ImM11[16]* d331 ;
    double Asum3Im_M_11 = HA_ImM_11[0]* d33_1+ HA_ImM_11[4]*d13_1+ HA_ImM_11[9]*d33_1 + HA_ImM_11[12]*d53_1 +HA_ImM_11[13]*d13_1+ HA_ImM_11[14]*d33_1 + HA_ImM_11[15]*d73_1+ HA_ImM_11[16]*d33_1;
    double Asum3Im_M1_1 = HA_ImM1_1[0]* d311 + HA_ImM1_1[4]*d111 + HA_ImM1_1[9]* d311 + HA_ImM1_1[12]*d511 + HA_ImM1_1[13]*d111 + HA_ImM1_1[14]* d311 + HA_ImM1_1[15]*d711+ HA_ImM1_1[16]* d311;
    double Asum3Im_M_1_1= HA_ImM_1_1[0]*d31_1+ HA_ImM_1_1[4]*d11_1+HA_ImM_1_1[9]*d31_1+ HA_ImM_1_1[12]*d51_1+HA_ImM_1_1[13]*d11_1+ HA_ImM_1_1[14]*d31_1+HA_ImM_1_1[15]*d71_1+ HA_ImM_1_1[16]*d31_1;

    double Asum3Im_P11  = HA_ImP11[0]*d31_1 + HA_ImP11[4]*d11_1 + HA_ImP11[9]*d31_1 + HA_ImP11[12]*d51_1 + HA_ImP11[13]*d11_1 + HA_ImP11[14]*d31_1 + HA_ImP11[15]*d71_1 + HA_ImP11[16]*d31_1;
    double Asum3Im_P_11 = HA_ImP_11[0]*d311 + HA_ImP_11[4]*d111 + HA_ImP_11[9]*d311 + HA_ImP_11[12]*d511 + HA_ImP_11[13]*d111 + HA_ImP_11[14]*d311 + HA_ImP_11[15]*d711+ HA_ImP_11[16]*d311;
    double Asum3Im_P1_1 = HA_ImP1_1[0]*d33_1+ HA_ImP1_1[4]*d13_1+ HA_ImP1_1[9]*d33_1+ HA_ImP1_1[12]*d53_1+ HA_ImP1_1[13]*d13_1+ HA_ImP1_1[14]*d33_1+ HA_ImP1_1[15]*d73_1+ HA_ImP1_1[16]*d33_1;
    double Asum3Im_P_1_1= HA_ImP_1_1[0]*d331+ HA_ImP_1_1[4]*d131+ HA_ImP_1_1[9]*d331+ HA_ImP_1_1[12]*d531+ HA_ImP_1_1[13]*d131+ HA_ImP_1_1[14]*d331+ HA_ImP_1_1[15]*d731+ HA_ImP_1_1[16]*d331;

    double Asum3Im_0M11  = HA0_ImM11[0]*d311  + HA0_ImM11[4]* d111 + HA0_ImM11[9]* d311 + HA0_ImM11[12]* d511 + HA0_ImM11[13]* d111 + HA0_ImM11[14]* d311 + HA0_ImM11[15]*d711 + HA0_ImM11[16]* d311;
    double Asum3Im_0M_11 = HA0_ImM_11[0]*d31_1+ HA0_ImM_11[4]*d11_1+ HA0_ImM_11[9]*d31_1+ HA0_ImM_11[12]*d51_1+ HA0_ImM_11[13]*d11_1+ HA0_ImM_11[14]*d31_1+ HA0_ImM_11[15]*d71_1+ HA0_ImM_11[16]*d31_1;
    double Asum3Im_0M1_1 = HA0_ImM1_1[0]*d31_1+ HA0_ImM1_1[4]*d11_1+ HA0_ImM1_1[9]*d31_1+ HA0_ImM1_1[12]*d51_1+ HA0_ImM1_1[13]*d11_1+ HA0_ImM1_1[14]*d31_1+ HA0_ImM1_1[15]*d71_1+ HA0_ImM1_1[16]*d31_1;
    double Asum3Im_0M_1_1= HA0_ImM_1_1[0]*d311+ HA0_ImM_1_1[4]*d111+ HA0_ImM_1_1[9]*d311+ HA0_ImM_1_1[12]*d511+ HA0_ImM_1_1[13]*d111+ HA0_ImM_1_1[14]*d311+ HA0_ImM_1_1[15]*d711+ HA0_ImM_1_1[16]*d311;

    double Asum3Im_0P11  = HA0_ImP11[0]*d311  + HA0_ImP11[4]* d111 + HA0_ImP11[9]* d311 + HA0_ImP11[12]* d511 + HA0_ImP11[13]* d111 + HA0_ImP11[14]* d311 + HA0_ImP11[15]*d711 + HA0_ImP11[16]* d311;
    double Asum3Im_0P_11 = HA0_ImP_11[0]*d31_1+ HA0_ImP_11[4]*d11_1+ HA0_ImP_11[9]*d31_1+ HA0_ImP_11[12]*d51_1+ HA0_ImP_11[13]*d11_1+ HA0_ImP_11[14]*d31_1+ HA0_ImP_11[15]*d71_1+ HA0_ImP_11[16]*d31_1;
    double Asum3Im_0P1_1 = HA0_ImP1_1[0]*d31_1+ HA0_ImP1_1[4]*d11_1+ HA0_ImP1_1[9]*d31_1+ HA0_ImP1_1[12]*d51_1+ HA0_ImP1_1[13]*d11_1+ HA0_ImP1_1[14]*d31_1+ HA0_ImP1_1[15]*d71_1+ HA0_ImP1_1[16]*d31_1;
    double Asum3Im_0P_1_1= HA0_ImP_1_1[0]*d311+ HA0_ImP_1_1[4]*d111+ HA0_ImP_1_1[9]*d311+ HA0_ImP_1_1[12]*d511+ HA0_ImP_1_1[13]*d111+ HA0_ImP_1_1[14]*d311+ HA0_ImP_1_1[15]*d711+ HA0_ImP_1_1[16]*d311;

    // HV-HA
    double sum3Re_M11  = Vsum3Re_M11  - Asum3Re_M11;
    double sum3Re_M_11 = Vsum3Re_M_11 - Asum3Re_M_11 ;
    double sum3Re_M1_1 = Vsum3Re_M1_1 - Asum3Re_M1_1;
    double sum3Re_M_1_1= Vsum3Re_M_1_1- Asum3Re_M_1_1 ;

    double sum3Re_P11  = Vsum3Re_P11  - Asum3Re_P11;
    double sum3Re_P_11 = Vsum3Re_P_11 - Asum3Re_P_11 ;
    double sum3Re_P1_1 = Vsum3Re_P1_1 - Asum3Re_P1_1;
    double sum3Re_P_1_1= Vsum3Re_P_1_1- Asum3Re_P_1_1 ;

    double sum3Re_0M11  = Vsum3Re_0M11  - Asum3Re_0M11;
    double sum3Re_0M_11 = Vsum3Re_0M_11 - Asum3Re_0M_11 ;
    double sum3Re_0M1_1 = Vsum3Re_0M1_1 - Asum3Re_0M1_1;
    double sum3Re_0M_1_1= Vsum3Re_0M_1_1- Asum3Re_0M_1_1 ;

    double sum3Re_0P11  = Vsum3Re_0P11  - Asum3Re_0P11;
    double sum3Re_0P_11 = Vsum3Re_0P_11 - Asum3Re_0P_11 ;
    double sum3Re_0P1_1 = Vsum3Re_0P1_1 - Asum3Re_0P1_1;
    double sum3Re_0P_1_1= Vsum3Re_0P_1_1- Asum3Re_0P_1_1 ;


    double sum3Im_M11  = Vsum3Im_M11  - Asum3Im_M11;
    double sum3Im_M_11 = Vsum3Im_M_11 - Asum3Im_M_11 ;
    double sum3Im_M1_1 = Vsum3Im_M1_1 - Asum3Im_M1_1;
    double sum3Im_M_1_1= Vsum3Im_M_1_1- Asum3Im_M_1_1 ;

    double sum3Im_P11  = Vsum3Im_P11  - Asum3Im_P11;
    double sum3Im_P_11 = Vsum3Im_P_11 - Asum3Im_P_11 ;
    double sum3Im_P1_1 = Vsum3Im_P1_1 - Asum3Im_P1_1;
    double sum3Im_P_1_1= Vsum3Im_P_1_1- Asum3Im_P_1_1 ;

    double sum3Im_0M11  = Vsum3Im_0M11  - Asum3Im_0M11;
    double sum3Im_0M_11 = Vsum3Im_0M_11 - Asum3Im_0M_11 ;
    double sum3Im_0M1_1 = Vsum3Im_0M1_1 - Asum3Im_0M1_1;
    double sum3Im_0M_1_1= Vsum3Im_0M_1_1- Asum3Im_0M_1_1 ;

    double sum3Im_0P11  = Vsum3Im_0P11  - Asum3Im_0P11;
    double sum3Im_0P_11 = Vsum3Im_0P_11 - Asum3Im_0P_11 ;
    double sum3Im_0P1_1 = Vsum3Im_0P1_1 - Asum3Im_0P1_1;
    double sum3Im_0P_1_1= Vsum3Im_0P_1_1- Asum3Im_0P_1_1 ;

    ///////////////////////////////////////////////////////////////////////////////////////////

    *diff_Xsec =  //d sigma/dQ2 dW dtheta dphi

      ((0.974)*(0.974)*(1./Gvcm)*(1./Gvcm)*(1e13)*G_F*G_F/(64.*2.*pi*pi*pi*pi)) * (abs_mom_q/(abs_mom_k_L*abs_mom_k_L))
      *(  
          C_T_plus_square*( 
            3.*(sum3Im_M11*sum3Im_M11 + sum3Im_M_11*sum3Im_M_11 + sum3Im_M1_1*sum3Im_M1_1 + sum3Im_M_1_1*sum3Im_M_1_1) 

            +  pow( ((F_minus11   - G_minus11  ) + sqrt_chk(3.)*sum3Re_M11) , 2)
            +  pow( ((F_minus_11  - G_minus_11 ) + sqrt_chk(3.)*sum3Re_M_11) , 2) 
            +  pow( ((F_minus1_1  - G_minus1_1 ) + sqrt_chk(3.)*sum3Re_M1_1 ) , 2)
            +  pow( ((F_minus_1_1 - G_minus_1_1) + sqrt_chk(3.)*sum3Re_M_1_1 ) , 2)
            )
          + C_T_minus_square*(
            3.*(sum3Im_P11*sum3Im_P11 + sum3Im_P_11*sum3Im_P_11 + sum3Im_P1_1*sum3Im_P1_1 + sum3Im_P_1_1*sum3Im_P_1_1) 

            +  pow( ((F_plus11   - G_plus11  ) + sqrt_chk(3.)*sum3Re_P11) , 2)
            +  pow( ((F_plus_11  - G_plus_11 ) + sqrt_chk(3.)*sum3Re_P_11) , 2) 
            +  pow( ((F_plus1_1  - G_plus1_1 ) + sqrt_chk(3.)*sum3Re_P1_1 ) , 2)
            +  pow( ((F_plus_1_1 - G_plus_1_1) + sqrt_chk(3.)*sum3Re_P_1_1 ) , 2)
            )
          + C_S_minus_square*(
            3.*(sum3Im_0M11*sum3Im_0M11 + sum3Im_0M_11*sum3Im_0M_11 + sum3Im_0M1_1*sum3Im_0M1_1 + sum3Im_0M_1_1*sum3Im_0M_1_1)

            +  pow( ((F_zero_minus11   - G_zero_minus11  ) + sqrt_chk(3.)*sum3Re_0M11) , 2)
            +  pow( ((F_zero_minus_11  - G_zero_minus_11 ) + sqrt_chk(3.)*sum3Re_0M_11) , 2) 
            +  pow( ((F_zero_minus1_1  - G_zero_minus1_1 ) + sqrt_chk(3.)*sum3Re_0M1_1 ) , 2)
            +  pow( ((F_zero_minus_1_1 - G_zero_minus_1_1) + sqrt_chk(3.)*sum3Re_0M_1_1 ) , 2)
            )
          + C_S_plus_square*(
              3.*(sum3Im_0P11*sum3Im_0P11 + sum3Im_0P_11*sum3Im_0P_11 + sum3Im_0P1_1*sum3Im_0P1_1 + sum3Im_0P_1_1*sum3Im_0P_1_1)

              +  pow( ((F_zero_plus11   - G_zero_plus11  ) + sqrt_chk(3.)*sum3Re_0P11) , 2)
              +  pow( ((F_zero_plus_11  - G_zero_plus_11 ) + sqrt_chk(3.)*sum3Re_0P_11) , 2) 
              +  pow( ((F_zero_plus1_1  - G_zero_plus1_1 ) + sqrt_chk(3.)*sum3Re_0P1_1 ) , 2)
              +  pow( ((F_zero_plus_1_1 - G_zero_plus_1_1) + sqrt_chk(3.)*sum3Re_0P_1_1 ) , 2)
              )

          + 2.*cos (phi)*(  

              + C_L_minus*sqrt_chk(C_S_minus_square)*( 

                3.*(sum3Im_M11*sum3Im_0M11 + sum3Im_M_11*sum3Im_0M_11 + sum3Im_M1_1*sum3Im_0M1_1 + sum3Im_M_1_1*sum3Im_0M_1_1) 

                +   ((F_minus11   - G_minus11  ) + sqrt_chk(3.)*sum3Re_M11)* ((F_zero_minus11   - G_zero_minus11  ) + sqrt_chk(3.)*sum3Re_0M11)
                +   ((F_minus_11  - G_minus_11 ) + sqrt_chk(3.)*sum3Re_M_11)* ((F_zero_minus_11  - G_zero_minus_11 ) + sqrt_chk(3.)*sum3Re_0M_11)
                +   ((F_minus1_1  - G_minus1_1 ) + sqrt_chk(3.)*sum3Re_M1_1 )*((F_zero_minus1_1  - G_zero_minus1_1 ) + sqrt_chk(3.)*sum3Re_0M1_1 ) 
                +   ((F_minus_1_1 - G_minus_1_1) + sqrt_chk(3.)*sum3Re_M_1_1)* ((F_zero_minus_1_1 - G_zero_minus_1_1) + sqrt_chk(3.)*sum3Re_0M_1_1 )
                )
              + C_R_minus*sqrt_chk(C_S_minus_square)  *( 

                3.*(sum3Im_P11*sum3Im_0M11 + sum3Im_P_11*sum3Im_0M_11 + sum3Im_P1_1*sum3Im_0M1_1 + sum3Im_P_1_1*sum3Im_0M_1_1) 

                +   ((F_plus11   - G_plus11  ) + sqrt_chk(3.)*sum3Re_P11)*((F_zero_minus11   - G_zero_minus11  ) + sqrt_chk(3.)*sum3Re_0M11)
                +   ((F_plus_11  - G_plus_11 ) + sqrt_chk(3.)*sum3Re_P_11)*((F_zero_minus_11  - G_zero_minus_11 ) + sqrt_chk(3.)*sum3Re_0M_11)
                +   ((F_plus1_1  - G_plus1_1 ) + sqrt_chk(3.)*sum3Re_P1_1)*((F_zero_minus1_1  - G_zero_minus1_1 ) + sqrt_chk(3.)*sum3Re_0M1_1 )
                +   ((F_plus_1_1 - G_plus_1_1) + sqrt_chk(3.)*sum3Re_P_1_1)*((F_zero_minus_1_1 - G_zero_minus_1_1) + sqrt_chk(3.)*sum3Re_0M_1_1 ) 
                )
              + C_L_plus*sqrt_chk(C_S_plus_square)  *( 

                3.*(sum3Im_M11*sum3Im_0P11 + sum3Im_M_11*sum3Im_0P_11 + sum3Im_M1_1*sum3Im_0P1_1 + sum3Im_M_1_1*sum3Im_0P_1_1) 

                +   ((F_minus11   - G_minus11  ) + sqrt_chk(3.)*sum3Re_M11)* ((F_zero_plus11   - G_zero_plus11  ) + sqrt_chk(3.)*sum3Re_0P11)
                +   ((F_minus_11  - G_minus_11 ) + sqrt_chk(3.)*sum3Re_M_11)* ((F_zero_plus_11  - G_zero_plus_11 ) + sqrt_chk(3.)*sum3Re_0P_11)
                +   ((F_minus1_1  - G_minus1_1 ) + sqrt_chk(3.)*sum3Re_M1_1 )*((F_zero_plus1_1  - G_zero_plus1_1 ) + sqrt_chk(3.)*sum3Re_0P1_1 ) 
                +   ((F_minus_1_1 - G_minus_1_1) + sqrt_chk(3.)*sum3Re_M_1_1)* ((F_zero_plus_1_1 - G_zero_plus_1_1) + sqrt_chk(3.)*sum3Re_0P_1_1 )
                )
              + C_R_plus*sqrt_chk(C_S_plus_square)*( 

                  3.*(sum3Im_P11*sum3Im_0P11 + sum3Im_P_11*sum3Im_0P_11 + sum3Im_P1_1*sum3Im_0P1_1 + sum3Im_P_1_1*sum3Im_0P_1_1) 

                  +   ((F_plus11   - G_plus11  ) + sqrt_chk(3.)*sum3Re_P11)*((F_zero_plus11   - G_zero_plus11  ) + sqrt_chk(3.)*sum3Re_0P11)
                  +   ((F_plus_11  - G_plus_11 ) + sqrt_chk(3.)*sum3Re_P_11)*((F_zero_plus_11  - G_zero_plus_11 ) + sqrt_chk(3.)*sum3Re_0P_11)
                  +   ((F_plus1_1  - G_plus1_1 ) + sqrt_chk(3.)*sum3Re_P1_1)*((F_zero_plus1_1  - G_zero_plus1_1 ) + sqrt_chk(3.)*sum3Re_0P1_1 )
                  +   ((F_plus_1_1 - G_plus_1_1) + sqrt_chk(3.)*sum3Re_P_1_1)*((F_zero_plus_1_1 - G_zero_plus_1_1) + sqrt_chk(3.)*sum3Re_0P_1_1 ) 
                  )

              )

              + 2.*sin (phi) *(       
                  -sqrt_chk(3.)*C_L_minus *sqrt_chk(C_S_minus_square)  *( 

                    sum3Im_M11  *((F_zero_minus11   - G_zero_minus11  ) + sqrt_chk(3.)*sum3Re_0M11)                   
                    + sum3Im_M_11 *((F_zero_minus_11  - G_zero_minus_11 ) + sqrt_chk(3.)*sum3Re_0M_11)                 
                    + sum3Im_M1_1 *((F_zero_minus1_1  - G_zero_minus1_1 ) + sqrt_chk(3.)*sum3Re_0M1_1 )                 
                    + sum3Im_M_1_1*((F_zero_minus_1_1 - G_zero_minus_1_1) + sqrt_chk(3.)*sum3Re_0M_1_1 )


                    -   ((F_minus11   - G_minus11  ) + sqrt_chk(3.)*sum3Re_M11)* sum3Im_0M11
                    -   ((F_minus_11  - G_minus_11 ) + sqrt_chk(3.)*sum3Re_M_11)*  sum3Im_0M_11
                    -   ((F_minus1_1  - G_minus1_1 ) + sqrt_chk(3.)*sum3Re_M1_1 )*sum3Im_0M1_1  
                    -   ((F_minus_1_1 - G_minus_1_1) + sqrt_chk(3.)*sum3Re_M_1_1)* sum3Im_0M_1_1
                    )

                  + sqrt_chk(3.)*C_R_minus*sqrt_chk(C_S_minus_square)*( 

                    sum3Im_P11*((F_zero_minus11   - G_zero_minus11  ) + sqrt_chk(3.)*sum3Re_0M11)                 
                    + sum3Im_P_11*((F_zero_minus_11  - G_zero_minus_11 ) + sqrt_chk(3.)*sum3Re_0M_11)                 
                    + sum3Im_P1_1*((F_zero_minus1_1  - G_zero_minus1_1 ) + sqrt_chk(3.)*sum3Re_0M1_1 )                 
                    + sum3Im_P_1_1*((F_zero_minus_1_1 - G_zero_minus_1_1) + sqrt_chk(3.)*sum3Re_0M_1_1 )

                    -   ((F_plus11   - G_plus11  ) + sqrt_chk(3.)*sum3Re_P11)*sum3Im_0M11 
                    -   ((F_plus_11  - G_plus_11 ) + sqrt_chk(3.)*sum3Re_P_11)*sum3Im_0M_11
                    -   ((F_plus1_1  - G_plus1_1 ) + sqrt_chk(3.)*sum3Re_P1_1)*sum3Im_0M1_1 
                    -   ((F_plus_1_1 - G_plus_1_1) + sqrt_chk(3.)*sum3Re_P_1_1)*sum3Im_0M_1_1 
                    )

                  - sqrt_chk(3.)*C_L_plus*sqrt_chk(C_S_plus_square)  *( 

                      sum3Im_M11*((F_zero_plus11   - G_zero_plus11  ) + sqrt_chk(3.)*sum3Re_0P11)                  
                      + sum3Im_M_11*((F_zero_plus_11  - G_zero_plus_11 ) + sqrt_chk(3.)*sum3Re_0P_11)                 
                      + sum3Im_M1_1*((F_zero_plus1_1  - G_zero_plus1_1 ) + sqrt_chk(3.)*sum3Re_0P1_1 )                   
                      + sum3Im_M_1_1*((F_zero_plus_1_1 - G_zero_plus_1_1) + sqrt_chk(3.)*sum3Re_0P_1_1 )                   

                      -   ((F_minus11   - G_minus11  ) + sqrt_chk(3.)*sum3Re_M11)* sum3Im_0P11 
                      -   ((F_minus_11  - G_minus_11 ) + sqrt_chk(3.)*sum3Re_M_11)*  sum3Im_0P_11 
                      -   ((F_minus1_1  - G_minus1_1 ) + sqrt_chk(3.)*sum3Re_M1_1 )*sum3Im_0P1_1 
                      -   ((F_minus_1_1 - G_minus_1_1) + sqrt_chk(3.)*sum3Re_M_1_1)* sum3Im_0P_1_1
                      )

                  + sqrt_chk(3.)*C_R_plus*sqrt_chk(C_S_plus_square)*( 

                      sum3Im_P11*((F_zero_plus11   - G_zero_plus11  ) + sqrt_chk(3.)*sum3Re_0P11)                  
                      + sum3Im_P_11*((F_zero_plus_11  - G_zero_plus_11 ) + sqrt_chk(3.)*sum3Re_0P_11)                 
                      + sum3Im_P1_1*((F_zero_plus1_1  - G_zero_plus1_1 ) + sqrt_chk(3.)*sum3Re_0P1_1 )                  
                      + sum3Im_P_1_1*((F_zero_plus_1_1 - G_zero_plus_1_1) + sqrt_chk(3.)*sum3Re_0P_1_1 )                 

                      -   ((F_plus11   - G_plus11  ) + sqrt_chk(3.)*sum3Re_P11)*sum3Im_0P11
                      -   ((F_plus_11  - G_plus_11 ) + sqrt_chk(3.)*sum3Re_P_11)*sum3Im_0P_11 
                      -   ((F_plus1_1  - G_plus1_1 ) + sqrt_chk(3.)*sum3Re_P1_1)*sum3Im_0P1_1
                      -   ((F_plus_1_1 - G_plus_1_1) + sqrt_chk(3.)*sum3Re_P_1_1)*sum3Im_0P_1_1 
                      )
                  )

                  + 2.*cos(2.*phi)* C_LR *( 

                      3.*(sum3Im_M11*sum3Im_P11 + sum3Im_M_11*sum3Im_P_11 + sum3Im_M1_1*sum3Im_P1_1 + sum3Im_M_1_1*sum3Im_P_1_1) 

                      +   ((F_minus11   - G_minus11  ) + sqrt_chk(3.)*sum3Re_M11)* ((F_plus11   - G_plus11  ) + sqrt_chk(3.)*sum3Re_P11)
                      +   ((F_minus_11  - G_minus_11 ) + sqrt_chk(3.)*sum3Re_M_11)* ((F_plus_11  - G_plus_11 ) + sqrt_chk(3.)*sum3Re_P_11)
                      +   ((F_minus1_1  - G_minus1_1 ) + sqrt_chk(3.)*sum3Re_M1_1 )*((F_plus1_1  - G_plus1_1 ) + sqrt_chk(3.)*sum3Re_P1_1 ) 
                      +   ((F_minus_1_1 - G_minus_1_1) + sqrt_chk(3.)*sum3Re_M_1_1)* ((F_plus_1_1 - G_plus_1_1) + sqrt_chk(3.)*sum3Re_P_1_1 )
                      )

                  - 2.* sin (2.*phi)*sqrt_chk(3.)* C_LR *( 

                      sum3Im_M11*((F_plus11   - G_plus11  ) + sqrt_chk(3.)*sum3Re_P11)                   
                      + sum3Im_M_11*((F_plus_11  - G_plus_11 ) + sqrt_chk(3.)*sum3Re_P_11)                  
                      + sum3Im_M1_1*((F_plus1_1  - G_plus1_1 ) + sqrt_chk(3.)*sum3Re_P1_1 )                  
                      + sum3Im_M_1_1*((F_plus_1_1 - G_plus_1_1) + sqrt_chk(3.)*sum3Re_P_1_1 )


                      -   ((F_minus11   - G_minus11  ) + sqrt_chk(3.)*sum3Re_M11)* sum3Im_P11
                      -   ((F_minus_11  - G_minus_11 ) + sqrt_chk(3.)*sum3Re_M_11)*  sum3Im_P_11
                      -   ((F_minus1_1  - G_minus1_1 ) + sqrt_chk(3.)*sum3Re_M1_1 )*sum3Im_P1_1  
                      -   ((F_minus_1_1 - G_minus_1_1) + sqrt_chk(3.)*sum3Re_M_1_1)* sum3Im_P_1_1
                      )
                  );

    return ;
  }

}




