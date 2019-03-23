//____________________________________________________________________________
/*!

\class    neut::rew::NSyst_t

\brief    An enumeration of systematic parameters

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

	  Patrick de Perio <pdeperio \at physics.utoronto.ca>
	  University of Toronto

\created  Aug 1, 2009

\cpright  Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _N_SYSTEMATIC_PARAM_H_
#define _N_SYSTEMATIC_PARAM_H_

#include <string>

using std::string;

namespace neut {
namespace rew   {

typedef enum ENSyst {

  kNullSystematic = 0,

  //
  // Neutrino cross section systematics
  // 
  // Note: This list comes from GENIE systematics, many of which are 
  //       not currently implemented in NEUT (commented out)
  //
  //

  // NCEL tweaking parameters:
  kXSecTwkDial_NormNCEL,          ///< tweak NCEL normalization (energy independent)
  //kXSecTwkDial_NormNCELenu,       ///< tweak NCEL normalization (maintains dependence on neutrino energy)
  kXSecTwkDial_MaNCELshape,       ///< tweak Ma NCEL, affects dsigma(NCEL)/dQ2 in shape only (normalized to constant integral). Warning: This overrides kXSecTwkDial_MaNCEL below
  kXSecTwkDial_1overMaNCEL2,       ///< tweak 1/MaNCEL^2, affects dsigma(NCEL)/dQ2. More symmetric response in cross section when assuming Gaussian errors on this parameter
  kXSecTwkDial_MaNCEL,            ///< tweak Ma NCEL, affects dsigma(NCEL)/dQ2 both in shape and normalization
  kXSecTwkDial_AxlFFNCEL,    ///< tweak elastic nucleon form factors (default/dipole -> BBBA07) NOT VALIDATED
  kXSecTwkDial_VecFFNCEL,    ///< tweak elastic nucleon form factors (default/dipole -> BBBA05)
  //kXSecTwkDial_VecFFNCELshape,    ///< tweak elastic nucleon form factors (BBA/default -> dipole) - shape only effect of dsigma(NCEL)/dQ2
  //kXSecTwkDial_EtaNCEL,           ///< tweak NCEL strange axial form factor eta, affects dsigma(NCEL)/dQ2 both in shape and normalization

  // CCQE tweaking parameters:
  kXSecTwkDial_NormCCQE,          ///< tweak CCQE normalization (energy independent)
  //kXSecTwkDial_NormCCQEenu,       ///< tweak CCQE normalization (maintains dependence on neutrino energy)
  kXSecTwkDial_MaCCQEshape,       ///< tweak Ma CCQE, affects dsigma(CCQE)/dQ2 in shape only (normalized to constant integral). Warning: This overrides kXSecTwkDial_MaCCQE below
  kXSecTwkDial_1overMaCCQE2,       ///< tweak 1/MaCCQE^2, affects dsigma(CCQE)/dQ2. More symmetric response in cross section when assuming Gaussian errors on this parameter
  kXSecTwkDial_MaCCQE,            ///< tweak Ma CCQE, affects dsigma(CCQE)/dQ2 both in shape and normalization
  kXSecTwkDial_AxlFFCCQE,    ///< tweak elastic nucleon form factors (default/dipole -> BBBA07) NOT VALIDATED
  kXSecTwkDial_VecFFCCQE,      ///< tweak the MDLQE used in the calculation (default is MDLQE = 402) 
/// Note that this affects the default calculation, so should be changed to whatever was used to generate the file.
  kXSecTwkDial_VecFFCCQE_out,  ///< tweak the MDLQE used to calculate the output. Use in conjunction with kXSecTwkDial_VecFFCCQE for all possibilities.
  //kXSecTwkDial_VecFFCCQEshape,    ///< tweak elastic nucleon form factors (BBA/default -> dipole) - shape only effect of dsigma(CCQE)/dQ2

  kXSecTwkDial_SCCVecQE,    ///< tweak vector 2nd class current FF
  kXSecTwkDial_SCCAxlQE,    ///< tweak axial 2nd class current FF
  kXSecTwkDial_PsFF,    ///< tweak pseudoscalar FF at Q2=0

  // Resonance neutrino-production tweaking parameters:

  kXSecTwkDial_NormRES,         ///< tweak CCRES normalization
  kXSecTwkDial_MaRES,           ///< tweak Ma CCRES, affects d2sigma(CCRES)/dWdQ2 both in shape and normalization
  kXSecTwkDial_MaRESshape,      ///< tweak Ma RES, affects d2sigma(RES)/dWdQ2 in shape only (normalized to constant integral)
  kXSecTwkDial_MvRES,           ///< tweak Mv CCRES, affects d2sigma(CCRES)/dWdQ2 both in shape and normalization

  kXSecTwkDial_FFRES,
  kXSecTwkDial_TypeRES,
  kXSecTwkDial_CA5RES,
  kXSecTwkDial_BgSclRES,
  kXSecTwkDial_MaNFFRES,
  kXSecTwkDial_MvNFFRES,
  kXSecTwkDial_MaRSRES,
  kXSecTwkDial_MvRSRES,

  kXSecTwkDial_NormCCRES,         ///< tweak CCRES normalization
  kXSecTwkDial_MaCCRESshape,      ///< tweak Ma CCRES, affects d2sigma(CCRES)/dWdQ2 in shape only (normalized to constant integral)
  //kXSecTwkDial_MvCCRESshape,      ///< tweak Mv CCRES, affects d2sigma(CCRES)/dWdQ2 in shape only (normalized to constant integral)
  kXSecTwkDial_MaCCRES,           ///< tweak Ma CCRES, affects d2sigma(CCRES)/dWdQ2 both in shape and normalization
  kXSecTwkDial_MvCCRES,           ///< tweak Mv CCRES, affects d2sigma(CCRES)/dWdQ2 both in shape and normalization

  kXSecTwkDial_FFCCRES,
  kXSecTwkDial_TypeCCRES,
  kXSecTwkDial_CA5CCRES,
  kXSecTwkDial_BgSclCCRES,
  kXSecTwkDial_MaNFFCCRES,
  kXSecTwkDial_MvNFFCCRES,
  kXSecTwkDial_MaRSCCRES,
  kXSecTwkDial_MvRSCCRES,

  kXSecTwkDial_NormNCRES,         ///< tweak NCRES normalization
  kXSecTwkDial_MaNCRESshape,      ///< tweak Ma NCRES, affects d2sigma(NCRES)/dWdQ2 in shape only (normalized to constant integral)
  //kXSecTwkDial_MvNCRESshape,      ///< tweak Mv NCRES, affects d2sigma(NCRES)/dWdQ2 in shape only (normalized to constant integral)
  kXSecTwkDial_MaNCRES,           ///< tweak Ma NCRES, affects d2sigma(NCRES)/dWdQ2 both in shape and normalization
  kXSecTwkDial_MvNCRES,           ///< tweak Mv NCRES, affects d2sigma(NCRES)/dWdQ2 both in shape and normalization

  kXSecTwkDial_FFNCRES,
  kXSecTwkDial_TypeNCRES,
  kXSecTwkDial_CA5NCRES,
  kXSecTwkDial_BgSclNCRES,
  kXSecTwkDial_MaNFFNCRES,
  kXSecTwkDial_MvNFFNCRES,
  kXSecTwkDial_MaRSNCRES,
  kXSecTwkDial_MvRSNCRES,

  // Coherent pion production tweaking parameters:
  kXSecTwkDial_NECOHEPI,          ///< tweak the model of coherent pion production
  kXSecTwkDial_MaCOHpi,           ///< tweak Ma for COH pion production
  kXSecTwkDial_R0COHpi,           ///< tweak R0 for COH pion production
  kXSecTwkDial_fA1COHpi,          ///< tweak A1 coef for COH pion production (Berger&Sehgal)
  kXSecTwkDial_fb1COHpi,          ///< tweak b1 coef for COH pion production (Berger&Sehgal)

  // Non-resonance background tweaking parameters:
  //kXSecTwkDial_RvpCC1pi,          ///< tweak the 1pi non-RES bkg in the RES region, for v+p CC
  //kXSecTwkDial_RvpCC2pi,          ///< tweak the 2pi non-RES bkg in the RES region, for v+p CC
  //kXSecTwkDial_RvpNC1pi,          ///< tweak the 1pi non-RES bkg in the RES region, for v+p NC
  //kXSecTwkDial_RvpNC2pi,          ///< tweak the 2pi non-RES bkg in the RES region, for v+p NC
  //kXSecTwkDial_RvnCC1pi,          ///< tweak the 1pi non-RES bkg in the RES region, for v+n CC
  //kXSecTwkDial_RvnCC2pi,          ///< tweak the 2pi non-RES bkg in the RES region, for v+n CC
  //kXSecTwkDial_RvnNC1pi,          ///< tweak the 1pi non-RES bkg in the RES region, for v+n NC
  //kXSecTwkDial_RvnNC2pi,          ///< tweak the 2pi non-RES bkg in the RES region, for v+n NC
  //kXSecTwkDial_RvbarpCC1pi,       ///< tweak the 1pi non-RES bkg in the RES region, for vbar+p CC
  //kXSecTwkDial_RvbarpCC2pi,       ///< tweak the 2pi non-RES bkg in the RES region, for vbar+p CC
  //kXSecTwkDial_RvbarpNC1pi,       ///< tweak the 1pi non-RES bkg in the RES region, for vbar+p NC
  //kXSecTwkDial_RvbarpNC2pi,       ///< tweak the 2pi non-RES bkg in the RES region, for vbar+p NC
  //kXSecTwkDial_RvbarnCC1pi,       ///< tweak the 1pi non-RES bkg in the RES region, for vbar+n CC
  //kXSecTwkDial_RvbarnCC2pi,       ///< tweak the 2pi non-RES bkg in the RES region, for vbar+n CC
  //kXSecTwkDial_RvbarnNC1pi,       ///< tweak the 1pi non-RES bkg in the RES region, for vbar+n NC
  //kXSecTwkDial_RvbarnNC2pi,       ///< tweak the 2pi non-RES bkg in the RES region, for vbar+n NC
  //// DIS tweaking parameters - applied for DIS events with (Q2>Q2o, W>Wo), typically Q2o=1GeV^2, Wo=1.7-2.0GeV
  //kXSecTwkDial_AhtBY,             ///< tweak the Bodek-Yang model parameter A_{ht} - incl. both shape and normalization effect
  //kXSecTwkDial_BhtBY,             ///< tweak the Bodek-Yang model parameter B_{ht} - incl. both shape and normalization effect 
  //kXSecTwkDial_CV1uBY,            ///< tweak the Bodek-Yang model parameter CV1u - incl. both shape and normalization effect 
  //kXSecTwkDial_CV2uBY,            ///< tweak the Bodek-Yang model parameter CV2u - incl. both shape and normalization effect 
  //kXSecTwkDial_AhtBYshape,        ///< tweak the Bodek-Yang model parameter A_{ht} - shape only effect to d2sigma(DIS)/dxdy
  //kXSecTwkDial_BhtBYshape,        ///< tweak the Bodek-Yang model parameter B_{ht} - shape only effect to d2sigma(DIS)/dxdy
  //kXSecTwkDial_CV1uBYshape,       ///< tweak the Bodek-Yang model parameter CV1u - shape only effect to d2sigma(DIS)/dxdy
  //kXSecTwkDial_CV2uBYshape,       ///< tweak the Bodek-Yang model parameter CV2u - shape only effect to d2sigma(DIS)/dxdy
  kXSecTwkDial_NormDIS,         ///< tweak the inclusive DIS CC normalization
  kXSecTwkDial_BYOnOffDIS,         ///< tweak the inclusive DIS CC normalization
  //kXSecTwkDial_RnubarnuCC,        ///< tweak the ratio of \sigma(\bar\nu CC) / \sigma(\nu CC)
  //kXSecTwkDial_DISNuclMod,        ///< tweak DIS nuclear modification (shadowing, anti-shadowing, EMC)
  //
  kXSecTwkDial_NC,                ///< tweak the inclusive NC normalization


  //
  // Hadronization (free nucleon target)
  // 

  //kHadrAGKYTwkDial_xF1pi,         ///< tweak xF distribution for low multiplicity (N + pi) DIS f/s produced by AGKY
  //kHadrAGKYTwkDial_pT1pi,         ///< tweak pT distribution for low multiplicity (N + pi) DIS f/s produced by AGKY


  //
  // Medium-effects to hadronization
  // 

  //kHadrNuclTwkDial_FormZone,         ///< tweak formation zone


  //
  // Intranuclear rescattering systematics.
  // There are 2 sets of parameters:
  // - parameters that control the total rescattering probability, P(total)
  // - parameters that control the fraction of each process (`fate'), given a total rescat. prob., P(fate|total)
  // These parameters are considered separately for pions and nucleons.
  //

  kCascTwkDial_FrAbs_pi,       ///< tweak absorption      probability for pions		    
  kCascTwkDial_FrInelLow_pi,   ///< tweak inelastic (QE in NEUT) probability for low energy pions 
  kCascTwkDial_FrCExLow_pi,    ///< tweak charge exchange probability for low energy pions	    
  kCascTwkDial_FrInelHigh_pi,  ///< tweak inelastic (QE in NEUT) probability for high energy pions
  kCascTwkDial_FrCExHigh_pi,   ///< tweak charge exchange probability for high energy pions	    
  kCascTwkDial_FrPiProd_pi,    ///< tweak pion (hadron) production (inelastic in NEUT) probability for pions
  //kINukeTwkDial_MFP_N,       ///< tweak mean free path for nucleons
  //kINukeTwkDial_FrCEx_N,     ///< tweak charge exchange probability for nucleons, for given total rescattering probability
  //kINukeTwkDial_FrElas_N,    ///< tweak elastic         probability for nucleons, for given total rescattering probability
  //kINukeTwkDial_FrInel_N,    ///< tweak inelastic       probability for nucleons, for given total rescattering probability
  //kINukeTwkDial_FrAbs_N,     ///< tweak absorption      probability for nucleons, for given total rescattering probability
  //kINukeTwkDial_FrPiProd_N,  ///< tweak pion production probability for nucleons, for given total rescattering probability

  //
  // Nuclear model
  // 

  kSystNucl_CCQEPauliSupViaKF,   ///<
  kSystNucl_CCQEFermiSurfMom,   ///<
  kSystNucl_CCQEBindingEnergy,   ///<
  //kSystNucl_CCQEMomDistroFGtoSF, ///<

  //
  // Resonance decays
  // 

  //kRDcyTwkDial_BR1gamma,        ///< tweak Resonance -> X + gamma branching ratio, eg Delta+(1232) -> p gamma
  //kRDcyTwkDial_BR1eta,          ///< tweak Resonance -> X + eta   branching ratio, eg N+(1440) -> p eta
  //kRDcyTwkDial_Theta_Delta2Npi  ///< distort pi angular distribution in Delta -> N + pi

  
  kSystNucl_PilessDcyRES

  //
  // Misc
  // 


} NSyst_t;


class NSyst {
public:
  //......................................................................................
  static string AsString(NSyst_t syst) 
  {
    switch(syst) {
    case ( kXSecTwkDial_NormNCEL         ) : return "NormNCEL";             break;
      //case ( kXSecTwkDial_NormNCELenu      ) : return "NormNCELenu";          break;
    case ( kXSecTwkDial_MaNCEL           ) : return "MaNCEL";               break;
    case ( kXSecTwkDial_MaNCELshape      ) : return "MaNCELshape";          break;
    case ( kXSecTwkDial_1overMaNCEL2      ) : return "1overMaNCEL2";          break;
    case ( kXSecTwkDial_AxlFFNCEL   ) : return "AxlFFNCEL";       break;
    case ( kXSecTwkDial_VecFFNCEL   ) : return "VecFFNCEL";       break;
      //case ( kXSecTwkDial_VecFFNCELshape   ) : return "VecFFNCELshape";       break;
      //case ( kXSecTwkDial_EtaNCEL          ) : return "EtaNCEL";              break;

    case ( kXSecTwkDial_NormCCQE         ) : return "NormCCQE";             break;
      //case ( kXSecTwkDial_NormCCQEenu      ) : return "NormCCQEenu";          break;
    case ( kXSecTwkDial_MaCCQE           ) : return "MaCCQE";               break;
    case ( kXSecTwkDial_MaCCQEshape      ) : return "MaCCQEshape";          break;
    case ( kXSecTwkDial_1overMaCCQE2      ) : return "1overMaCCQE2";          break;
    case ( kXSecTwkDial_AxlFFCCQE     ) : return "AxlFFCCQE";       break;
    case ( kXSecTwkDial_VecFFCCQE     ) : return "VecFFCCQE";       break;
    case ( kXSecTwkDial_VecFFCCQE_out ) : return "VecFFCCQE_out";   break;
      //case ( kXSecTwkDial_VecFFCCQEshape   ) : return "VecFFCCQEshape";       break;
    case ( kXSecTwkDial_SCCVecQE   ) : return "SCCVecQE";       break;
    case ( kXSecTwkDial_SCCAxlQE   ) : return "SCCAxlQE";       break;
    case ( kXSecTwkDial_PsFF   ) : return "PsFF";       break;


    case ( kXSecTwkDial_NormRES        ) : return "NormRES";            break;
    case ( kXSecTwkDial_MaRES          ) : return "MaRES";              break;
    case ( kXSecTwkDial_MaRESshape     ) : return "MaRESshape";         break;
    case ( kXSecTwkDial_MvRES          ) : return "MvRES";              break;

    case ( kXSecTwkDial_NormCCRES        ) : return "NormCCRES";            break;
    case ( kXSecTwkDial_MaCCRESshape     ) : return "MaCCRESshape";         break;
      //case ( kXSecTwkDial_MvCCRESshape     ) : return "MvCCRESshape";         break;
    case ( kXSecTwkDial_MaCCRES          ) : return "MaCCRES";              break;
    case ( kXSecTwkDial_MvCCRES          ) : return "MvCCRES";              break;

    case ( kXSecTwkDial_NormNCRES        ) : return "NormNCRES";            break;
    case ( kXSecTwkDial_MaNCRESshape     ) : return "MaNCRESshape";         break;
      //case ( kXSecTwkDial_MvNCRESshape     ) : return "MvNCRESshape";         break;
    case ( kXSecTwkDial_MaNCRES          ) : return "MaNCRES";              break;
    case ( kXSecTwkDial_MvNCRES          ) : return "MvNCRES";              break;

    case ( kXSecTwkDial_NECOHEPI         ) : return "NECOHEPI";             break;
    case ( kXSecTwkDial_MaCOHpi          ) : return "MaCOHpi";              break;
    case ( kXSecTwkDial_R0COHpi          ) : return "R0COHpi";              break;
    case ( kXSecTwkDial_fA1COHpi          ) : return "A1COHpi";              break;
    case ( kXSecTwkDial_fb1COHpi          ) : return "b1COHpi";              break;

      //case ( kXSecTwkDial_RvpCC1pi         ) : return "NonRESBGvpCC1pi";      break;
      //case ( kXSecTwkDial_RvpCC2pi         ) : return "NonRESBGvpCC2pi";      break;
      //case ( kXSecTwkDial_RvpNC1pi         ) : return "NonRESBGvpNC1pi";      break;
      //case ( kXSecTwkDial_RvpNC2pi         ) : return "NonRESBGvpNC2pi";      break;
      //case ( kXSecTwkDial_RvnCC1pi         ) : return "NonRESBGvnCC1pi";      break;
      //case ( kXSecTwkDial_RvnCC2pi         ) : return "NonRESBGvnCC2pi";      break;
      //case ( kXSecTwkDial_RvnNC1pi         ) : return "NonRESBGvnNC1pi";      break;
      //case ( kXSecTwkDial_RvnNC2pi         ) : return "NonRESBGvnNC2pi";      break;
      //case ( kXSecTwkDial_RvbarpCC1pi      ) : return "NonRESBGvbarpCC1pi";   break;
      //case ( kXSecTwkDial_RvbarpCC2pi      ) : return "NonRESBGvbarpCC2pi";   break;
      //case ( kXSecTwkDial_RvbarpNC1pi      ) : return "NonRESBGvbarpNC1pi";   break;
      //case ( kXSecTwkDial_RvbarpNC2pi      ) : return "NonRESBGvbarpNC2pi";   break;
      //case ( kXSecTwkDial_RvbarnCC1pi      ) : return "NonRESBGvbarnCC1pi";   break;
      //case ( kXSecTwkDial_RvbarnCC2pi      ) : return "NonRESBGvbarnCC2pi";   break;
      //case ( kXSecTwkDial_RvbarnNC1pi      ) : return "NonRESBGvbarnNC1pi";   break;
      //case ( kXSecTwkDial_RvbarnNC2pi      ) : return "NonRESBGvbarnNC2pi";   break;
      //case ( kXSecTwkDial_AhtBY            ) : return "AhtBY";                break;
      //case ( kXSecTwkDial_BhtBY            ) : return "BhtBY";                break;
      //case ( kXSecTwkDial_CV1uBY           ) : return "CV1uBY";               break;
      //case ( kXSecTwkDial_CV2uBY           ) : return "CV2uBY";               break;
      //case ( kXSecTwkDial_AhtBYshape       ) : return "AhtBYshape";           break;
      //case ( kXSecTwkDial_BhtBYshape       ) : return "BhtBYshape";           break;
      //case ( kXSecTwkDial_CV1uBYshape      ) : return "CV1uBYshape";          break;
      //case ( kXSecTwkDial_CV2uBYshape      ) : return "CV2uBYshape";          break;
    case ( kXSecTwkDial_NormDIS            ) : return "NormDIS";              break;
    case ( kXSecTwkDial_BYOnOffDIS            ) : return "ByOnOffDIS";              break;
      //case ( kXSecTwkDial_RnubarnuCC       ) : return "RnubarnuCC";           break;
      //case ( kXSecTwkDial_DISNuclMod       ) : return "DISNuclMod";           break;
    case ( kXSecTwkDial_NC                 ) : return "NC";            break;
      //case ( kHadrAGKYTwkDial_xF1pi        ) : return "AGKYxF1pi";            break;
      //case ( kHadrAGKYTwkDial_pT1pi        ) : return "AGKYpT1pi";            break;
      //case ( kHadrNuclTwkDial_FormZone     ) : return "FormZone";             break;
    case ( kCascTwkDial_FrAbs_pi           ) : return "FrAbs_pi";             break;
    case ( kCascTwkDial_FrInelLow_pi       ) : return "FrInelLow_pi";         break;
    case ( kCascTwkDial_FrCExLow_pi        ) : return "FrCExLow_pi";          break;
    case ( kCascTwkDial_FrInelHigh_pi      ) : return "FrInelHigh_pi";        break;
    case ( kCascTwkDial_FrCExHigh_pi       ) : return "FrCExHigh_pi";         break;
    case ( kCascTwkDial_FrPiProd_pi        ) : return "FrPiProd_pi";          break;
      //case ( kINukeTwkDial_MFP_N           ) : return "MFP_N";                break;
      //case ( kINukeTwkDial_FrCEx_N         ) : return "FrCEx_N";              break;
      //case ( kINukeTwkDial_FrElas_N        ) : return "FrElas_N";             break;
      //case ( kINukeTwkDial_FrInel_N        ) : return "FrInel_N";             break;
      //case ( kINukeTwkDial_FrAbs_N         ) : return "FrAbs_N";              break;
      //case ( kINukeTwkDial_FrPiProd_N      ) : return "FrPiProd_N";           break;
    case ( kSystNucl_CCQEPauliSupViaKF   ) : return "CCQEPauliSupViaKF";    break;
    case ( kSystNucl_CCQEFermiSurfMom   ) : return "CCQEFermiSurfMom";    break;
    case ( kSystNucl_CCQEBindingEnergy   ) : return "CCQEBindingEnergy";    break;
      //case ( kSystNucl_CCQEMomDistroFGtoSF ) : return "CCQEMomDistroFGtoSF";  break;
      //case ( kRDcyTwkDial_BR1gamma         ) : return "RDecBR1gamma";         break;
      //case ( kRDcyTwkDial_BR1eta           ) : return "RDecBR1eta";           break;
      //case ( kRDcyTwkDial_Theta_Delta2Npi  ) : return "Theta_Delta2Npi";      break;

    default: 
      return "-";
    }
    return "";
  }
 //......................................................................................
 static NSyst_t FromString(string syst)
 {
   NSyst_t systematics[] = 
     {
       kXSecTwkDial_NormNCEL,   
       //kXSecTwkDial_NormNCELenu,   
       kXSecTwkDial_MaNCEL,        
       kXSecTwkDial_MaNCELshape,   
       kXSecTwkDial_1overMaNCEL2,   
       kXSecTwkDial_AxlFFNCEL,
       kXSecTwkDial_VecFFNCEL,
       //kXSecTwkDial_VecFFNCELshape,
       //kXSecTwkDial_EtaNCEL,

       kXSecTwkDial_NormCCQE,   
       //kXSecTwkDial_NormCCQEenu,   
       kXSecTwkDial_MaCCQE,        
       kXSecTwkDial_MaCCQEshape,   
       kXSecTwkDial_1overMaCCQE2,   
       kXSecTwkDial_AxlFFCCQE,
       kXSecTwkDial_VecFFCCQE,
       kXSecTwkDial_VecFFCCQE_out,
       //kXSecTwkDial_VecFFCCQEshape,
       kXSecTwkDial_SCCVecQE,
       kXSecTwkDial_SCCAxlQE,
       kXSecTwkDial_PsFF,

       kXSecTwkDial_NormRES,
       kXSecTwkDial_MaRESshape,   
       kXSecTwkDial_MaRES,      
       kXSecTwkDial_MvRES,      

       kXSecTwkDial_NormCCRES,    
       kXSecTwkDial_MaCCRESshape, 
       //kXSecTwkDial_MvCCRESshape, 
       kXSecTwkDial_MaCCRES,      
       kXSecTwkDial_MvCCRES,      

       kXSecTwkDial_NormNCRES,    
       kXSecTwkDial_MaNCRESshape, 
       //kXSecTwkDial_MvNCRESshape, 
       kXSecTwkDial_MaNCRES,      
       kXSecTwkDial_MvNCRES, 
       
       kXSecTwkDial_NECOHEPI,
       kXSecTwkDial_MaCOHpi,      
       kXSecTwkDial_R0COHpi,    
       kXSecTwkDial_fA1COHpi,      
       kXSecTwkDial_fb1COHpi,      

       //kXSecTwkDial_RvpCC1pi,   
       //kXSecTwkDial_RvpCC2pi,   
       //kXSecTwkDial_RvpNC1pi,    
       //kXSecTwkDial_RvpNC2pi,     
       //kXSecTwkDial_RvnCC1pi,     
       //kXSecTwkDial_RvnCC2pi,     
       //kXSecTwkDial_RvnNC1pi,     
       //kXSecTwkDial_RvnNC2pi,     
       //kXSecTwkDial_RvbarpCC1pi,  
       //kXSecTwkDial_RvbarpCC2pi,  
       //kXSecTwkDial_RvbarpNC1pi,  
       //kXSecTwkDial_RvbarpNC2pi,  
       //kXSecTwkDial_RvbarnCC1pi,  
       //kXSecTwkDial_RvbarnCC2pi,  
       //kXSecTwkDial_RvbarnNC1pi,  
       //kXSecTwkDial_RvbarnNC2pi,  
       //kXSecTwkDial_AhtBY,        
       //kXSecTwkDial_BhtBY,        
       //kXSecTwkDial_CV1uBY,       
       //kXSecTwkDial_CV2uBY,       
       //kXSecTwkDial_AhtBYshape,   
       //kXSecTwkDial_BhtBYshape,   
       //kXSecTwkDial_CV1uBYshape,  
       //kXSecTwkDial_CV2uBYshape,  
       kXSecTwkDial_NormDIS,    
       kXSecTwkDial_BYOnOffDIS,    
       //kXSecTwkDial_RnubarnuCC,   
       //kXSecTwkDial_DISNuclMod, 
       kXSecTwkDial_NC,
       //kHadrAGKYTwkDial_xF1pi,    
       //kHadrAGKYTwkDial_pT1pi,    
       //kHadrNuclTwkDial_FormZone, 
       kCascTwkDial_FrAbs_pi,    
       kCascTwkDial_FrInelLow_pi,   
       kCascTwkDial_FrCExLow_pi,    
       kCascTwkDial_FrInelHigh_pi,   
       kCascTwkDial_FrCExHigh_pi,    
       kCascTwkDial_FrPiProd_pi, 
       //kINukeTwkDial_MFP_N,       
       //kINukeTwkDial_FrCEx_N,    
       //kINukeTwkDial_FrElas_N,   
       //kINukeTwkDial_FrInel_N,   
       //kINukeTwkDial_FrAbs_N,    
       //kINukeTwkDial_FrPiProd_N, 
       kSystNucl_CCQEPauliSupViaKF,   
       kSystNucl_CCQEFermiSurfMom,   
       kSystNucl_CCQEBindingEnergy,   
       //kSystNucl_CCQEMomDistroFGtoSF,
       //kRDcyTwkDial_BR1gamma,       
       //kRDcyTwkDial_BR1eta,         
       //kRDcyTwkDial_Theta_Delta2Npi,
       kSystNucl_PilessDcyRES,
       kNullSystematic
   };

   int isyst=0;
   while(systematics[isyst]!=kNullSystematic) {
     if( AsString(systematics[isyst]) == syst ) {
        return systematics[isyst];
     }
     isyst++;
   }

   return kNullSystematic;
 }
 /*
 //......................................................................................
 static bool IsINukeNuclFateSystematic(NSyst_t syst) 
 {
   switch(syst) {
     case ( kINukeTwkDial_FrCEx_N   ) : 
     case ( kINukeTwkDial_FrElas_N  ) : 
     case ( kINukeTwkDial_FrInel_N  ) : 
     case ( kINukeTwkDial_FrAbs_N   ) : 
     case ( kINukeTwkDial_FrPiProd_N) : 
        return true;
        break;
     default: 
        return false;
        break;
   }
   return false;
 }
 //......................................................................................
 static bool IsINukeFateSystematic(NSyst_t syst)
 {
   switch(syst) {
     case ( kINukeTwkDial_FrCEx_N     ) : 
     case ( kINukeTwkDial_FrElas_N    ) :
     case ( kINukeTwkDial_FrInel_N    ) :
     case ( kINukeTwkDial_FrAbs_N     ) :
     case ( kINukeTwkDial_FrPiProd_N  ) :
       return true;
       break;
     
     default:
       return false;
       break;
   }
   return false;
 }
 //......................................................................................
 static bool IsINukeNuclMeanFreePathSystematic(NSyst_t syst)
 {
   switch(syst) {
     case ( kINukeTwkDial_MFP_N  ) : 
       return true;
       break;
     
     default:
       return false;
       break;
   }
   return false;
 }
 //......................................................................................
 static bool IsINukeMeanFreePathSystematic(NSyst_t syst)
 {
   switch(syst) {
     case ( kINukeTwkDial_MFP_N  ) : 
       return true;
       break;
     
     default:
       return false;
       break;
   }
   return false;
 }
 //......................................................................................
 
 */
};

} // rew   namespace
} // neut namespace

#endif 

