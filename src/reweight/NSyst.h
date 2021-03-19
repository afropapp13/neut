#ifndef _N_SYSTEMATIC_PARAM_H_
#define _N_SYSTEMATIC_PARAM_H_

#ifdef X
#error NSyst.h uses X macros and as such 'X' cannot be defined before preprocessing.
#endif

// #define NSYST_DEBUG

#include <cmath>
#include <string>

// This semantically highlights that the uncertainty for the a dial has no
// meaning
#define NOUNCERT 0

// This uses an 'X macro' pattern to make sure various helper methods do not get
// out of synch, you must be careful to include a new-line escape if modifying
// the below list, and two-forward-slash style comments will break everything!
// The 'columns' are: [enum name, stringified version, neg-uncert, pos-uncert]
#define NEUTREWEIGHT_DIAL_LIST                                                 \
  /******************************** (Q)El dials *****************************/ \
  /*Tweak Ma NCEL, affects dsigma(NCEL)/dQ2 both in shape and normalization */ \
  /*X(kXSecTwkDial_MaNCEL, MaNCEL, 0.165289256, 0.165289256)    */             \
  /*Tweak elastic nucleon form factors (default/dipole -> BBBA07) */           \
  /*X(kXSecTwkDial_AxlFFNCEL, AxlFFNCEL, NOUNCERT, NOUNCERT) */                \
  /*Tweak Ma CCQE, affects dsigma(CCQE)/dQ2 both in shape and normalization*/  \
  X(kXSecTwkDial_MaCCQE, MaCCQE, 0.165289256, 0.165289256)                     \
  /*Tweak elastic nucleon form factors (default/dipole -> BBBA07)*/            \
  X(kXSecTwkDial_AxlFFCCQE, AxlFFCCQE, NOUNCERT, NOUNCERT)                     \
  /*Tweak vector 2nd class current FF*/                                        \
  X(kXSecTwkDial_SCCVecQE, SCCVecQE, 1.0, 1.0)                                 \
  /*Tweak axial 2nd class current FF*/                                         \
  X(kXSecTwkDial_SCCAxlQE, SCCAxlQE, 1.0, 1.0)                                 \
  /*Tweak pseudoscalar FF at Q2=0*/                                            \
  X(kXSecTwkDial_PsFF, PsFF, 0.03, 0.03)                                       \
  /*2 Component Alpha*/                                                        \
  X(kXSecTwkDial_FAxlCCQEAlpha, FAxlCCQEAlpha, 0.547374419586, 0.547374419586) \
  /*2 Component Gamma*/                                                        \
  X(kXSecTwkDial_FAxlCCQEGamma, FAxlCCQEGamma, 0.276619847554, 0.276619847554) \
  /*3 Component Beta*/                                                         \
  X(kXSecTwkDial_FAxlCCQEBeta, FAxlCCQEBeta, 0.164910087277, 0.164910087277)   \
  /*3 Component Theta*/                                                        \
  X(kXSecTwkDial_FAxlCCQETheta, FAxlCCQETheta, 0.179147035786, 0.179147035786) \
  /*N Coeff in ZEx*/                                                           \
  /* X(kXSecTwkDial_FAZExp_NTerms, FAZExp_NTerms, NOUNCERT, NOUNCERT) */       \
  /*Q2 Cut of in z definition*/                                                \
  /* X(kXSecTwkDial_FAZExp_TCut, FAZExp_TCut, 1.0, 1.0) */                     \
  /*Q2 Optimal Value in z definition*/                                         \
  /* X(kXSecTwkDial_FAZExp_T0, FAZExp_T0, 1.0, 1.0) */                         \
  /*Apply Q^{-4} constraint flag*/                                             \
  /* X(kXSecTwkDial_FAZExp_Q4Cut, FAZExp_Q4Cut, NOUNCERT, NOUNCERT)  */        \
  /*Z-Expansion A0 Co-eff (keep fixed for now)*/                               \
  /* X(kXSecTwkDial_FAZExp_A0, FAZExp_A0, 1.0, 1.0) */                         \
  /*Z-Expansion A1 Co-eff*/                                                    \
  X(kXSecTwkDial_FAZExp_A1, FAZExp_A1, 0.0809326542131, 0.0809326542131)       \
  /*Z-Expansion A2 Co-eff*/                                                    \
  X(kXSecTwkDial_FAZExp_A2, FAZExp_A2, 1.44337567298, 1.44337567298)           \
  /*Z-Expansion A3 Co-eff*/                                                    \
  X(kXSecTwkDial_FAZExp_A3, FAZExp_A3, 1.00947725152, 1.00947725152)           \
  /*Z-Expansion A4 Co-eff*/                                                    \
  X(kXSecTwkDial_FAZExp_A4, FAZExp_A4, 1.77410484897, 1.77410484897)           \
  /*Z-Expansion A5 Co-eff (keep fixed for now)*/                               \
  /* X(kXSecTwkDial_FAZExp_A5, FAZExp_A5, 1.0, 1.0) */                         \
  /*Z-Expansion A6 Co-eff (keep fixed for now)*/                               \
  /* X(kXSecTwkDial_FAZExp_A6, FAZExp_A6, 1.0, 1.0) */                         \
  /*Z-Expansion A7 Co-eff (keep fixed for now)*/                               \
  /* X(kXSecTwkDial_FAZExp_A7, FAZExp_A7, 1.0, 1.0) */                         \
  /*Z-Expansion A8 Co-eff (keep fixed for now)*/                               \
  /* X(kXSecTwkDial_FAZExp_A8, FAZExp_A8, 1.0, 1.0) */                         \
  /*Z-Expansion A9 Co-eff (keep fixed for now)*/                               \
  /* X(kXSecTwkDial_FAZExp_A9, FAZExp_A9, 1.0, 1.0) */                         \
  /****************************** SPP dials**** *****************************/ \
  /*Tweak Ma CCRES, affects d2sigma(CCRES)/dWdQ2 both in shape and*/           \
  /* normalization*/                                                           \
  X(kXSecTwkDial_MaRES, MaRES, 0.157894737, 0.157894737)                       \
  /*Tweak Mv CCRES, affects d2sigma(CCRES)/dWdQ2 both in shape and */          \
  /* normalization*/                                                           \
  X(kXSecTwkDial_MvRES, MvRES, 0.119047619, 0.119047619)                       \
  X(kXSecTwkDial_CA5RES, CA5RES, 0.247524752, 0.247524752)                     \
  X(kXSecTwkDial_BgSclRES, BgSclRES, 0.153846154, 0.153846154)                 \
  X(kXSecTwkDial_UseSeparateBgSclLMCPiBar, UseSeparateBgSclLMCPiBar, NOUNCERT, \
    NOUNCERT)                                                                  \
  /*A logically separate BgScl parameter for low momentum single pion */       \
  /* events from anti-neutrinos where the pion is below Cherenkov threshold.*/ \
  X(kXSecTwkDial_BgSclLMCPiBarRES, BgSclLMCPiBarRES, 0.153846154 * 3,          \
    0.153846154 * 3)                                                           \
  /*Tweak the model of coherent pion production*/                              \
  /* X(kXSecTwkDial_NECOHEPI, NECOHEPI, NOUNCERT, NOUNCERT)   */               \
  /*Tweak A1 coef for COH pion production (Berger&Sehgal)*/                    \
  /* X(kXSecTwkDial_fA1COHpi, fA1COHpi, 0.10, 0.10)         */                 \
  /*Tweak b1 coef for COH pion production (Berger&Sehgal)*/                    \
  /* X(kXSecTwkDial_fb1COHpi, fb1COHpi, 0.10, 0.10)   */                       \
  /*Tweak Ma for COH pion production*/                                         \
  /* X(kXSecTwkDial_MaCOHpi, MaCOHpi, 0.50, 0.50) */                           \
  /*Tweak R0 for COH pion production*/                                         \
  /* X(kXSecTwkDial_R0COHpi, R0COHpi, 0.10, 0.10) */                           \
  /*Tweak the inclusive DIS CC normalization*/                                 \
  /* X(kXSecTwkDial_BYOnOffDIS, BYOnOffDIS, NOUNCERT, NOUNCERT)   */           \
  /*Tweak absorption probability for low energy pions*/                        \
  X(kCascTwkDial_FrAbs_pi, FrAbs_pi, 0.50, 0.50)                               \
  /*Tweak inelastic (QE in NEUT) probability for low energy pions*/            \
  X(kCascTwkDial_FrInelLow_pi, FrInelLow_pi, 0.50, 0.50)                       \
  /*Tweak charge exchange probability for low energy pions*/                   \
  X(kCascTwkDial_FrCExLow_pi, FrCExLow_pi, 0.50, 0.50)                         \
  /*Tweak absorption probability for high energy pions*/                       \
  X(kCascTwkDial_FrInelHigh_pi, FrInelHigh_pi, 0.30, 0.30)                     \
  /*Tweak inelastic (QE in NEUT) probability for high energy pions*/           \
  X(kCascTwkDial_FrCExHigh_pi, FrCExHigh_pi, 0.30, 0.30)                       \
  /*Tweak charge exchange probability for high energy pions*/                  \
  X(kCascTwkDial_FrPiProd_pi, FrPiProd_pi, 0.50, 0.50)                         \
  /**/                                                                         \
  X(kCascTwkDial_TotalProb_N, TotalProb_N, 0.1, 0.1)                           \
  /**/                                                                         \
  X(kCascTwkDial_ElasticProb_N, ElasticProb_N, 0.1, 0.1)                       \
  /**/                                                                         \
  X(kCascTwkDial_SinglePiProb_N, SinglePiProb_N, 0.1, 0.1)                     \
  /**/                                                                         \
  X(kCascTwkDial_DoublePiProb_N, DoublePiProb_N, 0.1, 0.1)

namespace neut {
namespace rew {

#define X(A, B, C, D) A,

typedef enum ENSyst {

  kNullSystematic = 0,

  NEUTREWEIGHT_DIAL_LIST

} NSyst_t;

#undef X

class NSyst {
public:
#define X(A, B, C, D)                                                          \
  case A:                                                                      \
    return #B;

  static std::string AsString(NSyst_t syst) {
    switch (syst) {

      NEUTREWEIGHT_DIAL_LIST
    default:
      return "-";
    }
  }

#undef X

#define X(A, B, C, D)                                                          \
  else if (syst == #B) {                                                       \
    return A;                                                                  \
  }

  static NSyst_t FromString(std::string syst) {

    if (syst == "Null") {
      return kNullSystematic;
    }
    NEUTREWEIGHT_DIAL_LIST
    else {
      throw std::string("Unknown systematic: ") + syst;
    }
  }

#undef X
};

} // namespace rew
} // namespace neut

// Kind of awful set of macros but standardises the use of common reweight
// tooling across the reweighter implementations
#define NSYST_USESDIAL(syst, DN)                                               \
  if (syst == DN) {                                                            \
    return true;                                                               \
  }
#define NSYST_CATVARNAME(VN, DN) VN##DN
#define NSYST_CURRVAR(DN) NSYST_CATVARNAME(fCurr, DN)
#define NSYST_DEFVAR(DN) NSYST_CATVARNAME(fDef, DN)
#define NSYST_TWKVAR(DN) NSYST_CATVARNAME(fTwk, DN)
#define NSYST_ISTWKD(DN) bool(std::abs(NSYST_TWKVAR(DN)) > 1E-8)
#define NSYST_DECLAREDIALVARIABLES(DN)                                         \
  double NSYST_CURRVAR(DN);                                                    \
  double NSYST_DEFVAR(DN);                                                     \
  double NSYST_TWKVAR(DN);
#define NSYST_SETDEF(DN, def)                                                  \
  NSYST_DEFVAR(DN) = def;                                                      \
  NSYST_CURRVAR(DN) = NSYST_DEFVAR(DN);                                        \
  NSYST_TWKVAR(DN) = 0;
#define NSYST_UPDATEVALUE(DN, syst, val)                                       \
  if (syst == DN) {                                                            \
    NSYST_TWKVAR(DN) = val;                                                    \
    return;                                                                    \
  }
#define NSYST_SETTODEF(DN)                                                     \
  NSYST_CURRVAR(DN) = NSYST_DEFVAR(DN);                                        \
  NSYST_TWKVAR(DN) = 0;

#define NSYST_ASSIGNIFTWKD(TGT, DN)                                            \
  if (NSYST_ISTWKD(DN)) {                                                      \
    TGT = NSYST_CURRVAR(DN);                                                   \
  }

#ifndef NSYST_DEBUG
#define NSYST_RECONFCURRVALUE(DN, err)                                         \
  NSYST_CURRVAR(DN) =                                                          \
      NSYST_DEFVAR(DN) *                                                       \
      (1. + NSYST_TWKVAR(DN) *                                                 \
                err->OneSigmaErr(DN, ((NSYST_TWKVAR(DN) > 0) ? +1 : -1)));
#define NSYST_RECONFCURRVALUE_NOUNCERT(DN) NSYST_CURRVAR(DN) = NSYST_TWKVAR(DN);
#define NSYST_GETTWKFORABS(DN, syst, absval, err)                              \
  if (syst == DN) {                                                            \
    if (std::isnormal(                                                         \
            err->OneSigmaErr(DN, ((NSYST_TWKVAR(DN) > 0) ? +1 : -1)))) {       \
      return ((absval / NSYST_DEFVAR(DN)) - 1) /                               \
             err->OneSigmaErr(DN, ((NSYST_TWKVAR(DN) > 0) ? +1 : -1));         \
    } else {                                                                   \
      std::cout << "[WARN] Cannot determine twkval for abs = " << absval       \
                << " for dial " << #DN << std::endl;                           \
    }                                                                          \
  }
#else
#define NSYST_RECONFCURRVALUE(DN, err)                                         \
  NSYST_CURRVAR(DN) =                                                          \
      NSYST_DEFVAR(DN) *                                                       \
      (1. + NSYST_TWKVAR(DN) *                                                 \
                err->OneSigmaErr(DN, ((NSYST_TWKVAR(DN) > 0) ? +1 : -1)));     \
  std::cout << "dial: " << #DN << ", tweak: " << NSYST_TWKVAR(DN)              \
            << " (istwk?" << NSYST_ISTWKD(DN) << ")"                           \
            << ", def: " << NSYST_DEFVAR(DN) << " set to "                     \
            << NSYST_CURRVAR(DN) << " (err = "                                 \
            << err->OneSigmaErr(DN, ((NSYST_TWKVAR(DN) > 0) ? +1 : -1)) << ")" \
            << std::endl;
#define NSYST_RECONFCURRVALUE_NOUNCERT(DN)                                     \
  NSYST_CURRVAR(DN) = NSYST_TWKVAR(DN);                                        \
  std::cout << "dial: " << #DN << ", tweak: " << NSYST_TWKVAR(DN)              \
            << " (istwk?" << NSYST_ISTWKD(DN) << ")"                           \
            << ", def: " << NSYST_DEFVAR(DN) << " set to "                     \
            << NSYST_CURRVAR(DN) << std::endl;
#endif
// Unless explicitly asked for (for use in xmacros in other TUs) these should be
// undefined before the end of this header
#ifndef NEUTREWEIGHT_LEAVE_DIALS_DEFINED
#undef NOUNCERT
#undef NEUTREWEIGHT_DIAL_LIST
#endif

#endif
