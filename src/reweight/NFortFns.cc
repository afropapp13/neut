#include <iostream>

#include "NFortFns.h"

extern "C" {
  FUNCTION_RETURN evdifcrs_();
  FUNCTION_RETURN evpiprob_();
  void  nesetfgparams_();
  void  nefillmodel_();
  void  zexpconfig_();   // P.S (26.01.17) AxialFF Patch                  
}
  
using namespace neut;
using namespace neut::rew;


NFortFns * NFortFns::fInstance = 0;
//____________________________________________________________________________
NFortFns::NFortFns()
{
}
//____________________________________________________________________________
NFortFns::~NFortFns()
{
  fInstance = 0;
}
//____________________________________________________________________________

NFortFns * NFortFns::Instance()
{
  if(fInstance == 0) {
#ifdef NEUT_REWEIGHT_DEBUG
    std::cout << "NFortFns late initialization" << '\n';
#endif
    static NFortFns::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new NFortFns;
  
    // Default Parameters Used in MC Generation
    fInstance->XMAQEdef   = 1.21;
    fInstance->XMVQEdef   = 0.84;
    fInstance->KAPPdef    = 1.0;

	/*
    fInstance->XMASPIdef  = 1.21;
    fInstance->XMVSPIdef  = 0.84;
	*/
    fInstance->XMARESdef  = 1.21;
    fInstance->XMVRESdef  = 0.84;

    fInstance->XMASPIdef  = 0.95;
    fInstance->XMVSPIdef  = 0.84;
    fInstance->NEIFFdef   = 1;
    fInstance->NENRTYPEdef= 0;
    fInstance->RNECA5Idef = 1.01;
    fInstance->RNEBGSCLdef= 1.30;

    fInstance->XMANFFRESdef = 0.95;
    fInstance->XMVNFFRESdef = 0.84;
    fInstance->XMARSRESdef = 1.21;
    fInstance->XMVRSRESdef = 0.84;
	
    fInstance->XMACOHdef  = 1.0;
    fInstance->RAD0NUdef  = 1.0;

    fInstance->NEPDFdef     = 12;
    fInstance->NEBODEKdef   = 1;

    fInstance->NECOHEPIdef  = 0;
    fInstance->MDLQEAFdef   = 0;
    fInstance->MDLQEdef     = 402;

    fInstance->FEFQEdef   = 1.0;
    fInstance->FEFQEHdef  = 1.8;
    fInstance->FEFINELdef = 1.0;
    fInstance->FEFABSdef  = 1.1;
    fInstance->FEFCXdef   = 1.0;
    fInstance->FEFCXHdef  = 1.8;

    fInstance->SCCFVdef = 0.0;
    fInstance->SCCFAdef = 0.0;
    fInstance->FPQEdef = 1.0;

    fInstance->SetDefaults();
  }
  return fInstance;
}
//____________________________________________________________________________

FUNCTION_RETURN NFortFns::evdifcrs()      {  return evdifcrs_(); }
FUNCTION_RETURN NFortFns::evpiprob()      {  return evpiprob_(); }
void  NFortFns::nesetfgparams() { nesetfgparams_(); }
void  NFortFns::nefillmodel()   { nefillmodel_(); }
void  NFortFns::zexpconfig() { zexpconfig(); } // P.S (26.01.17) AxialFF Patch

void NFortFns::Reconfigure() {
  nefillmodel_();
  nesetfgparams_();
  if (nemdls_.mdlqeaf == 4) zexpconfig_(); // P.S (26.01.17) AxialFF Patch
}

void NFortFns::SetMCDefaultVal(NSyst_t syst, double val) {

# ifdef NEUT_REWEIGHT_DEBUG
  std::cout << "NFortFns::SetMCDefaultVal(): Changing MC generated default value for systematic: " << NSyst::AsString(syst).c_str() << " from ";
# endif

  switch(syst) {

  case ( kXSecTwkDial_MaCCQEshape ) :
  case ( kXSecTwkDial_MaCCQE ) :
#   ifdef NEUT_REWEIGHT_DEBUG
    std::cout << XMAQEdef;
#   endif
    XMAQEdef = val;

  break;
  
  case ( kXSecTwkDial_AxlFFCCQE ) :
#   ifdef NEUT_REWEIGHT_DEBUG
    std::cout << MDLQEAFdef;
#   endif
    MDLQEAFdef = val;
    break;
    
  case ( kXSecTwkDial_VecFFCCQE ) :
#   ifdef NEUT_REWEIGHT_DEBUG
    std::cout << MDLQEdef;
#   endif
    MDLQEdef = val;
    break;
    
  case ( kXSecTwkDial_SCCVecQE ) :
#   ifdef NEUT_REWEIGHT_DEBUG
    std::cout << SCCFVdef;
#   endif
    SCCFVdef = val;
    break;

  case ( kXSecTwkDial_SCCAxlQE ) :
#   ifdef NEUT_REWEIGHT_DEBUG
    std::cout << SCCFAdef;
#   endif
    SCCFAdef = val;
    break;
    
  case ( kXSecTwkDial_PsFF ) :
#   ifdef NEUT_REWEIGHT_DEBUG
    std::cout << FPQEdef;
#   endif
    FPQEdef = val;
    break;

  case ( kSystNucl_CCQEPauliSupViaKF ) :
#   ifdef NEUT_REWEIGHT_DEBUG
    std::cout << KAPPdef;
#   endif
    KAPPdef = val;
    break;

  case ( kXSecTwkDial_BYOnOffDIS ) :
#   ifdef NEUT_REWEIGHT_DEBUG
    std::cout << NEBODEKdef;
#   endif
    NEBODEKdef = val;
    break;

  case ( kXSecTwkDial_FFRES ) :
#   ifdef NEUT_REWEIGHT_DEBUG
    std::cout << NEIFFdef;
#   endif
    NEIFFdef = val;
    break;

  case ( kXSecTwkDial_TypeRES ) :
#   ifdef NEUT_REWEIGHT_DEBUG
    std::cout << NENRTYPEdef;
#   endif
    NENRTYPEdef = val;
    break;

  case ( kXSecTwkDial_CA5RES ) :
#   ifdef NEUT_REWEIGHT_DEBUG
    std::cout << RNECA5Idef;
#   endif
    RNECA5Idef = val;
    break;

  case ( kXSecTwkDial_BgSclRES ) :
#   ifdef NEUT_REWEIGHT_DEBUG
    std::cout << RNEBGSCLdef;
#   endif
    RNEBGSCLdef = val;
    break;

    //if (NEIFFdef!=0) {
    case ( kXSecTwkDial_MaRES ) :
#   ifdef NEUT_REWEIGHT_DEBUG
      std::cout << XMASPIdef;
#   endif
      XMASPIdef = val;
      break;
    case ( kXSecTwkDial_MvRES ) :
#   ifdef NEUT_REWEIGHT_DEBUG
      std::cout << XMVSPIdef;
#   endif
      XMVSPIdef = val;
      break;
    //} else {
    //case ( kXSecTwkDial_MaRES ) :
    //std::cout << XMARESdef;
    //XMASPIdef = val;
    //break;
    //case ( kXSecTwkDial_MvRES ) :
    //std::cout << XMVRESdef;
    //XMVSPIdef = val;
    //break;
    //}

  default:
#   ifdef NEUT_REWEIGHT_DEBUG
    std::cout << "NFortFns::SetMCDefaultVal() Warning: " << NSyst::AsString(syst).c_str() << " not implemented." << '\n';
#   endif
    return;
    break;
  }
  
#   ifdef NEUT_REWEIGHT_DEBUG
  std::cout << " to " << val << '\n';
#   endif
  
}

void NFortFns::SetDefaults() {
  // Documentation: See neutmodel.h and necard.h 

  neutcoh_.necohepi = NECOHEPIdef;
  nemdls_.mdlqeaf = MDLQEAFdef;
  nemdls_.mdlqe = MDLQEdef;

  nemdls_.xmaqe = XMAQEdef;
  nemdls_.xmvqe = XMVQEdef;
  nemdls_.kapp  = KAPPdef;

  nemdls_.sccfv = SCCFVdef;
  nemdls_.sccfa = SCCFAdef;
  nemdls_.fpqe = FPQEdef;

  //  nemdls_.xmaspi = XMASPIdef;
  //  nemdls_.xmvspi = XMVSPIdef;
  //nemdls_.xmares = XMARESdef;
  //nemdls_.xmvres = XMVRESdef;

  //nemdls_.xmaspi    = XMASPIdef;
  //nemdls_.xmvspi    = XMVSPIdef;
  //nemdls_.iffspi    = NEIFFdef;
  //nemdls_.nrtypespi = NENRTYPEdef;
  //nemdls_.rca5ispi  = RNECA5Idef;
  //nemdls_.rbgsclspi = RNEBGSCLdef;
  
  neut1pi_.xmanffres = XMANFFRESdef;
  neut1pi_.xmvnffres = XMVNFFRESdef;
  neut1pi_.xmarsres = XMARSRESdef;
  neut1pi_.xmvrsres = XMVRSRESdef;
  neut1pi_.neiff    = NEIFFdef;
  neut1pi_.nenrtype = NENRTYPEdef;
  neut1pi_.rneca5i  = RNECA5Idef;
  neut1pi_.rnebgscl = RNEBGSCLdef;

  nemdls_.xmacoh = XMACOHdef;
  nemdls_.rad0nu = RAD0NUdef;

  neutdis_.nepdf = NEPDFdef;
  neutdis_.nebodek = NEBODEKdef;

  // Fixed parameters
  neutcard_.nefrmflg  = 0;
  neutcard_.nepauflg  = 0;
  neutcard_.nenefo16  = 0;
  neutcard_.nemodflg  = 0;
  neutcard_.nenefmodl = 1;
  neutcard_.nenefmodh = 1;
  neutcard_.nenefkinh = 1;
  neutpiabs_.neabspiemit = 1;

  //if (MCID==piscat)
  //neutcard_.nusim = 0;
  //else
  neutcard_.nusim = 1;
  
  // Nuclear parameter //
  // Note: Nucleus dependant so set using nesetfgparams_() after nucleus is determined
  //nenupr_.pfsurf      = 0.225;
  //nenupr_.pfmax       = 0.225;
  //nenupr_.vnuini      = -1. * (sqrt(0.9396^2 + pdfsurf^2) - 0.9396);
  //nenupr_.vnufin      =  0;
  nenupr_.iformlen    =  1;  // On by default

  // Pion FSI Tuned Parameters
  neffpr_.fefqe   = FEFQEdef;
  neffpr_.fefqeh  = FEFQEHdef;
  neffpr_.fefinel = FEFINELdef;
  neffpr_.fefabs  = FEFABSdef;
  neffpr_.fefcx   = FEFCXdef;
  neffpr_.fefcxh  = FEFCXHdef;
  
  // Fixed values
  neffpr_.fefcoh =  1;
  neffpr_.fefqehf = 1;
  neffpr_.fefcxhf = 0;
  neffpr_.fefcohf = 0;
}

void NFortFns::print_fsihist() {
  
  std::cout << '\n' << "-------------------------------" << '\n';
  std::cout << "FSI History Common Block" << '\n';
  std::cout << '\n' << "Nvert: " << fsihist_.nvert << '\n';
  std::cout << "ivert  iflgvert posvert[3]" << '\n';
  for (int ivert=0; ivert<fsihist_.nvert; ivert++) {
    std::cout << fsihist_.iflgvert[ivert] << " ";
    for (int j=0; j<3; j++)
      std::cout << fsihist_.posvert[ivert][j] << " ";
    std::cout << '\n';
  }
    
  
  std::cout << '\n' << "Nvcvert: " << fsihist_.nvcvert << '\n';
  std::cout << "ip ipvert abspvert iverti ivertf dirvert[3]" << '\n';
  for (int ip=0; ip<fsihist_.nvcvert; ip++) {
    std::cout << fsihist_.ipvert[ip] << " " << fsihist_.abspvert[ip] << " " << fsihist_.iverti[ip] << " " << fsihist_.ivertf[ip] << " ";
    for (int j=0; j<3; j++)
      std::cout << fsihist_.dirvert[ip][j] << " ";
    std::cout << '\n';
  }

  std::cout << '\n' << "Ibound = " << posinnuc_.ibound << '\n';
  std::cout << '\n' << "Fsiprob = " << fsihist_.fsiprob << '\n';
 
  std::cout << '\n';
}


void NFortFns::print_nework() {
  
  std::cout << '\n' << "-------------------------------" << '\n';
  std::cout << "NEWORK Common Block" << '\n';
  std::cout << '\n' << "Mode = " << nework_.modene << ", Numne = " << nework_.numne << '\n';
  std::cout << "ipne pne[3]" << '\n';
  for (int ip=0; ip<nework_.numne; ip++) {
    std::cout << nework_.ipne[ip] << " ";
    for (int j=0; j<3; j++)
      std::cout << nework_.pne[ip][j] << " ";
    std::cout << '\n';
  }
  
  std::cout << '\n';
}

void NFortFns::print_neutcrs() {
  
  std::cout << '\n' << "-------------------------------" << '\n';
  std::cout << "NEUTCRS Common Block" << '\n';
  std::cout << "Crs[3] = " << neutcrscom_.crsx << " " << neutcrscom_.crsy << " " << neutcrscom_.crsz << " " << '\n';
  std::cout << "Crsphi = " << neutcrscom_.crsphi << '\n';
  std::cout << '\n';
}

void NFortFns::print_allevent() {
  print_nework();
  print_neutcrs();
  print_fsihist();
  std::cout << '\n' << "-------------------------------" << '\n';
}


void NFortFns::print_allparams() {

  std::cout << '\n' << '\n';

  // Input parameters
  std::cout << '\n' << "-- Variable Model Parameters --" << '\n';

  std::cout << "nemdls_.xmaqe = " <<         nemdls_.xmaqe   << '\n';    
  std::cout << "nemdls_.xmvqe = " <<        nemdls_.xmvqe     << '\n';  
  std::cout << "nemdls_.kapp = " <<          nemdls_.kapp      << '\n';  
  std::cout << "nemdls_.xmaspi = " <<      nemdls_.xmaspi     << '\n';
  std::cout << "nemdls_.xmvspi = " <<      nemdls_.xmvspi   << '\n';
  std::cout << "nemdls_.xmares = " <<      nemdls_.xmares     << '\n';
  std::cout << "nemdls_.xmvres = " <<      nemdls_.xmvres   << '\n';
  std::cout << "nemdls_.xmacoh = " <<         nemdls_.xmacoh   << '\n';    
  std::cout << "nemdls_.rad0nu = " <<         nemdls_.rad0nu   << '\n';    
  std::cout << "neut1pi_.xmanffres = " <<   neut1pi_.xmanffres << '\n';
  std::cout << "neut1pi_.xmvnffres = " <<   neut1pi_.xmvnffres << '\n';
  std::cout << "neut1pi_.xmarsres = " <<   neut1pi_.xmarsres << '\n';
  std::cout << "neut1pi_.xmvrsres = " <<   neut1pi_.xmvrsres << '\n';
  std::cout << "neut1pi_.neiff = " <<   neut1pi_.neiff << '\n';
  std::cout << "neut1pi_.nenrtype = " <<   neut1pi_.nenrtype << '\n';
  std::cout << "neut1pi_.rneca5i = " <<   neut1pi_.rneca5i << '\n';
  std::cout << "neut1pi_.rnebgscl = " <<   neut1pi_.rnebgscl << '\n';
  std::cout << "neutdis_.nepdf = " <<       neutdis_.nepdf   << '\n';  
  std::cout << "neutdis_.nebodek = " <<     neutdis_.nebodek  << '\n';  
  std::cout << "neffpr_.fefqe = " <<     neffpr_.fefqe   << '\n';
  std::cout << "neffpr_.fefqeh = " <<    neffpr_.fefqeh  << '\n';
  std::cout << "neffpr_.fefinel = " <<   neffpr_.fefinel << '\n';
  std::cout << "neffpr_.fefabs = " <<    neffpr_.fefabs  << '\n';
  std::cout << "neffpr_.fefcx = " <<     neffpr_.fefcx   << '\n';
  std::cout << "neffpr_.fefcxh = " <<    neffpr_.fefcxh  << '\n';

  // Fixed parameters
  std::cout << '\n' << "-- These depend on target nucleus --" << '\n';

  std::cout << "neuttarget_.numbndn = " << neuttarget_.numbndn << '\n';
  std::cout << "neuttarget_.numbndp = " <<   neuttarget_.numbndp << '\n';
  std::cout << "neuttarget_.numfrep = " <<  neuttarget_.numfrep<< '\n';
  std::cout << "neuttarget_.numatom = " <<  neuttarget_.numatom<< '\n';
  
  std::cout << "nenupr_.pfsurf = " <<  nenupr_.pfsurf<< '\n';
  std::cout << "nenupr_.pfmax = " <<  nenupr_.pfmax<< '\n';
  std::cout << "nenupr_.vnuini = " << nenupr_.vnuini<< '\n';
  std::cout << "nenupr_.vnufin = " << nenupr_.vnufin<< '\n';

  std::cout << '\n' << "-- These are fixed parameters [should = ()] --" << '\n';

  std::cout << "nemdls_.mdlqeaf    = (0)" << nemdls_.mdlqeaf << '\n';
  std::cout << "nemdls_.mdlqe      = (402)" << nemdls_.mdlqe   << '\n';
  std::cout << "nemdls_.sccfv      = (0)" << nemdls_.sccfv   << '\n';
  std::cout << "nemdls_.sccfa      = (0)" << nemdls_.sccfa   << '\n';
  std::cout << "nemdls_.fpqe       = (1)" << nemdls_.fpqe   << '\n';

  std::cout << "nenupr_.iformlen = (1)" <<     nenupr_.iformlen   << '\n';
  std::cout << "neutcard_.nefrmflg = (0)" <<  neutcard_.nefrmflg<< '\n'; 
  std::cout << "neutcard_.nepauflg = (0)" <<   neutcard_.nepauflg<< '\n'; 
  std::cout << "neutcard_.nenefo16 = (0)" <<  neutcard_.nenefo16 << '\n';
  std::cout << "neutcard_.nemodflg = (0)" <<  neutcard_.nemodflg << '\n';
  std::cout << "neutcard_.nenefmodl = (1)" <<  neutcard_.nenefmodl << '\n';
  std::cout << "neutcard_.nenefmodh = (1)" <<  neutcard_.nenefmodh << '\n';
  std::cout << "neutcard_.nenefkinh = (1)" <<  neutcard_.nenefkinh << '\n';
  std::cout << "neutpiabs_.neabspiemit = (1)" <<  neutpiabs_.neabspiemit << '\n';
  std::cout << "neutcard_.nusim = (1)" <<  neutcard_.nusim << '\n';
  
  std::cout << "neutcoh_.necohepi  = (0)" <<  neutcoh_.necohepi << '\n';

  std::cout << "neffpr_.fefcoh = (1)" <<   neffpr_.fefcoh << '\n';
  std::cout << "neffpr_.fefqehf = (1)" <<  neffpr_.fefqehf<< '\n';
  std::cout << "neffpr_.fefcxhf = (0)" <<  neffpr_.fefcxhf<< '\n';
  std::cout << "neffpr_.fefcohf  = (0)" << neffpr_.fefcohf<< '\n';

  std::cout << '\n';
}
