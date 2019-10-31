#include <iostream>

#include "NFortFns.h"

extern "C" {
  FUNCTION_RETURN evdifcrs_();
  FUNCTION_RETURN evpiprob_();
  void  nesetfgparams_();
  void  nefillmodel_();
  void  zexpconfig_();
  void necard_();
  void necardev_();
}
  
using namespace neut;
using namespace neut::rew;

using std::cout;
using std::endl;

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
    //cout << "NFortFns late initialization" << endl;
    static NFortFns::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new NFortFns;
  
    // Default Parameters Used in MC Generation
	/*
    fInstance->XMAQEdef   = 1.21;
    fInstance->XMVQEdef   = 0.84;
    fInstance->KAPPdef    = 1.0;

    fInstance->XMASPIdef  = 1.21;
    fInstance->XMVSPIdef  = 0.84;

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
	
    fInstance->NECOHEPIdef  = 0;  
    fInstance->XMACOHdef  = 1.0;
    fInstance->RAD0NUdef  = 1.0;
    fInstance->fA1COHdef  = 0.0;
    fInstance->fb1COHdef  = 0.0;

    fInstance->NEPDFdef     = 12;
    fInstance->NEBODEKdef   = 1;

    fInstance->MDLQEAFdef   = 1;
    fInstance->MDLQEdef     = 402;

    fInstance->FEFQEdef   = 1.0;
    fInstance->FEFQEHdef  = 1.8;
    fInstance->FEFINELdef = 1.0;
    fInstance->FEFABSdef  = 1.1;
    fInstance->FEFCXdef   = 1.0;
    fInstance->FEFCXHdef  = 1.8;
    fInstance->FEFALLdef  = 1.0;
	*/

    //These two kept for historical reasons
    fInstance->XMARESdef  = 1.21;
    fInstance->XMVRESdef  = 0.84;

    //All below can not be found in NEUT ROOT file at moment
    fInstance->NEPDFdef     = 12;
    fInstance->NEBODEKdef   = 1;

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
void  NFortFns::zexpconfig() { zexpconfig(); }

void NFortFns::Reconfigure() {
  nefillmodel_();
  nesetfgparams_();
  if (nemdls_.mdlqeaf == 5) zexpconfig_();
}

void NFortFns::SetMCDefaultVal(NSyst_t syst, double val) {

  cout << endl << "NFortFns::SetMCDefaultVal(): Changing MC generated default value for systematic: " << NSyst::AsString(syst).c_str() << " from ";

  switch(syst) {

  case ( kXSecTwkDial_MaCCQEshape ) :
  case ( kXSecTwkDial_MaCCQE ) :
    cout << XMAQEdef;
    XMAQEdef = val;

  break;
  
  case ( kXSecTwkDial_AxlFFCCQE ) :
    cout << MDLQEAFdef;
    MDLQEAFdef = val;
    break;
    
  case ( kXSecTwkDial_VecFFCCQE ) :
    cout << MDLQEdef;
    MDLQEdef = val;
    break;
    
  case ( kXSecTwkDial_SCCVecQE ) :
    cout << SCCFVdef;
    SCCFVdef = val;
    break;

  case ( kXSecTwkDial_SCCAxlQE ) :
    cout << SCCFAdef;
    SCCFAdef = val;
    break;
    
  case ( kXSecTwkDial_PsFF ) :
    cout << FPQEdef;
    FPQEdef = val;
    break;

  case ( kSystNucl_CCQEPauliSupViaKF ) :
    cout << KAPPdef;
    KAPPdef = val;
    break;

    //case ( kSystNucl_CCQEFermiSurfMom ) :
    //fPfTwkDial = twk_dial;
    //break;

    //case ( kSystNucl_CCQEBindingEnergy ) :
    //fEbTwkDial = twk_dial;
    //break;
    
  case ( kXSecTwkDial_BYOnOffDIS ) :
    cout << NEBODEKdef;
    NEBODEKdef = val;
    break;

  case ( kXSecTwkDial_FFRES ) :
    cout << NEIFFdef;
    NEIFFdef = val;
    break;

  case ( kXSecTwkDial_TypeRES ) :
    cout << NENRTYPEdef;
    NENRTYPEdef = val;
    break;

  case ( kXSecTwkDial_CA5RES ) :
    cout << RNECA5Idef;
    RNECA5Idef = val;
    break;

  case ( kXSecTwkDial_BgSclRES ) :
    cout << RNEBGSCLdef;
    RNEBGSCLdef = val;
    break;

    //if (NEIFFdef!=0) {
    case ( kXSecTwkDial_MaRES ) :
      cout << XMASPIdef;
      XMASPIdef = val;
      break;
    case ( kXSecTwkDial_MvRES ) :
      cout << XMVSPIdef;
      XMVSPIdef = val;
      break;

    case ( kXSecTwkDial_NECOHEPI ) :
      cout << NECOHEPIdef;
      NECOHEPIdef = val;
    break;
      //} else {
      //case ( kXSecTwkDial_MaRES ) :
      //cout << XMARESdef;
      //XMASPIdef = val;
      //break;
      //case ( kXSecTwkDial_MvRES ) :
      //cout << XMVRESdef;
      //XMVSPIdef = val;
      //break;
      //}

  default:
    cout << "NFortFns::SetMCDefaultVal() Warning: " << NSyst::AsString(syst).c_str() << " not implemented." << endl;
    return;
    break;
  }
  
  cout << " to " << val << endl;
  
}

void NFortFns::SetDefaults_FromCard() {

  // Set up the defaults from a card file
  // See NEUT D card file in https://t2k.org/nd280/datacomp/production006/mcp/neutcards
  necard_();
  necardev_();
  nefillmodel_();

  /*
  // Hack job to get NEUT working for Prod 6T

  // MEC model default is Nieves lepton tensor
  fInstance->MECModeldef = 1;

  // Default Parameters Used in MC Generation
  fInstance->XMAQEdef   = 1.21;
  fInstance->XMVQEdef   = 0.84;
  fInstance->KAPPdef    = 1.0;

  fInstance->XMASPIdef  = 1.21;
  fInstance->XMVSPIdef  = 0.84;

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

  fInstance->NECOHEPIdef  = 0;  
  fInstance->XMACOHdef  = 1.0;
  fInstance->RAD0NUdef  = 1.0;
  fInstance->fA1COHdef  = 0.0;
  fInstance->fb1COHdef  = 0.0;

  fInstance->NEPDFdef     = 12;
  fInstance->NEBODEKdef   = 1;

  fInstance->MDLQEAFdef   = 1;
  fInstance->MDLQEdef     = 402;

  fInstance->FEFQEdef   = 1.0;
  fInstance->FEFQEHdef  = 1.8;
  fInstance->FEFINELdef = 1.0;
  fInstance->FEFABSdef  = 1.1;
  fInstance->FEFCXdef   = 1.0;
  fInstance->FEFCXHdef  = 1.8;
  fInstance->FEFALLdef  = 1.0;

  //These two kept for historical reasons
  fInstance->XMARESdef  = 1.21;
  fInstance->XMVRESdef  = 0.84;

  //All below can not be found in NEUT ROOT file at moment
  fInstance->NEPDFdef     = 12;
  fInstance->NEBODEKdef   = 1;

  fInstance->SCCFVdef = 0.0;
  fInstance->SCCFAdef = 0.0;
  fInstance->FPQEdef = 1.0;

  nemdls_.sccfv = SCCFVdef;
  nemdls_.sccfa = SCCFAdef;
  nemdls_.fpqe = FPQEdef;

  nemdls_.mdlqeaf = MDLQEAFdef;
  nemdls_.mdlqe = MDLQEdef;
  // Should be 0 as the default model is 0 (Rein and Sehgal)
  neutcoh_.necohepi = NECOHEPIdef;

  nemdls_.xmaqe = XMAQEdef;
  nemdls_.xmvqe = XMVQEdef;
  nemdls_.kapp  = KAPPdef;

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
  nemdls_.fa1coh = fA1COHdef;
  nemdls_.fb1coh = fb1COHdef;

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
  neffpr_.fefall  = FEFALLdef;

  // Fixed values
  neffpr_.fefcoh =  1;
  neffpr_.fefqehf = 1;
  neffpr_.fefcxhf = 0;
  neffpr_.fefcohf = 0;
  */
}

void NFortFns::SetDefaults() {
  // Documentation: See neutmodel.h and necard.h 

  /*
     nemdls_.sccfv = SCCFVdef;
     nemdls_.sccfa = SCCFAdef;
     nemdls_.fpqe = FPQEdef;

     nemdls_.mdlqeaf = MDLQEAFdef;
     nemdls_.mdlqe = MDLQEdef;
  // Should be 0 as the default model is 0 (Rein and Sehgal)
  neutcoh_.necohepi = NECOHEPIdef;

  nemdls_.xmaqe = XMAQEdef;
  nemdls_.xmvqe = XMVQEdef;
  nemdls_.kapp  = KAPPdef;

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
  nemdls_.fa1coh = fA1COHdef;
  nemdls_.fb1coh = fb1COHdef;

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
  neffpr_.fefall  = FEFALLdef;

  // Fixed values
  neffpr_.fefcoh =  1;
  neffpr_.fefqehf = 1;
  neffpr_.fefcxhf = 0;
  neffpr_.fefcohf = 0;
  */

    //defaults that still need to be passed through
    neutdis_.nepdf = NEPDFdef;
  neutdis_.nebodek = NEBODEKdef;

  nemdls_.sccfv = SCCFVdef;
  nemdls_.sccfa = SCCFAdef;
  nemdls_.fpqe = FPQEdef;

  neutcard_.nenefmodl = 1;
  neutcard_.nenefmodh = 1;
  neutcard_.nenefkinh = 1;
  neutpiabs_.neabspiemit = 1;

  neutcard_.nusim = 1;
}

void NFortFns::print_fsihist() {

  cout << endl << "-------------------------------" << endl;
  cout << "FSI History Common Block" << endl;
  cout << endl << "Nvert: " << fsihist_.nvert << endl;
  cout << "ivert  iflgvert posvert[3]" << endl;
  for (int ivert=0; ivert<fsihist_.nvert; ivert++) {
    cout << fsihist_.iflgvert[ivert] << " ";
    for (int j=0; j<3; j++)
      cout << fsihist_.posvert[ivert][j] << " ";
    cout << endl;
  }


  cout << endl << "Nvcvert: " << fsihist_.nvcvert << endl;
  cout << "ip ipvert abspvert iverti ivertf dirvert[3]" << endl;
  for (int ip=0; ip<fsihist_.nvcvert; ip++) {
    cout << fsihist_.ipvert[ip] << " " << fsihist_.abspvert[ip] << " " << fsihist_.iverti[ip] << " " << fsihist_.ivertf[ip] << " ";
    for (int j=0; j<3; j++)
      cout << fsihist_.dirvert[ip][j] << " ";
    cout << endl;
  }

  cout << endl << "Ibound = " << posinnuc_.ibound << endl;
  cout << endl << "Fsiprob = " << fsihist_.fsiprob << endl;

  cout << endl;
}


void NFortFns::print_nework() {

  cout << endl << "-------------------------------" << endl;
  cout << "NEWORK Common Block" << endl;
  cout << endl << "Mode = " << nework_.modene << ", Numne = " << nework_.numne << endl;
  cout << "ipne pne[3]" << endl;
  for (int ip=0; ip<nework_.numne; ip++) {
    cout << nework_.ipne[ip] << " ";
    for (int j=0; j<3; j++)
      cout << nework_.pne[ip][j] << " ";
    cout << endl;
  }

  cout << endl;
}

void NFortFns::print_neutcrs() {

  cout << endl << "-------------------------------" << endl;
  cout << "NEUTCRS Common Block" << endl;
  cout << "Crs[3] = " << neutcrscom_.crsx << " " << neutcrscom_.crsy << " " << neutcrscom_.crsz << " " << endl;
  cout << "Crsphi = " << neutcrscom_.crsphi << endl;
  cout << endl;
}

void NFortFns::print_allevent() {
  print_nework();
  print_neutcrs();
  print_fsihist();
  cout << endl << "-------------------------------" << endl;
}


void NFortFns::print_allparams() {

  cout << endl << endl;

  // Input parameters
  cout << endl << "-- Variable Model Parameters --" << endl;

  cout << "nemdls_.xmaqe = " <<         nemdls_.xmaqe   << endl;    
  cout << "nemdls_.xmvqe = " <<        nemdls_.xmvqe     << endl;  
  cout << "nemdls_.kapp = " <<          nemdls_.kapp      << endl;  
  cout << "nemdls_.xmaspi = " <<      nemdls_.xmaspi     << endl;
  cout << "nemdls_.xmvspi = " <<      nemdls_.xmvspi   << endl;
  cout << "nemdls_.xmares = " <<      nemdls_.xmares     << endl;
  cout << "nemdls_.xmvres = " <<      nemdls_.xmvres   << endl;
  cout << "nemdls_.xmacoh = " <<         nemdls_.xmacoh   << endl;    
  cout << "nemdls_.rad0nu = " <<         nemdls_.rad0nu   << endl;    
  cout << "nemdls_.fa1coh = " <<         nemdls_.fa1coh   << endl;    
  cout << "nemdls_.fb1coh = " <<         nemdls_.fb1coh   << endl; 
  cout << "nemdls_.mdlcoh = " <<         nemdls_.mdlcoh   << endl; 
  cout << "neut1pi_.xmanffres = " <<   neut1pi_.xmanffres << endl;
  cout << "neut1pi_.xmvnffres = " <<   neut1pi_.xmvnffres << endl;
  cout << "neut1pi_.xmarsres = " <<   neut1pi_.xmarsres << endl;
  cout << "neut1pi_.xmvrsres = " <<   neut1pi_.xmvrsres << endl;
  cout << "neut1pi_.neiff = " <<   neut1pi_.neiff << endl;
  cout << "neut1pi_.nenrtype = " <<   neut1pi_.nenrtype << endl;
  cout << "neut1pi_.rneca5i = " <<   neut1pi_.rneca5i << endl;
  cout << "neut1pi_.rnebgscl = " <<   neut1pi_.rnebgscl << endl;
  cout << "neutdis_.nepdf = " <<       neutdis_.nepdf   << endl;  
  cout << "neutdis_.nebodek = " <<     neutdis_.nebodek  << endl;  
  cout << "neffpr_.fefqe = " <<     neffpr_.fefqe   << endl;
  cout << "neffpr_.fefqeh = " <<    neffpr_.fefqeh  << endl;
  cout << "neffpr_.fefinel = " <<   neffpr_.fefinel << endl;
  cout << "neffpr_.fefabs = " <<    neffpr_.fefabs  << endl;
  cout << "neffpr_.fefcx = " <<     neffpr_.fefcx   << endl;
  cout << "neffpr_.fefcxh = " <<    neffpr_.fefcxh  << endl;
  cout << "neffpr_.fefall = " <<    neffpr_.fefall  << endl;

  // Fixed parameters
  cout << endl << "-- These depend on target nucleus --" << endl;

  cout << "neuttarget_.numbndn = " << neuttarget_.numbndn << endl;
  cout << "neuttarget_.numbndp = " <<   neuttarget_.numbndp << endl;
  cout << "neuttarget_.numfrep = " <<  neuttarget_.numfrep<< endl;
  cout << "neuttarget_.numatom = " <<  neuttarget_.numatom<< endl;

  cout << "nenupr_.pfsurf = " <<  nenupr_.pfsurf<< endl;
  cout << "nenupr_.pfmax = " <<  nenupr_.pfmax<< endl;
  cout << "nenupr_.vnuini = " << nenupr_.vnuini<< endl;
  cout << "nenupr_.vnufin = " << nenupr_.vnufin<< endl;

  cout << endl << "-- These are fixed parameters [should = ()] --" << endl;

  cout << "nemdls_.mdlqeaf    = (1)" << nemdls_.mdlqeaf << endl;
  cout << "nemdls_.mdlqe      = (402)" << nemdls_.mdlqe   << endl;
  cout << "nemdls_.sccfv      = (0)" << nemdls_.sccfv   << endl;
  cout << "nemdls_.sccfa      = (0)" << nemdls_.sccfa   << endl;
  cout << "nemdls_.fpqe       = (1)" << nemdls_.fpqe   << endl;

  cout << "nenupr_.iformlen = (1)" <<     nenupr_.iformlen   << endl;
  cout << "neutcard_.nefrmflg = (0)" <<  neutcard_.nefrmflg<< endl; 
  cout << "neutcard_.nepauflg = (0)" <<   neutcard_.nepauflg<< endl; 
  cout << "neutcard_.nenefo16 = (0)" <<  neutcard_.nenefo16 << endl;
  cout << "neutcard_.nemodflg = (0)" <<  neutcard_.nemodflg << endl;
  cout << "neutcard_.nenefmodl = (1)" <<  neutcard_.nenefmodl << endl;
  cout << "neutcard_.nenefmodh = (1)" <<  neutcard_.nenefmodh << endl;
  cout << "neutcard_.nenefkinh = (1)" <<  neutcard_.nenefkinh << endl;
  cout << "neutpiabs_.neabspiemit = (1)" <<  neutpiabs_.neabspiemit << endl;
  cout << "neutcard_.nusim = (1)" <<  neutcard_.nusim << endl;

  cout << "neutcoh_.necohepi  = (0)" <<  neutcoh_.necohepi << endl;

  cout << "neffpr_.fefcoh = (1)" <<   neffpr_.fefcoh << endl;
  cout << "neffpr_.fefqehf = (1)" <<  neffpr_.fefqehf<< endl;
  cout << "neffpr_.fefcxhf = (0)" <<  neffpr_.fefcxhf<< endl;
  cout << "neffpr_.fefcohf  = (0)" << neffpr_.fefcohf<< endl;

  cout << endl;
}
