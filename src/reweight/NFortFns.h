
#include "necardC.h"
#include "neutmodelC.h"
#include "neutparamsC.h"
#include "neworkC.h"
#include "fsihistC.h"
#include "nucleonfsihistC.h"
#include "vcworkC.h"
#include "posinnucC.h"
#include "neutcrsC.h"

#include "neutvect.h"
#include "neutpart.h"
#include "neutfsipart.h"
#include "neutfsivert.h"
#include "neutrootTreeSingleton.h"

#include "NModeDefn.h"

#include "NSyst.h"

/*  for 64 bits OS & g77 compiler -- M. Fechner*/
#if defined(f2cFortran)&&!defined(gFortran)
#define FUNCTION_RETURN double
#else
#define FUNCTION_RETURN float
#endif

#ifndef _N_FORTFNS_H_
#define _N_FORTFNS_H_

namespace neut {
namespace rew   {

 class NFortFns
 { 
 public:
  static NFortFns* Instance (void);
  
  NFortFns();
  //NFortFns(const NFortFns & err);
  ~NFortFns();

  static FUNCTION_RETURN evdifcrs();
  static FUNCTION_RETURN evpiprob();
  static void  nesetfgparams();
  static void  nefillmodel();
  
  static void print_nework();
  static void print_neutcrs();
  static void print_fsihist();
  static void print_allevent();

  static void print_allparams();

  void SetDefaults();
  void Reconfigure();
  
  void SetMCDefaultVal(NSyst_t syst, double val);

  // Default Parameters Used in MC Generation
  float NECOHEPIdef;
  float MDLQEAFdef;
  float MDLQEdef;

  float SCCFVdef;
  float SCCFAdef;
  float FPQEdef;

  float XMAQEdef  ;
  float XMVQEdef  ;
  float KAPPdef   ;
  float EBINDdef  ;
  float PFSURFdef ;
  	      
  /*
  float XMASPIdef ;
  float XMVSPIdef ;
  */
  	      
  float XMASPIdef;
  float XMVSPIdef;
  float XMARESdef;
  float XMVRESdef;
              
  float NEIFFdef;
  float NENRTYPEdef;
  float RNECA5Idef;
  float RNEBGSCLdef;

  float XMANFFRESdef;
  float XMVNFFRESdef;
  float XMARSRESdef;
  float XMVRSRESdef;

  float XMACOHdef ;
  float RAD0NUdef ;
  float fA1COHdef ;
  float fb1COHdef ;
	      
  int NEPDFdef    ;
  int NEBODEKdef  ;
  	         		      
  double FEFQEdef  ;
  double FEFQEHdef ;
  double FEFINELdef;
  double FEFABSdef ;
  double FEFCXdef  ;
  double FEFCXHdef ;
  double FEFALLdef ;

  static NFortFns * fInstance;

  struct Cleaner {
      void DummyMethodAndSilentCompiler() { }
      ~Cleaner() {
         if (NFortFns::fInstance !=0) {
            delete NFortFns::fInstance;
            NFortFns::fInstance = 0;
         }
      }
  };
  friend struct Cleaner;

 };

} // rew   namespace
} // neut namespace

#endif

