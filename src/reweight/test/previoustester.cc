// NEUT ReWeight
#include <iostream>
#include <sstream>
#include "NReWeight.h"
#include "neutvect.h"
#include "NReWeightINuke.h"
#include "NReWeightCasc.h"
#include "NReWeightNuXSecCCQE.h"
#include "NReWeightNuXSecCCRES.h"
#include "NReWeightNuXSecCOH.h"
#include "NReWeightNuXSecDIS.h"
#include "NReWeightNuXSecNC.h"
#include "NReWeightNuXSecNCEL.h"
#include "NReWeightNuXSecNCRES.h"
#include "NReWeightNuXSecRES.h"
#include "NReWeightNuclPiless.h"
#include "NSyst.h"
#include "NSystUncertainty.h"
#include "NReWeightCascNucleon.h"


// Common block variables for NEUT
#include "nefillverC.h"
#include "necardC.h"
#include "neutmodelC.h"
#include "neutparamsC.h"
#include "neworkC.h"
#include "fsihistC.h"
#include "neutcrsC.h"

// NEUT class headers
#include "neutvect.h"
#include "neutpart.h"
#include "neutfsipart.h"
#include "neutfsivert.h"
#include "neutrootTreeSingleton.h"


void FillNeutCommons(NeutVect* nvect) {
  
  // WARNING: This has only been implemented for a neuttree and not GENIE                                                                                            
  // This should be kept in sync with T2KNIWGUtils::GetNIWGEvent(TTree)                                                                                              

  // NEUT version info.  Can't get it to compile properly with this yet                                                                                              
  // neutversion_.corev  =   nvect->COREVer;                                                                                                                         
  // neutversion_.nucev  =   nvect->NUCEVer;                                                                                                                         
  // neutversion_.nuccv  =   nvect->NUCCVer;                                                                                                                         

  // Documentation: See nework.h                                                   // ah Tobes + Clarence 4eva
  //  nrint_.pcascprob = 0.5;
  // nrint_.pnuccounter = 2;
  //not sure about this initialisation yet



  nework_.modene = nvect->Mode;
  nework_.numne = nvect->Npart();

  //  nemdls_.mdlqeaf = nvect->QEVForm;
  nemdls_.mdlqe = nvect->QEModel;
  nemdls_.mdlspi = nvect->SPIModel;
  nemdls_.mdldis = nvect->DISModel;
  nemdls_.mdlcoh = nvect->COHModel;
  neutcoh_.necohepi = nvect->COHModel;

  nemdls_.xmaqe = nvect->QEMA;
  nemdls_.xmvqe = nvect->QEMV;
  nemdls_.kapp = nvect->KAPPA;

  // nemdls_.sccfv = SCCFVdef;                                                                                                                                       
  // nemdls_.sccfa = SCCFAdef;                                                                                                                                       
  // nemdls_.fpqe = FPQEdef;                                                                                                                                         
  
  nemdls_.xmaspi = nvect->SPIMA;
  nemdls_.xmvspi = nvect->SPIMV;
  nemdls_.xmares = nvect->RESMA;
  nemdls_.xmvres = nvect->RESMV;
  
  neut1pi_.xmanffres = nvect->SPIMA;
  neut1pi_.xmvnffres = nvect->SPIMV;
  neut1pi_.xmarsres = nvect->RESMA;
  neut1pi_.xmvrsres = nvect->RESMV;
  neut1pi_.neiff = nvect->SPIForm;
  neut1pi_.nenrtype = nvect->SPINRType;
  neut1pi_.rneca5i = nvect->SPICA5I;
  neut1pi_.rnebgscl = nvect->SPIBGScale;

  nemdls_.xmacoh = nvect->COHMA;
  nemdls_.rad0nu = nvect->COHR0;
  // nemdls_.fa1coh = nvect->COHA1err;                                                                                                                               
  // nemdls_.fb1coh = nvect->COHb1err;                                                                                                                               

  // neutdis_.nepdf = NEPDFdef;                                                                                                                                      
  // neutdis_.nebodek = NEBODEKdef;                                                                                                                                  

  neutcard_.nefrmflg = nvect->FrmFlg;
  neutcard_.nepauflg = nvect->PauFlg;
  neutcard_.nenefo16 = nvect->NefO16;
  neutcard_.nemodflg = nvect->ModFlg;
  // neutcard_.nenefmodl = 1;                                                                                                                                        
  // neutcard_.nenefmodh = 1;                                                                                                                                        
  // neutcard_.nenefkinh = 1;                                                                                                                                        
  // neutpiabs_.neabspiemit = 1;  

  nenupr_.iformlen = nvect->FormLen;

  neutpiless_.ipilessdcy = nvect->IPilessDcy;
  neutpiless_.rpilessdcy = nvect->RPilessDcy;

  neutpiless_.ipilessdcy = nvect->IPilessDcy;
  neutpiless_.rpilessdcy = nvect->RPilessDcy;

  neffpr_.fefqe = nvect->NuceffFactorPIQE;
  neffpr_.fefqeh = nvect->NuceffFactorPIQEH;
  neffpr_.fefinel = nvect->NuceffFactorPIInel;
  neffpr_.fefabs = nvect->NuceffFactorPIAbs;
  neffpr_.fefcx = nvect->NuceffFactorPICX;
  neffpr_.fefcxh = nvect->NuceffFactorPICXH;

  neffpr_.fefcoh = nvect->NuceffFactorPICoh;
  neffpr_.fefqehf = nvect->NuceffFactorPIQEHKin;
  neffpr_.fefcxhf = nvect->NuceffFactorPICXKin;
  neffpr_.fefcohf = nvect->NuceffFactorPIQELKin;

    nrint_.pcascprob = nvect->NrintNucleonCascadeProb;

  /*nucres_.xnucfact = nvect->NucresTotalNucleon;
  nucres_.xnucelafact = nvect->NucresElasticNucleon;
  nucres_.xnucspifact = nvect->NucresSinglePiNucleon;
  nucres_.xnucdpifact = nvect->NucresDoublePiNucleon;
  */
  //nrint_.pcascprob = nvect->NrintNucleonCascadeProb ;
  //nrint_.pnuccounter = nvect->NucleonCascadeInteractionCounter ;
  
  /*  if (nvect->NucFsiStepInfo(0))
	      {
		throw;
	      }
  */
  nucleonfsihist_.nfnstep = nvect->NnucFsiStep();
    nucleonfsihist_.nfnvert = nvect->NnucFsiVert();
    //    std::cout << "nfnstep = " << nucleonfsihist_.nfnstep << " nfnvert = " << nucleonfsihist_.nfnvert << std::endl;
    std::cout << "nfnstep prevtest = " << nucleonfsihist_.nfnstep << std::endl;
    for(int k =0; k<nucleonfsihist_.nfnstep; k++)
      {
	nucleonfsihist_.nfecms2[k]= nvect->NucFsiStepInfo(k)->fECMS2;


	//	std::cout <<"fecms2 = " << nvect->NucFsiStepInfo(k)->fECMS2 << std::endl;
	//	nucleonfsihist_.nfiflag[j] = nvect->NucFsiVertInfo(j)->fVertFlag;


	nucleonfsihist_.nfxstep[k] = (float)nvect->NucFsiStepInfo(k)->fPosStep.X();
	nucleonfsihist_.nfystep[k] = (float)nvect->NucFsiStepInfo(k)->fPosStep.Y();
	nucleonfsihist_.nfzstep[k] = (float)nvect->NucFsiStepInfo(k)->fPosStep.Z(); //divide by 1000?                                                               
        nucleonfsihist_.nfpxstep[k] = (float)nvect->NucFsiStepInfo(k)->fMomStep.Px();
        nucleonfsihist_.nfpystep[k] = (float)nvect->NucFsiStepInfo(k)->fMomStep.Py();
        nucleonfsihist_.nfpzstep[k] = (float)nvect->NucFsiStepInfo(k)->fMomStep.Pz(); //divide by 1000?   
	nucleonfsihist_.nfestep[k] = (float)nvect->NucFsiStepInfo(k)->fMomStep.E();
	nucleonfsihist_.nfiflagstep[k] = (float)nvect->NucFsiStepInfo(k)->fVertFlagStep;
	nucleonfsihist_.nfirhon[k]= (float)nvect->NucFsiStepInfo(k)->fVertFsiRhon;
	nucleonfsihist_.nfipel[k]= (float)nvect->NucFsiStepInfo(k)->fStepPel;
	nucleonfsihist_.nfipsp[k]= (float)nvect->NucFsiStepInfo(k)->fStepPsp;
	nucleonfsihist_.nfipdp[k]= (float)nvect->NucFsiStepInfo(k)->fStepPdp;

	

	std::cout << "step number = " << k << std::endl;
	
	std::cout <<  "nfipel = " << (float)nvect->NucFsiStepInfo(k)->fStepPel << std::endl;
	std::cout <<  "nfipsp = " << (float)nvect->NucFsiStepInfo(k)->fStepPsp << std::endl;
	std::cout <<  "nfipdp = " << (float)nvect->NucFsiStepInfo(k)->fStepPdp << std::endl;

	
	std::cout << "nfxstep = " << (float)nvect->NucFsiStepInfo(k)->fPosStep.X() << " nfiflagstep = " << (float)nvect->NucFsiStepInfo(k)->fVertFlagStep << std::endl;
	std::cout <<  "nfiRHON = " << (float)nvect->NucFsiStepInfo(k)->fVertFsiRhon << std::endl;
	/*
	nucleonfsihist_.nfpx[j] = (float)nvect->NucFsiVertInfo(j)->fMom.X();
	nucleonfsihist_.nfpy[j] = (float)nvect->NucFsiVertInfo(j)->fMom.Y();
	nucleonfsihist_.nfpz[j] = (float)nvect->NucFsiVertInfo(j)->fMom.Z();

	nucleonfsihist_.nfe[j] = (float)nvect->NucFsiVertInfo(j)->fMom.E();
	nucleonfsihist_.nffirststep[j] = (float)nvect->NucFsiVertInfo(j)->fVertFirstStep;
	*/

      }
    for (int j =0; j<nucleonfsihist_.nfnvert; j++)
    {
      nucleonfsihist_.nfiflag[j] = nvect->NucFsiVertInfo(j)->fVertFlag;
      nucleonfsihist_.nfx[j] = (float)nvect->NucFsiVertInfo(j)->fPos.X();
      nucleonfsihist_.nfy[j] = (float)nvect->NucFsiVertInfo(j)->fPos.Y();
      nucleonfsihist_.nfz[j] = (float)nvect->NucFsiVertInfo(j)->fPos.Z(); //divide by 1000?
      
      nucleonfsihist_.nfpx[j] = (float)nvect->NucFsiVertInfo(j)->fMom.X();
      nucleonfsihist_.nfpy[j] = (float)nvect->NucFsiVertInfo(j)->fMom.Y();
      nucleonfsihist_.nfpz[j] = (float)nvect->NucFsiVertInfo(j)->fMom.Z();
     
      nucleonfsihist_.nfe[j] = (float)nvect->NucFsiVertInfo(j)->fMom.E();
      nucleonfsihist_.nffirststep[j] = (float)nvect->NucFsiVertInfo(j)->fVertFirstStep;
      //      nucleonfsihist_.nfecms2[j] = (float)nvect->NucFsiStepInfo(j)->fECMS2;

      //      nucleonfsihist_.nfreweightnucleonflag = 1;
      //      std::cout << "fecms2 =" << nucleonfsihist_.nfecms2[j] << std::endl;
}


  for (int i = 0; i < nework_.numne; i++) {
    nework_.ipne[i] = nvect->PartInfo(i)->fPID;
    nework_.pne[i][0] =
      (float)nvect->PartInfo(i)->fP.X() / 1000;  // VC(NE)WORK in M(G)eV                                                                                             
    nework_.pne[i][1] =
      (float)nvect->PartInfo(i)->fP.Y() / 1000;  // VC(NE)WORK in M(G)eV                                                                                             
    nework_.pne[i][2] =
      (float)nvect->PartInfo(i)->fP.Z() / 1000;  // VC(NE)WORK in M(G)eV                                                                                             
  }
  // fsihist.h                                                                                                                                                       

  // neutroot fills a dummy object for events with no FSI to prevent memory leak                                                                                     
  // when                                                                                                                                                            
  // reading the TTree, so check for it here                                                                                                                         
  std::cout<<"nvect->NFsiVert = "<<nvect->NfsiVert() << std::endl;
  std::cout<<"fsiprob = "<<  nvect->Fsiprob << std::endl;
  //I believe only one vertex is required for nucleon FSI
  //Pion FSI requires 2
  if ((int)nvect->NfsiVert() ==      0) 
    {
      //    if (nvect->NfsiPart()!=1 || nvect->Fsiprob!=-1)                                                                                                            
      //      ERR(WRN) << "T2KNeutUtils::fill_neut_commons(TTree) NfsiPart!=1 or                                                                                       
      //      Fsiprob!=-1 when NfsiVert==1" << std::endl;      
      
      fsihist_.nvert = 0;
      fsihist_.nvcvert = 0;
      fsihist_.fsiprob = 1;
    } else {  // Real FSI event                                                                                                                                        
    fsihist_.nvert = (int)nvect->NfsiVert();
    for (int ivert = 0; ivert < fsihist_.nvert; ivert++) {
      fsihist_.iflgvert[ivert] = nvect->FsiVertInfo(ivert)->fVertID;
      fsihist_.posvert[ivert][0] = (float)nvect->FsiVertInfo(ivert)->fPos.X();
      fsihist_.posvert[ivert][1] = (float)nvect->FsiVertInfo(ivert)->fPos.Y();
      fsihist_.posvert[ivert][2] = (float)nvect->FsiVertInfo(ivert)->fPos.Z();
    }

    fsihist_.nvcvert = nvect->NfsiPart();
    for (int ip = 0; ip < fsihist_.nvcvert; ip++) {
      fsihist_.abspvert[ip] = (float)nvect->FsiPartInfo(ip)->fMomLab;
      fsihist_.abstpvert[ip] = (float)nvect->FsiPartInfo(ip)->fMomNuc;
      fsihist_.ipvert[ip] = nvect->FsiPartInfo(ip)->fPID;
      fsihist_.iverti[ip] = nvect->FsiPartInfo(ip)->fVertStart;
      fsihist_.ivertf[ip] = nvect->FsiPartInfo(ip)->fVertEnd;
      fsihist_.dirvert[ip][0] = (float)nvect->FsiPartInfo(ip)->fDir.X();
      fsihist_.dirvert[ip][1] = (float)nvect->FsiPartInfo(ip)->fDir.Y();
      fsihist_.dirvert[ip][2] = (float)nvect->FsiPartInfo(ip)->fDir.Z();
    }
    fsihist_.fsiprob = nvect->Fsiprob;
  }

  neutcrscom_.crsx = nvect->Crsx;
  neutcrscom_.crsy = nvect->Crsy;
  neutcrscom_.crsz = nvect->Crsz;
  neutcrscom_.crsphi = nvect->Crsphi;
  neutcrscom_.crsq2 = nvect->Crsq2;

  neuttarget_.numbndn = nvect->TargetA - nvect->TargetZ;
  neuttarget_.numbndp = nvect->TargetZ;
  neuttarget_.numfrep = nvect->TargetH;
  neuttarget_.numatom = nvect->TargetA;
  posinnuc_.ibound = nvect->Ibound;

  // put empty nucleon FSI history (since it is not saved in the NeutVect                                                                                            
  // format)                                                                                                                                                         
  // Comment out as NEUT does not have the necessary proton FSI information yet                                                                                      
  //  nucleonfsihist_.nfnvert = 0;                                                                                                                                   
  //  nucleonfsihist_.nfnstep = 0;                                                                                                                                  
}


int main(int argc, char *argv[]) {
  

  if (argc != 2 && argc != 3) {
    std::cerr << "Syntax is: ./tester OUTPUT_NEUT_ROOT_FILE.root (variation)" << std::endl;
    std::cerr << "Variation 0, Variation 1, Variation 2, Variation 3" << std::endl;
    return 1;
  }


  int variation = -1;
  if (argc == 3) variation = int(std::atoi(argv[2]));

  std::string InputName = std::string(argv[1]);

  // List of weight engines
  //
  // NReWeightNuXSecCasc
  // NReWeightNuXSecINuke
  //
  // NReWeightNuXSecCCQE
  // NReWeightNuXSecCCRES
  // NReWeightNuXSecCCCOH
  // NReWeightNuXSecCCDIS
  //
  // NReWeightNuXSecNC
  // NReWeightNuXSecNCEL
  // NReWeightNuXSecNCRES
  // NReWeightNuXSecRES
  // NReWeightNuXSecPiless

  neut::rew::NReWeight rw;
  // rw.AdoptWghtCalc("FSI", new neut::rew::NReWeightINuke);
  //trying pion atm
  //rw.AdoptWghtCalc("FSI", new neut::rew::NReWeightCasc);
  rw.AdoptWghtCalc("FSI", new neut::rew::NReWeightCascNucleon);

  std::vector<neut::rew::NSyst_t> SystematicVector;
  //   SystematicVector.push_back(kINukeTwkDial_MFP_N);
   //SystematicVector.push_back(kXSecTwkDial_CA5CCRES);
  //SystematicVector.push_back(kXSecTwkDial_MaBKGM);
  //for pion, let's use:
  //SystematicVector.push_back(kCascTwkDial_All_pi);
    SystematicVector.push_back(kCascTwkDial_TotalProb); 
    SystematicVector.push_back(kCascTwkDial_ElasticProb);
    SystematicVector.push_back(kCascTwkDial_SinglePiProb);
    SystematicVector.push_back(kCascTwkDial_DoublePiProb);
  // Minoo background
  //kXSecTwkDial_MaBKGM;
  // Pion ejection method
  //kXSecTwkDial_RSPiEj;
  // Single pion model
  //kXSecTwkDial_SPPModel;

  // Initialise the systematic
  for (std::vector<neut::rew::NSyst_t>::iterator it = SystematicVector.begin(); it != SystematicVector.end(); ++it) {
    rw.Systematics().Init(*it);
  }


  // Number of reweights
  const int nSysts = 3;
  double ReWeightVals[nSysts] = {0.};
  
  if (variation == 0) {                                                                                                                                          
    ReWeightVals[0] = 0.5;                                                                                                                                 
    ReWeightVals[1] = 0.5;
    ReWeightVals[2] = 0.5;
    ReWeightVals[3] = 0.5;

  }
  else if (variation == 1) {                                                                                                                                      
    ReWeightVals[0] = 1;                                                                                                                                 

  }
  else if (variation == 2) {                                                                                                                                      
    ReWeightVals[0] = 0;                                                                                                                                  

  }
  else if (variation == 3) {                                                                                                                                      
    ReWeightVals[0] = -0.5;                                                                                                                              

  }                                                                                                                                                                  

  /*   // Number of reweights                                                                                                                                             
  const int nSysts = 1;
  double ReWeightVals[nSysts] = {0., 0., 0.};
    if (variation == 0) {
    ReWeightVals[0] = -1.295074385;
    ReWeightVals[1] = 0.28169004;
    ReWeightVals[2] = 0.476190476;
  } else if (variation == 1) {
    ReWeightVals[0] = -1.295074385;
    ReWeightVals[1] = -0.504066908;
    ReWeightVals[2] = 0.476190476;
  } else if (variation == 2) {
    ReWeightVals[0] = 0.720229192;
    ReWeightVals[1] = -0.30762618;
    ReWeightVals[2] = -0.952380952;
  } else if (variation == 3) {
    ReWeightVals[0] = 0.720229192; 
    ReWeightVals[1] = 0.28169004;
    ReWeightVals[2] = 1.428571429;
  }
  
  */
  // Make some 1D distributions for each mode
  const int nNEUTmodes = 53;
  TH1D *Weighted[nNEUTmodes];
  TH1D *UnWeighted[nNEUTmodes];
  for (int i = 0; i < nNEUTmodes; ++i) {
    std::stringstream ss;
    ss << "_" << i;
    Weighted[i] = new TH1D((std::string("Weighted")+ss.str()).c_str(), (std::string("Weighted")+ss.str()+";p_{#mu} (GeV);d#sigma/dp_{#mu} (cm^2/nucleon/GeV)").c_str(), 100, 0, 2);
    UnWeighted[i] = new TH1D((std::string("UnWeighted")+ss.str()).c_str(), (std::string("UnWeighted")+ss.str()+";p_{#mu} (GeV);d#sigma/dp_{#mu} (cm^2/nucleon/GeV)").c_str(), 100, 0, 2);
  }


  // Set the systematics
  int SystCounter = 0;
  for (std::vector<neut::rew::NSyst_t>::iterator it = SystematicVector.begin(); it != SystematicVector.end(); ++it, ++SystCounter) {
    rw.Systematics().Set(*it, ReWeightVals[SystCounter]);
  }




  // Reconfigure the reweight engine
  rw.Reconfigure();



  

  // Setup the output variables
  double Q2, EnuTrue, W, ppi, Eout, pinit, pn, costhlep, costhpinu, costhpimu, costhpin, weight, ScaleFactor = -999;
  int mode, leppdg, nupdg, pipdg, nucpdg = -999;
  int npart, npi, nnuc, nbad = 0;
  
  std::stringstream ss;
  ss << "var" << variation;
  TFile *OutputFile = new TFile((InputName+"_to_"+ss.str()+".root").c_str(), "RECREATE");
  

  OutputFile->cd();
  //  std::cout << "above treedef";
    TTree *OutputTree = new TTree("VARS", "VARS");
    //  std::cout << "below treedef";
  

  // Now set up the branches
  OutputTree->Branch("EnuTrue", &EnuTrue, "EnuTrue/D");
  OutputTree->Branch("Eout", &Eout, "Eout/D");
  OutputTree->Branch("Q2", &Q2, "Q2/D");
  OutputTree->Branch("W", &W, "W/D");
  OutputTree->Branch("ppi", &ppi, "ppi/D");
  OutputTree->Branch("pn", &pn, "pn/D");
  OutputTree->Branch("pinit", &pinit, "pinit/D");
  OutputTree->Branch("costhlep", &costhlep, "costhlep/D");
  OutputTree->Branch("costhpinu", &costhpinu, "costhpinu/D");
  OutputTree->Branch("costhpimu", &costhpimu, "costhpimu/D");
  OutputTree->Branch("costhpin", &costhpin, "costhpin/D");

  OutputTree->Branch("mode", &mode, "mode/I");
  OutputTree->Branch("nupdg", &nupdg, "nupdg/I");
  OutputTree->Branch("leppdg", &leppdg, "leppdg/I");
  OutputTree->Branch("pipdg", &pipdg, "pipdg/I");
  OutputTree->Branch("nucpdg", &nucpdg, "nucpdg/I");

  OutputTree->Branch("npart", &npart, "npart/I");
  OutputTree->Branch("npi", &npi, "npi/I");
  OutputTree->Branch("nnuc", &nnuc, "nnuc/I");
  OutputTree->Branch("nbad", &nbad, "nbad/I");

  OutputTree->Branch("weight", &weight, "weight/D");
  OutputTree->Branch("scalefactor", &ScaleFactor, "ScaleFactor/D");

  // Open the input
  TFile *Infile = new TFile(InputName.c_str(), "OPEN");
  if(!Infile){
    std::cerr << "Cannot open neutroot file!" << std::endl;
    throw;
  }

  // Open a TTree
  TTree * Tree = (TTree*) Infile->Get("neuttree");
  if(!Tree){
    std::cerr << "Cannot find neuttree!" << std::endl;
    throw;
  }

  NeutVect *nvect = NULL;
  Infile->cd();
  Tree->SetBranchAddress("vectorbranch", &nvect);
  int nEvents = Tree->GetEntriesFast();
  // Print width
  //  int PrintWidth = nEvents/20;
  int PrintWidth = 1;
 std::cout << "survived this long" << std::endl ;
  // Get the scaling factor
  //  TH1D* Flux = (TH1D*)(Infile->Get("flux_numu")->Clone());
  //TH1D* Eventrate = (TH1D*)(Infile->Get("evtrt_numu")->Clone());
  //ScaleFactor = ((Eventrate->Integral("width")*1.E-38)/(nEvents))/Flux->Integral("width");
  ScaleFactor = 100;
  // And each event
  //  for (int j = 0; j < nEvents; ++j) {
  for (int j = 0; j < nEvents; ++j) {
    std::cout << "Event " << j << std::endl;
          

    // Reset variables
    Q2 = -999;
    EnuTrue = -999;
    W = -999;
    ppi = -999;
    Eout = -999;
    pn = -999;
    pinit = -999;
    costhlep = -999;
    costhpinu = -999;
    costhpimu = -999;
    costhpin = -999;
    mode = -999;
    leppdg = -999;
    nupdg = -999;
    pipdg = -999;
    nucpdg = -999;
    weight = -999;

    npart = 0;
    npi = 0;
    nnuc = 0;
    nbad = 0;

    if (j % PrintWidth == 0) {
      std::cout << "On event " << j << "/" << nEvents << " (" << int(double(j)/double(nEvents)*100.0) << "%)" << std::endl;
    }

    Tree->GetEntry(j);
    // Need to fill nework first

    FillNeutCommons(nvect);
    //rw.PrintNeutAll();
    //    rw.nrnuc(j);
    weight = rw.CalcWeight();
    std::cout << "weight = " <<    weight << std::endl;
    // Count how many bad events
    if (weight == 1.0) nbad++;
    
    
    if (weight !=1.0)
      {    std::cout << "weight = " << weight << std::endl;
      }
    
    // Now get some distributions
    TLorentzVector Pnu = nvect->PartInfo(0)->fP;
    TLorentzVector Pinit = nvect->PartInfo(1)->fP;
    TLorentzVector Pout = nvect->PartInfo(2)->fP;
    TLorentzVector Ppi;
    TLorentzVector Pn;

    npart = nvect->Npart();
    EnuTrue = Pnu.E()/1.E3;
    Eout = Pout.E()/1.E3;
    costhlep = cos(Pnu.Vect().Angle(Pout.Vect()));
    mode = nvect->Mode;
    nupdg = nvect->PartInfo(0)->fPID;
    leppdg = nvect->PartInfo(2)->fPID;
    Q2 = -1*(Pnu-Pout)*(Pnu-Pout)/1.E6;
    W = sqrt((Pnu+Pinit-Pout)*(Pnu+Pinit-Pout))/1.E3;
    pinit = Pinit.Vect().Mag();
    std::cout << "Eout = " << Eout << std::endl;
    
    UnWeighted[mode]->Fill(Eout);
    // std::cout << "weight = " << weight << std::endl;
   
    Weighted[mode]->Fill(Eout, weight);
    //over events
    for (int k = 2; k < nvect->Npart(); ++k) {
      int PID = nvect->PartInfo(k)->fPID;
      if (!(nvect->PartInfo(k))->fIsAlive && (nvect->PartInfo(k))->fStatus != 0) continue;
      npart++;
      // Pion
      if (abs(PID) == 211 || PID == 111) {
        npi++;
        // Check highest energy and set
        if (nvect->PartInfo(k)->fP.E() > Ppi.E()) {
          Ppi = nvect->PartInfo(k)->fP;
          pipdg = nvect->PartInfo(k)->fPID;
        }
      } else if (PID == 2112 || PID == 2212) {
        nnuc++;
        if (nvect->PartInfo(k)->fP.E() > Pn.E()) {
          Pn = nvect->PartInfo(k)->fP;
          nucpdg = nvect->PartInfo(k)->fPID;
        }
      }
    } // Finished scanning for pions and nucleons
    
  
    if (npi > 0) {
      ppi = Ppi.Vect().Mag()/1.E3; 
      costhpinu = cos(Ppi.Vect().Angle(Pnu.Vect()));
      costhpimu = cos(Ppi.Vect().Angle(Pout.Vect()));
    }

    if (nnuc > 0) {
      pn = Pn.Vect().Mag()/1.E3;
      costhpin = cos(Ppi.Vect().Angle(Pn.Vect()));
    }

     OutputTree->Fill();
     
  } // Finished looping over events
  
  OutputFile->cd();
  OutputTree->Write();
  for (int j = 0; j < nNEUTmodes; ++j) {
    // Set Poisson errors
    for (int i = 0; i < Weighted[j]->GetNbinsX()+1; ++i) {
      double content1 = Weighted[j]->GetBinContent(i+1);
      double error1 = sqrt(content1);
      Weighted[j]->SetBinError(i+1, error1);
      double content2 = UnWeighted[j]->GetBinContent(i+1);
      double error2 = sqrt(content2);
      UnWeighted[j]->SetBinError(i+1, error2);
    }
    Weighted[j]->Scale(ScaleFactor, "width");
    UnWeighted[j]->Scale(ScaleFactor, "width");
    if (Weighted[j]->Integral() > 0) Weighted[j]->Write();
    if (UnWeighted[j]->Integral() > 0) UnWeighted[j]->Write();
  }

  std::cout << "Finished looping events" << std::endl;

  return 0;
}
