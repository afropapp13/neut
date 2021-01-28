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
  //  TH1::SetDefaultSumw2(true);
  



  nework_.modene = nvect->Mode;
  nework_.numne = nvect->Npart();

  //  nemdls_.mdlqeaf = nvect->QEVForm;
  /*  nemdls_.mdlqe = nvect->QEModel;
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
  */
  //Total cascade probability for nuckeon reweight
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
    //set number of steps and vertices
    nucleonfsihist_.nfnstep = nvect->NnucFsiStep();
    nucleonfsihist_.nfnvert = nvect->NnucFsiVert();


    for(int k =0; k<nucleonfsihist_.nfnstep; k++)
      {
	//set step level variables for reweight
      		nucleonfsihist_.nfiflagstep[k] = (float)nvect->NucFsiStepInfo(k)->fVertFlagStep;
	nucleonfsihist_.nfirhon[k]= (float)nvect->NucFsiStepInfo(k)->fVertFsiRhon;
	nucleonfsihist_.nfipel[k]= (float)nvect->NucFsiStepInfo(k)->fStepPel;
	nucleonfsihist_.nfipsp[k]= (float)nvect->NucFsiStepInfo(k)->fStepPsp;
	nucleonfsihist_.nfipdp[k]= (float)nvect->NucFsiStepInfo(k)->fStepPdp;

	


      }

    for (int j =0; j<nucleonfsihist_.nfnvert; j++)
    {
      //set vertex level info for reweight
      nucleonfsihist_.nfiflag[j] = nvect->NucFsiVertInfo(j)->fVertFlag;
      nucleonfsihist_.nfx[j] = (float)nvect->NucFsiVertInfo(j)->fPos.X();
      nucleonfsihist_.nfy[j] = (float)nvect->NucFsiVertInfo(j)->fPos.Y();
      nucleonfsihist_.nfz[j] = (float)nvect->NucFsiVertInfo(j)->fPos.Z(); //divide by 1000?
      
      nucleonfsihist_.nfpx[j] = (float)nvect->NucFsiVertInfo(j)->fMom.X();
      nucleonfsihist_.nfpy[j] = (float)nvect->NucFsiVertInfo(j)->fMom.Y();
      nucleonfsihist_.nfpz[j] = (float)nvect->NucFsiVertInfo(j)->fMom.Z();
     
      nucleonfsihist_.nfe[j] = (float)nvect->NucFsiVertInfo(j)->fMom.E();
      nucleonfsihist_.nffirststep[j] = (float)nvect->NucFsiVertInfo(j)->fVertFirstStep;

}

    /*    //set some other stuff, not sure if this is strictly needed
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


  if ((int)nvect->NfsiVert() ==      0) 
    {

      
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
    */
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

  TH1::SetDefaultSumw2(true);
  if (argc != 2 && argc != 3 && argc != 4) {
    std::cerr << "Syntax is: ./tester OUTPUT_NEUT_ROOT_FILE.root (variation) (nEvts)" << std::endl;
    std::cerr << "Variation 0, Variation 1, Variation 2, Variation 3" << std::endl;
    return 1;
  }


  int variation = -1;
  variation = int(std::atoi(argv[2]));

  int nEvtsToProcess = -1;
  if(argc==4) nEvtsToProcess = int(std::atoi(argv[3]));


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

  rw.AdoptWghtCalc("FSI", new neut::rew::NReWeightCascNucleon);
  //Set NReWeightCascNucleon as the reweighting engine


  std::vector<neut::rew::NSyst_t> SystematicVector;

    SystematicVector.push_back(kCascTwkDial_TotalProb); 
    SystematicVector.push_back(kCascTwkDial_ElasticProb);
    SystematicVector.push_back(kCascTwkDial_SinglePiProb);
    SystematicVector.push_back(kCascTwkDial_DoublePiProb);
    //prepare the nucleon dials

  // Initialise the systematic
  for (std::vector<neut::rew::NSyst_t>::iterator it = SystematicVector.begin(); it != SystematicVector.end(); ++it) {
    rw.Systematics().Init(*it);
  }


  // Number of reweights
  const int nSysts = 4;
  double ReWeightVals[nSysts] = {0.};
  //set a bunch of dial options
  if (variation == 0) {                                                                                                                                          
    ReWeightVals[0] = 0.5;                                                                                                                                 
    ReWeightVals[1] = 0.5;
    ReWeightVals[2] = 0.5;
    ReWeightVals[3] = 0.5;

  }
  else if (variation == 1) {                                                                                                                                      
    ReWeightVals[0] = 0;                                                                                                                                 

  }
  else if (variation == 2) {                                                                                                                                      
    ReWeightVals[0] = 10;                                                                                                                                  

  }
  //total dial 1.44
  else if (variation == 3) {                                                                                                                                      
    ReWeightVals[0] = 2.0;                                                                                                                              

  }                                                                                                                                                                  


  //total dial 0.64
  else if (variation == 4) {
    ReWeightVals[0] = -2.0;
  }
  //elastic dial 1.44
  else if (variation == 5) {
    ReWeightVals[1] = 2.0;
  }
  //elastic dial 0.64
  else if (variation == 6) {
    ReWeightVals[1] = -2.0;
  }
  //spi dial 1.44
  else if (variation == 7) {
    ReWeightVals[2] = 2.0;
  }
  //spi dial 0.64
  else if (variation == 8) {
    ReWeightVals[2] = -2.0;

  }
  //dpi dial 1.44
  else if (variation == 9) {
    ReWeightVals[3] = 2.0;

  }
  //dpi dial 0.64
  else if (variation == 10) {
    ReWeightVals[3] = -2.0;

  }





  
  
  // Make some 1D distributions for each mode
  //set up a bunch of validation plots for each neut mode
  //obviously this is optional based on what you want to plot
  //would probably be a good idea to just make root trees of the outputs

  //NB All neut modes combined correspond to "neut mode 7"
  const int nNEUTmodes = 53;
  TH1D *Weights = new TH1D("Weights", "Weights", 100, 0, 3);
  TH1D *Weighted[nNEUTmodes];
  TH1D *UnWeighted[nNEUTmodes];
  TH1D *WeightedAllParticles[nNEUTmodes];
  TH1D *UnWeightedAllParticles[nNEUTmodes];
  TH1D *WeightedAngle[nNEUTmodes];
  TH1D *UnWeightedAngle[nNEUTmodes];
  TH1D *WeightedAngleAllParticles[nNEUTmodes];
  TH1D *UnWeightedAngleAllParticles[nNEUTmodes];

  TH1D *NeutronWeighted[nNEUTmodes];
  TH1D *NeutronUnWeighted[nNEUTmodes];
  TH1D *NeutronWeightedAllParticles[nNEUTmodes];
  TH1D *NeutronUnWeightedAllParticles[nNEUTmodes];
  TH1D *NeutronWeightedAngle[nNEUTmodes];
  TH1D *NeutronUnWeightedAngle[nNEUTmodes];
  TH1D *NeutronWeightedAngleAllParticles[nNEUTmodes];
  TH1D *NeutronUnWeightedAngleAllParticles[nNEUTmodes];

  TH1D *ProtonWeighted[nNEUTmodes];
  TH1D *ProtonUnWeighted[nNEUTmodes];
  TH1D *ProtonWeightedAllParticles[nNEUTmodes];
  TH1D *ProtonUnWeightedAllParticles[nNEUTmodes];
  TH1D *ProtonWeightedAngle[nNEUTmodes];
  TH1D *ProtonUnWeightedAngle[nNEUTmodes];
  TH1D *ProtonWeightedAngleAllParticles[nNEUTmodes];
  TH1D *ProtonUnWeightedAngleAllParticles[nNEUTmodes];


  TH1D *PionWeighted[nNEUTmodes];
  TH1D *PionUnWeighted[nNEUTmodes];
  TH1D *PionWeightedAllParticles[nNEUTmodes];
  TH1D *PionUnWeightedAllParticles[nNEUTmodes];
  TH1D *PionWeightedAngle[nNEUTmodes];
  TH1D *PionUnWeightedAngle[nNEUTmodes];
  TH1D *PionWeightedAngleAllParticles[nNEUTmodes];
  TH1D *PionUnWeightedAngleAllParticles[nNEUTmodes];


  TH1D *LeptonWeighted[nNEUTmodes];
  TH1D *LeptonUnWeighted[nNEUTmodes];
  TH1D *LeptonWeightedAllParticles[nNEUTmodes];
  TH1D *LeptonUnWeightedAllParticles[nNEUTmodes];
  TH1D *LeptonWeightedAngle[nNEUTmodes];
  TH1D *LeptonUnWeightedAngle[nNEUTmodes];
  TH1D *LeptonWeightedAngleAllParticles[nNEUTmodes];
  TH1D *LeptonUnWeightedAngleAllParticles[nNEUTmodes];


  TH1D *NucleonMultiplicityWeighted[nNEUTmodes];
  TH1D *NucleonMultiplicityUnWeighted[nNEUTmodes];
  TH1D *ProtonMultiplicityWeighted[nNEUTmodes];
  TH1D *ProtonMultiplicityUnWeighted[nNEUTmodes];
  TH1D *NeutronMultiplicityWeighted[nNEUTmodes];
  TH1D *NeutronMultiplicityUnWeighted[nNEUTmodes];
  TH1D *PionMultiplicityWeighted[nNEUTmodes];
  TH1D *PionMultiplicityUnWeighted[nNEUTmodes];
  //make some 2D plots as well
  TH2D *AngleVSEnergyWeighted[nNEUTmodes];
  TH2D *AngleVSEnergyUnWeighted[nNEUTmodes];
  TH2D *LeadingNucleonAngleVSEnergyWeighted[nNEUTmodes];
  TH2D *LeadingNucleonAngleVSEnergyUnWeighted[nNEUTmodes];

  
  //iterate over neutmode
  for (int i = 0; i < nNEUTmodes; ++i) {
    std::stringstream ss;
    ss << "_" << i;
    //start with the 2D plots
    AngleVSEnergyWeighted[i]= new TH2D( (std::string("AngleVSEnergyWeighted")+ss.str()).c_str(),"",10,0,2,10,-1,1);
    AngleVSEnergyUnWeighted[i]= new TH2D( (std::string("AngleVSEnergyUnWeighted")+ss.str()).c_str(),"",10,0,2,10,-1,1);
    LeadingNucleonAngleVSEnergyWeighted[i]= new TH2D( (std::string("LeadingNucleonAngleVSEnergyWeighted")+ss.str()).c_str(),"",10,0,2,10,-1,1);
    LeadingNucleonAngleVSEnergyUnWeighted[i]= new TH2D( (std::string("LeadingNucleonAngleVSEnergyUnWeighted")+ss.str()).c_str(),"",10,0,2,10,-1,1);
   
    
    //then all the 1D plots
    Weighted[i] = new TH1D((std::string("Weighted")+ss.str()).c_str(), (std::string("Weighted")+ss.str()+";E_{nucleon} (GeV);d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, 0, 2);
    UnWeighted[i] = new TH1D((std::string("UnWeighted")+ss.str()).c_str(), (std::string("UnWeighted")+ss.str()+";E_{nucleon} (GeV);d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, 0, 2);
    WeightedAllParticles[i] = new TH1D((std::string("WeightedAllParticles")+ss.str()).c_str(), (std::string("WeightedAllParticles")+ss.str()+";E_{nucleon} (GeV);d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, 0, 2);
    UnWeightedAllParticles[i] = new TH1D((std::string("UnWeightedAllParticles")+ss.str()).c_str(), (std::string("UnWeightedAllParticles")+ss.str()+";E_{nucleon} (GeV);d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, 0, 2);
    
    WeightedAngle[i] = new TH1D((std::string("WeightedAngle")+ss.str()).c_str(), (std::string("WeightedAngle")+ss.str()+";costheta_{nucleon} ;d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, -1, 1);
    UnWeightedAngle[i] = new TH1D((std::string("UnWeightedAngle")+ss.str()).c_str(), (std::string("UnWeightedAngle")+ss.str()+";costheta_{nucleon} ;d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, -1, 1);
    WeightedAngleAllParticles[i] = new TH1D((std::string("WeightedAngleAllParticles")+ss.str()).c_str(), (std::string("WeightedAngleAllParticles")+ss.str()+";costheta_{nucleon} ;d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, -1, 1);    
    UnWeightedAngleAllParticles[i] = new TH1D((std::string("UnWeightedAngleAllParticles")+ss.str()).c_str(), (std::string("UnWeightedAngleAllParticles")+ss.str()+";costheta_{nucleon} ;d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, -1, 1);

    NeutronWeighted[i] = new TH1D((std::string("NeutronWeighted")+ss.str()).c_str(), (std::string("NeutronWeighted")+ss.str()+";E_{Neutron} (GeV);d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, 0, 2);
    NeutronUnWeighted[i] = new TH1D((std::string("NeutronUnWeighted")+ss.str()).c_str(), (std::string("NeutronUnWeighted")+ss.str()+";E_{Neutron} (GeV);d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, 0, 2);
    NeutronWeightedAllParticles[i] = new TH1D((std::string("NeutronWeightedAllParticles")+ss.str()).c_str(), (std::string("NeutronWeightedAllParticles")+ss.str()+";E_{Neutron} (GeV);d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, 0, 2);
    NeutronUnWeightedAllParticles[i] = new TH1D((std::string("NeutronUnWeightedAllParticles")+ss.str()).c_str(), (std::string("NeutronUnWeightedAllParticles")+ss.str()+";E_{Neutron} (GeV);d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, 0, 2);

    NeutronWeightedAngle[i] = new TH1D((std::string("NeutronWeightedAngle")+ss.str()).c_str(), (std::string("NeutronWeightedAngle")+ss.str()+";costheta_{Neutron} ;d#sigma/dp (cm^2/nucleon/GeV)").c_str(),50, -1, 1);   
    NeutronUnWeightedAngle[i] = new TH1D((std::string("NeutronUnWeightedAngle")+ss.str()).c_str(), (std::string("NeutronUnWeightedAngle")+ss.str()+";costheta_{Neutron} ;d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, -1, 1);
    NeutronWeightedAngleAllParticles[i] = new TH1D((std::string("NeutronWeightedAngleAllParticles")+ss.str()).c_str(), (std::string("NeutronWeightedAngleAllParticles")+ss.str()+";costheta_{Neutron} ;d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, -1, 1);
    NeutronUnWeightedAngleAllParticles[i] = new TH1D((std::string("NeutronUnWeightedAngleAllParticles")+ss.str()).c_str(), (std::string("NeutronUnWeightedAngleAllParticles")+ss.str()+";costheta_{Neutron} ;d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, -1, 1);



    ProtonWeighted[i] = new TH1D((std::string("ProtonWeighted")+ss.str()).c_str(), (std::string("ProtonWeighted")+ss.str()+";E_{Proton} (GeV);d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, 0, 2);
    ProtonUnWeighted[i] = new TH1D((std::string("ProtonUnWeighted")+ss.str()).c_str(), (std::string("ProtonUnWeighted")+ss.str()+";E_{Proton} (GeV);d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50,0, 2);
    ProtonWeightedAllParticles[i] = new TH1D((std::string("ProtonWeightedAllParticles")+ss.str()).c_str(), (std::string("ProtonWeightedAllParticles")+ss.str()+";E_{Proton} (GeV);d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, 0, 2);
    ProtonUnWeightedAllParticles[i] = new TH1D((std::string("ProtonUnWeightedAllParticles")+ss.str()).c_str(), (std::string("ProtonUnWeightedAllParticles")+ss.str()+";E_{Proton} (GeV);d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, 0, 2);

    ProtonWeightedAngle[i] = new TH1D((std::string("ProtonWeightedAngle")+ss.str()).c_str(), (std::string("ProtonWeightedAngle")+ss.str()+";costheta_{Proton} ;d#sigma/dp (cm^2/nucleon/GeV)").c_str(),50, -1, 1);
    ProtonUnWeightedAngle[i] = new TH1D((std::string("ProtonUnWeightedAngle")+ss.str()).c_str(), (std::string("ProtonUnWeightedAngle")+ss.str()+";costheta_{Proton} ;d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, -1, 1);
    ProtonWeightedAngleAllParticles[i] = new TH1D((std::string("ProtonWeightedAngleAllParticles")+ss.str()).c_str(), (std::string("ProtonWeightedAngleAllParticles")+ss.str()+";costheta_{Proton} ;d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, -1, 1);
    ProtonUnWeightedAngleAllParticles[i] = new TH1D((std::string("ProtonUnWeightedAngleAllParticles")+ss.str()).c_str(), (std::string("ProtonUnWeightedAngleAllParticles")+ss.str()+";costheta_{Proton} ;d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, -1, 1);

    PionWeighted[i] = new TH1D((std::string("PionWeighted")+ss.str()).c_str(), (std::string("PionWeighted")+ss.str()+";E_{Pion} (GeV);d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, 0, 1);
    PionUnWeighted[i] = new TH1D((std::string("PionUnWeighted")+ss.str()).c_str(), (std::string("PionUnWeighted")+ss.str()+";E_{Pion} (GeV);d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, 0, 1);
    PionWeightedAllParticles[i] = new TH1D((std::string("PionWeightedAllParticles")+ss.str()).c_str(), (std::string("PionWeightedAllParticles")+ss.str()+";E_{Pion} (GeV);d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, 0, 1);
    PionUnWeightedAllParticles[i] = new TH1D((std::string("PionUnWeightedAllParticles")+ss.str()).c_str(), (std::string("PionUnWeightedAllParticles")+ss.str()+";E_{Pion} (GeV);d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, 0, 1);
    PionWeightedAngle[i] = new TH1D((std::string("PionWeightedAngle")+ss.str()).c_str(), (std::string("PionWeightedAngle")+ss.str()+";costheta_{Pion} ;d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, -1, 1);
    PionUnWeightedAngle[i] = new TH1D((std::string("PionUnWeightedAngle")+ss.str()).c_str(), (std::string("PionUnWeightedAngle")+ss.str()+";costheta_{Pion} ;d#sigma/dp (cm^2/nucleon/GeV)").c_str(),50, -1, 1);
    PionWeightedAngleAllParticles[i] = new TH1D((std::string("PionWeightedAngleAllParticles")+ss.str()).c_str(), (std::string("PionWeightedAngleAllParticles")+ss.str()+";costheta_{Pion} ;d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, -1, 1);
    PionUnWeightedAngleAllParticles[i] = new TH1D((std::string("PionUnWeightedAngleAllParticles")+ss.str()).c_str(), (std::string("PionUnWeightedAngleAllParticles")+ss.str()+";costheta_{Pion} ;d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, -1, 1);








    LeptonWeighted[i] = new TH1D((std::string("LeptonWeighted")+ss.str()).c_str(), (std::string("LeptonWeighted")+ss.str()+";E_{Lepton} (GeV);d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, 0, 3);
    LeptonUnWeighted[i] = new TH1D((std::string("LeptonUnWeighted")+ss.str()).c_str(), (std::string("LeptonUnWeighted")+ss.str()+";E_{Lepton} (GeV);d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50,0, 3);
    LeptonWeightedAllParticles[i] = new TH1D((std::string("LeptonWeightedAllParticles")+ss.str()).c_str(), (std::string("LeptonWeightedAllParticles")+ss.str()+";E_{Lepton} (GeV);d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, 0, 3);
    LeptonUnWeightedAllParticles[i] = new TH1D((std::string("LeptonUnWeightedAllParticles")+ss.str()).c_str(), (std::string("LeptonUnWeightedAllParticles")+ss.str()+";E_{Lepton} (GeV);d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, 0, 3);

    LeptonWeightedAngle[i] = new TH1D((std::string("LeptonWeightedAngle")+ss.str()).c_str(), (std::string("LeptonWeightedAngle")+ss.str()+";costheta_{Lepton} ;d#sigma/dp (cm^2/nucleon/GeV)").c_str(),50, -1, 1);
    LeptonUnWeightedAngle[i] = new TH1D((std::string("LeptonUnWeightedAngle")+ss.str()).c_str(), (std::string("LeptonUnWeightedAngle")+ss.str()+";costheta_{Lepton} ;d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, -1, 1);
    LeptonWeightedAngleAllParticles[i] = new TH1D((std::string("LeptonWeightedAngleAllParticles")+ss.str()).c_str(), (std::string("LeptonWeightedAngleAllParticles")+ss.str()+";costheta_{Lepton} ;d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, -1, 1);
    LeptonUnWeightedAngleAllParticles[i] = new TH1D((std::string("LeptonUnWeightedAngleAllParticles")+ss.str()).c_str(), (std::string("LeptonUnWeightedAngleAllParticles")+ss.str()+";costheta_{Lepton} ;d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, -1, 1);




    NucleonMultiplicityWeighted[i] = new TH1D((std::string("NucleonMultiplicityWeighted")+ss.str()).c_str(), (std::string("NucleonMultiplicityWeighted")+ss.str()+";Number of nucleons;d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, 0, 20);
    NucleonMultiplicityUnWeighted[i] = new TH1D((std::string("NucleonMultiplicityUnWeighted")+ss.str()).c_str(), (std::string("NucleonMultiplicityUnWeighted")+ss.str()+";Number of nucleons;d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, 0, 20);

    ProtonMultiplicityWeighted[i] = new TH1D((std::string("ProtonMultiplicityWeighted")+ss.str()).c_str(), (std::string("ProtonMultiplicityWeighted")+ss.str()+";Number of protons;d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, 0, 20);
    ProtonMultiplicityUnWeighted[i] = new TH1D((std::string("ProtonMultiplicityUnWeighted")+ss.str()).c_str(), (std::string("ProtonMultiplicityUnWeighted")+ss.str()+";Number of protons;d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, 0, 20);

    NeutronMultiplicityWeighted[i] = new TH1D((std::string("NeutronMultiplicityWeighted")+ss.str()).c_str(), (std::string("NeutronMultiplicityWeighted")+ss.str()+";Number of neutrons;d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, 0, 20);
    NeutronMultiplicityUnWeighted[i] = new TH1D((std::string("NeutronMultiplicityUnWeighted")+ss.str()).c_str(), (std::string("NeutronMultiplicityUnWeighted")+ss.str()+";Number of neutrons;d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, 0, 20);




    PionMultiplicityWeighted[i] = new TH1D((std::string("PionMultiplicityWeighted")+ss.str()).c_str(), (std::string("PionMultiplicityWeighted")+ss.str()+";Number of pions;d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, 0, 20);
    PionMultiplicityUnWeighted[i] = new TH1D((std::string("PionMultiplicityUnWeighted")+ss.str()).c_str(), (std::string("PionMultiplicityUnWeighted")+ss.str()+";Number of pions;d#sigma/dp (cm^2/nucleon/GeV)").c_str(), 50, 0, 20);


  }


  // Set the systematics
  int SystCounter = 0;
  for (std::vector<neut::rew::NSyst_t>::iterator it = SystematicVector.begin(); it != SystematicVector.end(); ++it, ++SystCounter) {
    rw.Systematics().Set(*it, ReWeightVals[SystCounter]);
  }




  // Reconfigure the reweight engine
  rw.Reconfigure();



  

  // Setup the output variables
  double Q2, EnuTrue, W, ppi, Eout,lastbiggestEout,pionlastbiggestEout,leptonlastbiggestEout,costheta,lastbiggestcostheta,pionlastbiggestcostheta,leptonlastbiggestcostheta,pinit, pn, costhlep, costhpinu, costhpimu, costhpin, weight, ScaleFactor, Unity = -999;
  int mode, bestmode,pionbestmode,leptonbestmode, leppdg, nupdg, pipdg, nucpdg,protonmultiplicity,neutronmultiplicity,nucleonmultiplicity,pionmultiplicity = -999;
  int npart, npi, nnuc, nbad = 0;
  bool isproton, isneutron;
  std::stringstream ss;
  ss << "var" << variation;
  //set and name output file
  TFile *OutputFile = new TFile((InputName+"_to_"+ss.str()+"sumw2working.root").c_str(), "RECREATE");
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

  // custom number of events
  if(nEvtsToProcess!=-1) nEvents=nEvtsToProcess;

  int PrintWidth = 1;

 
  Tree->GetEntry(0);

  

  //have to scale to the cross-section for mono-energetic
  //this is done with a scaling "unity"
  //this is not done with a modified weighting "scalefactor" which has been set to 1
  Unity = nvect->Totcrs * 1E-38 /double(Tree->GetEntries());                                                                                                                        
   ScaleFactor = 1;
   

 // And each event

  for (int j = 0; j < nEvents; ++j) {

    Q2 = -999;
    EnuTrue = -999;
    W = -999;
    ppi = -999;
    Eout = -999;
    lastbiggestEout = -999;
    costheta = -999;
    lastbiggestcostheta = -999;
    pionlastbiggestEout = -999;
    pionlastbiggestcostheta = -999;
    leptonlastbiggestEout = -999;
    leptonlastbiggestcostheta = -999;
    


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

    isproton = false;
    isneutron = false;
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

    weight = rw.CalcWeight();
    //calculate the weight - the magic happens hear

    if (weight == 1.0) nbad++;
    
    
    if (weight !=1.0)
      {    //std::cout << "weight (not equal to 1) = " << weight << std::endl;
      }
    

    TLorentzVector Pnu = nvect->PartInfo(0)->fP;
    TLorentzVector Pinit = nvect->PartInfo(1)->fP;
    TLorentzVector Pout;
    TLorentzVector Ppi;
    TLorentzVector Pn;

    npart = nvect->Npart();
    EnuTrue = Pnu.E()/1.E3;

    costhlep = cos(Pnu.Vect().Angle(Pout.Vect()));
    mode = nvect->Mode;
    nupdg = nvect->PartInfo(0)->fPID;
    leppdg = nvect->PartInfo(2)->fPID;
    Q2 = -1*(Pnu-Pout)*(Pnu-Pout)/1.E6;
    W = sqrt((Pnu+Pinit-Pout)*(Pnu+Pinit-Pout))/1.E3;
    pinit = Pinit.Vect().Mag();


    //over events
    lastbiggestEout = -999;
    lastbiggestcostheta = -999;
    pionlastbiggestEout = -999;
    pionlastbiggestcostheta = -999;
    leptonlastbiggestEout = -999;
    leptonlastbiggestcostheta = -999;

    protonmultiplicity = 0;
    neutronmultiplicity = 0;
    nucleonmultiplicity = 0;
    pionmultiplicity = 0;
    
    //produce a plot of just the weights
    Weights->Fill(weight);
    for (int k = 2; k < nvect->Npart(); ++k) {
      
      Pout = nvect->PartInfo(k)->fP;
      Eout = sqrt((Pout.Px())*(Pout.Px()) + (Pout.Py())*(Pout.Py()) +(Pout.Pz())*(Pout.Pz())) /1.E3;
      costheta = nvect->PartInfo(k)->fP.CosTheta();
      //All final state nucleons
      //NB All neut modes combined correspond to "neut mode 7"
      if( ((nvect->PartInfo(k))->fPID == 2212 ||(nvect->PartInfo(k))->fPID == 2112  ) && (nvect->PartInfo(k))->fIsAlive == 1)
	{
	  
	  nucleonmultiplicity++;
	  
	  UnWeightedAllParticles[mode]->Fill(Eout,ScaleFactor);                                                                                                                                                       // std::cout << "weight = " << weight << std::endl;                                                                                                      
	  WeightedAllParticles[mode]->Fill(Eout, weight*ScaleFactor);                                                                                                                          
	  UnWeightedAllParticles[7]->Fill(Eout,ScaleFactor);                                                                                                                         
	  WeightedAllParticles[7]->Fill(Eout, weight*ScaleFactor);

	  UnWeightedAngleAllParticles[mode]->Fill(costheta,ScaleFactor);
	  // std::cout << "weight = " << weight << std::endl;                                                                                                                            
          WeightedAngleAllParticles[mode]->Fill(costheta, weight*ScaleFactor);
          UnWeightedAngleAllParticles[7]->Fill(costheta,ScaleFactor);
          WeightedAngleAllParticles[7]->Fill(costheta, weight*ScaleFactor);
	  AngleVSEnergyWeighted[7]->Fill(Eout,costheta,weight*ScaleFactor);
	  AngleVSEnergyUnWeighted[7]->Fill(Eout,costheta,ScaleFactor); 
	  AngleVSEnergyWeighted[mode]->Fill(Eout,costheta,weight*ScaleFactor);
          AngleVSEnergyUnWeighted[mode]->Fill(Eout,costheta,ScaleFactor);



	  //of which neutrons
	  if( (nvect->PartInfo(k))->fPID == 2112)
	    {
	      
	      neutronmultiplicity++;

	      NeutronUnWeightedAllParticles[mode]->Fill(Eout,ScaleFactor);                                                                                                
	      // std::cout << "weight = " << weight << std::endl;                                                                                                                           
	      NeutronWeightedAllParticles[mode]->Fill(Eout, weight*ScaleFactor);
	      NeutronUnWeightedAllParticles[7]->Fill(Eout,ScaleFactor);
	      NeutronWeightedAllParticles[7]->Fill(Eout, weight*ScaleFactor);
	      
	      NeutronUnWeightedAngleAllParticles[mode]->Fill(costheta,ScaleFactor);
	      // std::cout << "weight = " << weight << std::endl;                                                                                                                          
	      NeutronWeightedAngleAllParticles[mode]->Fill(costheta, weight*ScaleFactor);
	      NeutronUnWeightedAngleAllParticles[7]->Fill(costheta,ScaleFactor);
	      NeutronWeightedAngleAllParticles[7]->Fill(costheta, weight*ScaleFactor);
	  
	  
	    }
	  //of which protons
	  if( (nvect->PartInfo(k))->fPID == 2212)
            {
	      protonmultiplicity++;
	      
	      ProtonUnWeightedAllParticles[mode]->Fill(Eout,ScaleFactor);
	      // std::cout << "weight = " << weight << std::endl;                                                                                                                          
	      ProtonWeightedAllParticles[mode]->Fill(Eout, weight*ScaleFactor);
	      ProtonUnWeightedAllParticles[7]->Fill(Eout,ScaleFactor);
	      ProtonWeightedAllParticles[7]->Fill(Eout, weight*ScaleFactor);

	      ProtonUnWeightedAngleAllParticles[mode]->Fill(costheta,ScaleFactor);
	      // std::cout << "weight = " << weight << std::endl;                                                                                                                           
	      ProtonWeightedAngleAllParticles[mode]->Fill(costheta, weight*ScaleFactor);
	      ProtonUnWeightedAngleAllParticles[7]->Fill(costheta,ScaleFactor);
	      ProtonWeightedAngleAllParticles[7]->Fill(costheta, weight*ScaleFactor);
	    }
	  //of which highest energy final state nucleon
	  if(lastbiggestEout<Eout)
	    {

	      bestmode = mode;
	      lastbiggestEout = Eout;
	      lastbiggestcostheta = costheta;
	      //of which protons
	      if( (nvect->PartInfo(k))->fPID == 2212)
		{


		  isproton = true;
		  isneutron = false;
		}
	      //of which neutrons
	      if( (nvect->PartInfo(k))->fPID == 2112)
		{

		  isproton = false;
		  isneutron = true;
		}
	      
	    }
	}
      //all final state pions
      //NB All neut modes combined correspond to "neut mode 7"
      if( (abs((nvect->PartInfo(k))->fPID) == 211 ||(nvect->PartInfo(k))->fPID == 111 ) && (nvect->PartInfo(k))->fIsAlive == 1)
	{
	  
	  pionmultiplicity++;
	  
	  PionUnWeightedAllParticles[mode]->Fill(Eout,ScaleFactor);                                                                                                                                  
	  // std::cout << "weight = " << weight << std::endl;                                                                                                                         
          PionWeightedAllParticles[mode]->Fill(Eout, weight*ScaleFactor);
          PionUnWeightedAllParticles[7]->Fill(Eout,ScaleFactor);
          PionWeightedAllParticles[7]->Fill(Eout, weight*ScaleFactor);
          PionUnWeightedAngleAllParticles[mode]->Fill(costheta,ScaleFactor);
          // std::cout << "weight = " << weight << std::endl;                                                                                                                               
          PionWeightedAngleAllParticles[mode]->Fill(costheta, weight*ScaleFactor);
          PionUnWeightedAngleAllParticles[7]->Fill(costheta,ScaleFactor);
          PionWeightedAngleAllParticles[7]->Fill(costheta, weight*ScaleFactor);
	  //of which highest energy pion
	  if(pionlastbiggestEout<Eout)
            {

              pionbestmode = mode;
              pionlastbiggestEout = Eout;
              pionlastbiggestcostheta = costheta;
              

            }



	}
      //all final state leptons
      //NB All neut modes combined correspond to "neut mode 7"
      if( (abs((nvect->PartInfo(k))->fPID) == 11 || abs((nvect->PartInfo(k))->fPID) == 13 ||abs((nvect->PartInfo(k))->fPID) == 15 ) && (nvect->PartInfo(k))->fIsAlive == 1)
        {

          LeptonUnWeightedAllParticles[mode]->Fill(Eout,ScaleFactor);
          // std::cout << "weight = " << weight << std::endl;                                                                                                                               
          LeptonWeightedAllParticles[mode]->Fill(Eout, weight*ScaleFactor);
          LeptonUnWeightedAllParticles[7]->Fill(Eout,ScaleFactor);
          LeptonWeightedAllParticles[7]->Fill(Eout, weight*ScaleFactor);
          LeptonUnWeightedAngleAllParticles[mode]->Fill(costheta,ScaleFactor);
          // std::cout << "weight = " << weight << std::endl;                                                                                                                               
          LeptonWeightedAngleAllParticles[mode]->Fill(costheta, weight*ScaleFactor);
          LeptonUnWeightedAngleAllParticles[7]->Fill(costheta,ScaleFactor);
          LeptonWeightedAngleAllParticles[7]->Fill(costheta, weight*ScaleFactor);
	  //highest energy leptons
          if(leptonlastbiggestEout<Eout)
            {

              leptonbestmode = mode;
              leptonlastbiggestEout = Eout;
              leptonlastbiggestcostheta = costheta;


            }


	
        }
    

      /*
	UnWeighted[mode]->Fill(Eout);
	// std::cout << "weight = " << weight << std::endl;                                                                                                                               
	Weighted[mode]->Fill(Eout, weight*ScaleFactor);
	UnWeighted[7]->Fill(Eout);
	Weighted[7]->Fill(Eout, weight*ScaleFactor);
      */

	
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
    //Fill, with weights, scale factor just equals 1 so can be ignored
    NucleonMultiplicityUnWeighted[mode]->Fill(nucleonmultiplicity,ScaleFactor);
    NucleonMultiplicityWeighted[mode]->Fill(nucleonmultiplicity,weight*ScaleFactor);
    NucleonMultiplicityUnWeighted[7]->Fill(nucleonmultiplicity,ScaleFactor);
    NucleonMultiplicityWeighted[7]->Fill(nucleonmultiplicity,weight*ScaleFactor);

    ProtonMultiplicityUnWeighted[mode]->Fill(protonmultiplicity,ScaleFactor);
    ProtonMultiplicityWeighted[mode]->Fill(protonmultiplicity,weight*ScaleFactor);
    ProtonMultiplicityUnWeighted[7]->Fill(protonmultiplicity,ScaleFactor);
    ProtonMultiplicityWeighted[7]->Fill(protonmultiplicity,weight*ScaleFactor);

    NeutronMultiplicityUnWeighted[mode]->Fill(neutronmultiplicity,ScaleFactor);
    NeutronMultiplicityWeighted[mode]->Fill(neutronmultiplicity,weight*ScaleFactor);
    NeutronMultiplicityUnWeighted[7]->Fill(neutronmultiplicity,ScaleFactor);
    NeutronMultiplicityWeighted[7]->Fill(neutronmultiplicity,weight*ScaleFactor);

    PionMultiplicityUnWeighted[mode]->Fill(pionmultiplicity,ScaleFactor);
    PionMultiplicityWeighted[mode]->Fill(pionmultiplicity,weight*ScaleFactor);
    PionMultiplicityUnWeighted[7]->Fill(pionmultiplicity,ScaleFactor);
    PionMultiplicityWeighted[7]->Fill(pionmultiplicity,weight*ScaleFactor);



    if(lastbiggestEout >= 0)
      {
	UnWeighted[bestmode]->Fill(lastbiggestEout,ScaleFactor);
	// std::cout << "weight = " << weight << std::endl;                                                                                                                               
	Weighted[bestmode]->Fill(lastbiggestEout, weight*ScaleFactor);
	UnWeighted[7]->Fill(lastbiggestEout,ScaleFactor);
	Weighted[7]->Fill(lastbiggestEout, weight*ScaleFactor);

	UnWeightedAngle[bestmode]->Fill(lastbiggestcostheta,ScaleFactor);
	// std::cout << "weight = " << weight << std::endl;                                                                                                                                    
	WeightedAngle[bestmode]->Fill(lastbiggestcostheta, weight*ScaleFactor);
	UnWeightedAngle[7]->Fill(lastbiggestcostheta,ScaleFactor);
	WeightedAngle[7]->Fill(lastbiggestcostheta, weight*ScaleFactor);
	
	LeadingNucleonAngleVSEnergyWeighted[7]->Fill(lastbiggestEout,lastbiggestcostheta,weight*ScaleFactor);
	LeadingNucleonAngleVSEnergyUnWeighted[7]->Fill(lastbiggestEout,lastbiggestcostheta,ScaleFactor);
	LeadingNucleonAngleVSEnergyWeighted[bestmode]->Fill(lastbiggestEout,lastbiggestcostheta,weight*ScaleFactor);
	LeadingNucleonAngleVSEnergyUnWeighted[bestmode]->Fill(lastbiggestEout,lastbiggestcostheta,ScaleFactor);

	//Handle all the highest energy particle only plots
	if( isneutron == true)
	  {

	    NeutronUnWeighted[bestmode]->Fill(lastbiggestEout,ScaleFactor);
	    // std::cout << "weight = " << weight << std::endl;                                                                                                                                 
	    NeutronWeighted[bestmode]->Fill(lastbiggestEout, weight*ScaleFactor);
	    NeutronUnWeighted[7]->Fill(lastbiggestEout,ScaleFactor);
	    NeutronWeighted[7]->Fill(lastbiggestEout, weight*ScaleFactor);
	
	    NeutronUnWeightedAngle[bestmode]->Fill(lastbiggestcostheta,ScaleFactor);
	    // std::cout << "weight = " << weight << std::endl;                                                                                                                                 
	    NeutronWeightedAngle[bestmode]->Fill(lastbiggestcostheta, weight*ScaleFactor);
	    NeutronUnWeightedAngle[7]->Fill(lastbiggestcostheta,ScaleFactor);
	    NeutronWeightedAngle[7]->Fill(lastbiggestcostheta, weight*ScaleFactor);
	
	  }
      
	if( isproton == true)
	  {


	    ProtonUnWeighted[bestmode]->Fill(lastbiggestEout,ScaleFactor);
	    // std::cout << "weight = " << weight << std::endl;                                                                                                                                 
	    ProtonWeighted[bestmode]->Fill(lastbiggestEout, weight*ScaleFactor);
	    ProtonUnWeighted[7]->Fill(lastbiggestEout,ScaleFactor);
	    ProtonWeighted[7]->Fill(lastbiggestEout, weight*ScaleFactor);
	    ProtonUnWeightedAngle[bestmode]->Fill(lastbiggestcostheta,ScaleFactor);
	    // std::cout << "weight = " << weight << std::endl;                                                                                                                                 
	    ProtonWeightedAngle[bestmode]->Fill(lastbiggestcostheta, weight*ScaleFactor);
	    ProtonUnWeightedAngle[7]->Fill(lastbiggestcostheta,ScaleFactor);
	    ProtonWeightedAngle[7]->Fill(lastbiggestcostheta, weight*ScaleFactor);

	  }

      }

    if(pionlastbiggestEout >= 0)
      {
	PionUnWeighted[pionbestmode]->Fill(pionlastbiggestEout,ScaleFactor);
	// std::cout << "weight = " << weight << std::endl;                                                                                                                               
	PionWeighted[pionbestmode]->Fill(pionlastbiggestEout, weight*ScaleFactor);
	PionUnWeighted[7]->Fill(pionlastbiggestEout,ScaleFactor);
	PionWeighted[7]->Fill(pionlastbiggestEout, weight*ScaleFactor);
	PionUnWeightedAngle[pionbestmode]->Fill(pionlastbiggestcostheta,ScaleFactor);
	// std::cout << "weight = " << weight << std::endl;                                                                                                                                
	PionWeightedAngle[pionbestmode]->Fill(pionlastbiggestcostheta, weight*ScaleFactor);
	PionUnWeightedAngle[7]->Fill(pionlastbiggestcostheta,ScaleFactor);
	PionWeightedAngle[7]->Fill(pionlastbiggestcostheta, weight*ScaleFactor);
      }



    if(leptonlastbiggestEout >= 0)
      {
        LeptonUnWeighted[leptonbestmode]->Fill(leptonlastbiggestEout,ScaleFactor);
        // std::cout << "weight = " << weight << std::endl;                                                                                                                                  
        LeptonWeighted[leptonbestmode]->Fill(leptonlastbiggestEout, weight*ScaleFactor);
        LeptonUnWeighted[7]->Fill(leptonlastbiggestEout,ScaleFactor);
        LeptonWeighted[7]->Fill(leptonlastbiggestEout, weight*ScaleFactor);
        LeptonUnWeightedAngle[leptonbestmode]->Fill(leptonlastbiggestcostheta,ScaleFactor);
        // std::cout << "weight = " << weight << std::endl;                                                                                                                                  
        LeptonWeightedAngle[leptonbestmode]->Fill(leptonlastbiggestcostheta, weight*ScaleFactor);
        LeptonUnWeightedAngle[7]->Fill(leptonlastbiggestcostheta,ScaleFactor);
        LeptonWeightedAngle[7]->Fill(leptonlastbiggestcostheta, weight*ScaleFactor);
      }
  

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
    //requires all have same number of bins!
    
    /*for (int i = 0; i < Weighted[j]->GetNbinsX()+1; ++i) {
      double content1 = Weighted[j]->GetBinContent(i+1);
      double error1 = sqrt(content1);
      Weighted[j]->SetBinError(i+1, error1);
      double content2 = UnWeighted[j]->GetBinContent(i+1);
      double error2 = sqrt(content2);
      UnWeighted[j]->SetBinError(i+1, error2);

      double content3 = WeightedAllParticles[j]->GetBinContent(i+1);
      double error3 = sqrt(content3);
      WeightedAllParticles[j]->SetBinError(i+1, error3);
      double content4 = UnWeightedAllParticles[j]->GetBinContent(i+1);
      double error4 = sqrt(content4);
      UnWeightedAllParticles[j]->SetBinError(i+1, error4);

      double content5 = WeightedAngle[j]->GetBinContent(i+1);
      double error5 = sqrt(content5);
      WeightedAngle[j]->SetBinError(i+1, error5);
      double content6 = UnWeighted[j]->GetBinContent(i+1);
      double error6 = sqrt(content6);
      UnWeightedAngle[j]->SetBinError(i+1, error6);

      double content7 = WeightedAngleAllParticles[j]->GetBinContent(i+1);
      double error7 = sqrt(content7);
      WeightedAngleAllParticles[j]->SetBinError(i+1, error7);
      double content8 = UnWeightedAngleAllParticles[j]->GetBinContent(i+1);
      double error8 = sqrt(content8);
      UnWeightedAngleAllParticles[j]->SetBinError(i+1, error8);


      double content9 = NeutronWeighted[j]->GetBinContent(i+1);
      double error9 = sqrt(content9);
      NeutronWeighted[j]->SetBinError(i+1, error9);
      double content10 = NeutronUnWeighted[j]->GetBinContent(i+1);
      double error10 = sqrt(content10);
      NeutronUnWeighted[j]->SetBinError(i+1, error10);

      double content11 = NeutronWeightedAllParticles[j]->GetBinContent(i+1);
      double error11 = sqrt(content11);
      NeutronWeightedAllParticles[j]->SetBinError(i+1, error11);
      double content12 = NeutronUnWeightedAllParticles[j]->GetBinContent(i+1);
      double error12 = sqrt(content12);
      NeutronUnWeightedAllParticles[j]->SetBinError(i+1, error12);

      double content13 = NeutronWeightedAngle[j]->GetBinContent(i+1);
      double error13 = sqrt(content13);
      NeutronWeightedAngle[j]->SetBinError(i+1, error13);
      double content14 = NeutronUnWeighted[j]->GetBinContent(i+1);
      double error14 = sqrt(content14);
      NeutronUnWeightedAngle[j]->SetBinError(i+1, error14);

      double content15 = NeutronWeightedAngleAllParticles[j]->GetBinContent(i+1);
      double error15 = sqrt(content15);
      NeutronWeightedAngleAllParticles[j]->SetBinError(i+1, error15);
      double content16 = NeutronUnWeightedAngleAllParticles[j]->GetBinContent(i+1);
      double error16 = sqrt(content16);
      NeutronUnWeightedAngleAllParticles[j]->SetBinError(i+1, error16);


      
      double content17 = ProtonWeighted[j]->GetBinContent(i+1);
      double error17 = sqrt(content17);
      ProtonWeighted[j]->SetBinError(i+1, error17);
      double content18 = ProtonUnWeighted[j]->GetBinContent(i+1);
      double error18 = sqrt(content18);
      ProtonUnWeighted[j]->SetBinError(i+1, error18);

      double content19 = ProtonWeightedAllParticles[j]->GetBinContent(i+1);
      double error19 = sqrt(content19);
      ProtonWeightedAllParticles[j]->SetBinError(i+1, error19);
      double content20 = ProtonUnWeightedAllParticles[j]->GetBinContent(i+1);
      double error20 = sqrt(content20);
      ProtonUnWeightedAllParticles[j]->SetBinError(i+1, error20);

      double content21 = ProtonWeightedAngle[j]->GetBinContent(i+1);
      double error21 = sqrt(content21);
      ProtonWeightedAngle[j]->SetBinError(i+1, error21);
      double content22 = ProtonUnWeighted[j]->GetBinContent(i+1);
      double error22 = sqrt(content22);
      ProtonUnWeightedAngle[j]->SetBinError(i+1, error22);

      double content23 = ProtonWeightedAngleAllParticles[j]->GetBinContent(i+1);
      double error23 = sqrt(content23);
      ProtonWeightedAngleAllParticles[j]->SetBinError(i+1, error23);
      double content24 = ProtonUnWeightedAngleAllParticles[j]->GetBinContent(i+1);
      double error24 = sqrt(content24);
      ProtonUnWeightedAngleAllParticles[j]->SetBinError(i+1, error24);

      double content25 = PionWeighted[j]->GetBinContent(i+1);
      double error25 = sqrt(content25);
      PionWeighted[j]->SetBinError(i+1, error25);
      double content26 = PionUnWeighted[j]->GetBinContent(i+1);
      double error26 = sqrt(content26);
      PionUnWeighted[j]->SetBinError(i+1, error26);

      double content27 = PionWeightedAllParticles[j]->GetBinContent(i+1);
      double error27 = sqrt(content27);
      PionWeightedAllParticles[j]->SetBinError(i+1, error27);
      double content28 = PionUnWeightedAllParticles[j]->GetBinContent(i+1);
      double error28 = sqrt(content28);
      PionUnWeightedAllParticles[j]->SetBinError(i+1, error28);

      double content29 = PionWeightedAngle[j]->GetBinContent(i+1);
      double error29 = sqrt(content29);
      PionWeightedAngle[j]->SetBinError(i+1, error29);
      double content30 = PionUnWeighted[j]->GetBinContent(i+1);
      double error30 = sqrt(content30);
      PionUnWeightedAngle[j]->SetBinError(i+1, error30);

      double content31 = PionWeightedAngleAllParticles[j]->GetBinContent(i+1);
      double error31 = sqrt(content31);
      PionWeightedAngleAllParticles[j]->SetBinError(i+1, error31);
      double content32 = PionUnWeightedAngleAllParticles[j]->GetBinContent(i+1);
      double error32 = sqrt(content32);
      PionUnWeightedAngleAllParticles[j]->SetBinError(i+1, error32);




      double content33 = LeptonWeighted[j]->GetBinContent(i+1);
      double error33 = sqrt(content33);
      LeptonWeighted[j]->SetBinError(i+1, error33);
      double content34 = LeptonUnWeighted[j]->GetBinContent(i+1);
      double error34 = sqrt(content34);
      LeptonUnWeighted[j]->SetBinError(i+1, error34);

      double content35 = LeptonWeightedAllParticles[j]->GetBinContent(i+1);
      double error35 = sqrt(content35);
      LeptonWeightedAllParticles[j]->SetBinError(i+1, error35);
      double content36 = LeptonUnWeightedAllParticles[j]->GetBinContent(i+1);
      double error36 = sqrt(content36);
      LeptonUnWeightedAllParticles[j]->SetBinError(i+1, error36);

      double content37 = LeptonWeightedAngle[j]->GetBinContent(i+1);
      double error37 = sqrt(content37);
      LeptonWeightedAngle[j]->SetBinError(i+1, error37);
      double content38 = LeptonUnWeighted[j]->GetBinContent(i+1);
      double error38 = sqrt(content38);
      LeptonUnWeightedAngle[j]->SetBinError(i+1, error38);

      double content39 = LeptonWeightedAngleAllParticles[j]->GetBinContent(i+1);
      double error39 = sqrt(content39);
      LeptonWeightedAngleAllParticles[j]->SetBinError(i+1, error39);
      double content40 = LeptonUnWeightedAngleAllParticles[j]->GetBinContent(i+1);
      double error40 = sqrt(content40);
      LeptonUnWeightedAngleAllParticles[j]->SetBinError(i+1, error40);

      






      double content41 = NucleonMultiplicityWeighted[j]->GetBinContent(i+1);
      double error41 = sqrt(content41);
      NucleonMultiplicityWeighted[j]->SetBinError(i+1, error41);
      double content42 = NucleonMultiplicityUnWeighted[j]->GetBinContent(i+1);
      double error42 = sqrt(content42);
      NucleonMultiplicityUnWeighted[j]->SetBinError(i+1, error42);
      
      double content43 = ProtonMultiplicityWeighted[j]->GetBinContent(i+1);
      double error43 = sqrt(content43);
      ProtonMultiplicityWeighted[j]->SetBinError(i+1, error43);
      double content44 = ProtonMultiplicityUnWeighted[j]->GetBinContent(i+1);
      double error44 = sqrt(content44);
      ProtonMultiplicityUnWeighted[j]->SetBinError(i+1, error44);

      double content45 = NeutronMultiplicityWeighted[j]->GetBinContent(i+1);
      double error45 = sqrt(content45);
      NeutronMultiplicityWeighted[j]->SetBinError(i+1, error45);
      double content46 = NeutronMultiplicityUnWeighted[j]->GetBinContent(i+1);
      double error46 = sqrt(content46);
      NeutronMultiplicityUnWeighted[j]->SetBinError(i+1, error46);

      double content47 = PionMultiplicityWeighted[j]->GetBinContent(i+1);
      double error47 = sqrt(content47);
      PionMultiplicityWeighted[j]->SetBinError(i+1, error47);
      double content48 = PionMultiplicityUnWeighted[j]->GetBinContent(i+1);
      double error48 = sqrt(content48);
      PionMultiplicityUnWeighted[j]->SetBinError(i+1, error48);
}

    */

    //perform scaling to cross-section
    //NB "Unity" != 1
    //only write if there are more than 0 events
    Weighted[j]->Scale(Unity, "width");
    UnWeighted[j]->Scale(Unity, "width");
    if (Weighted[j]->Integral() > 0) Weighted[j]->Write();
    if (UnWeighted[j]->Integral() > 0) UnWeighted[j]->Write();

    WeightedAllParticles[j]->Scale(Unity, "width");
    UnWeightedAllParticles[j]->Scale(Unity, "width");
    if (WeightedAllParticles[j]->Integral() > 0) WeightedAllParticles[j]->Write();
    if (UnWeightedAllParticles[j]->Integral() > 0) UnWeightedAllParticles[j]->Write();

    WeightedAngle[j]->Scale(Unity, "width");
    UnWeightedAngle[j]->Scale(Unity, "width");
    if (WeightedAngle[j]->Integral() > 0) WeightedAngle[j]->Write();
    if (UnWeightedAngle[j]->Integral() > 0) UnWeightedAngle[j]->Write();

    WeightedAngleAllParticles[j]->Scale(Unity, "width");
    UnWeightedAngleAllParticles[j]->Scale(Unity, "width");
    if (WeightedAngleAllParticles[j]->Integral() > 0) WeightedAngleAllParticles[j]->Write();
    if (UnWeightedAngleAllParticles[j]->Integral() > 0) UnWeightedAngleAllParticles[j]->Write();


    AngleVSEnergyWeighted[j]->Scale(Unity, "width");
    AngleVSEnergyUnWeighted[j]->Scale(Unity, "width");
    if (AngleVSEnergyWeighted[j]->Integral() > 0) AngleVSEnergyWeighted[j]->Write();
    if (AngleVSEnergyUnWeighted[j]->Integral() > 0) AngleVSEnergyUnWeighted[j]->Write();

    LeadingNucleonAngleVSEnergyWeighted[j]->Scale(Unity, "width");
    LeadingNucleonAngleVSEnergyUnWeighted[j]->Scale(Unity, "width");
    if (LeadingNucleonAngleVSEnergyWeighted[j]->Integral() > 0) LeadingNucleonAngleVSEnergyWeighted[j]->Write();
    if (LeadingNucleonAngleVSEnergyUnWeighted[j]->Integral() > 0) LeadingNucleonAngleVSEnergyUnWeighted[j]->Write();




    NeutronWeighted[j]->Scale(Unity, "width");
    NeutronUnWeighted[j]->Scale(Unity, "width");
    if (NeutronWeighted[j]->Integral() > 0) NeutronWeighted[j]->Write();
    if (NeutronUnWeighted[j]->Integral() > 0) NeutronUnWeighted[j]->Write();

    NeutronWeightedAllParticles[j]->Scale(Unity, "width");
    NeutronUnWeightedAllParticles[j]->Scale(Unity, "width");
    if (NeutronWeightedAllParticles[j]->Integral() > 0) NeutronWeightedAllParticles[j]->Write();
    if (NeutronUnWeightedAllParticles[j]->Integral() > 0) NeutronUnWeightedAllParticles[j]->Write();

    NeutronWeightedAngle[j]->Scale(Unity, "width");
    NeutronUnWeightedAngle[j]->Scale(Unity, "width");
    if (NeutronWeightedAngle[j]->Integral() > 0) NeutronWeightedAngle[j]->Write();
    if (NeutronUnWeightedAngle[j]->Integral() > 0) NeutronUnWeightedAngle[j]->Write();

    NeutronWeightedAngleAllParticles[j]->Scale(Unity, "width");
    NeutronUnWeightedAngleAllParticles[j]->Scale(Unity, "width");
    if (NeutronWeightedAngleAllParticles[j]->Integral() > 0) NeutronWeightedAngleAllParticles[j]->Write();
    if (NeutronUnWeightedAngleAllParticles[j]->Integral() > 0) NeutronUnWeightedAngleAllParticles[j]->Write();
  



    //    std::cout << "protonweighted integral = " << ProtonWeighted[j]->Integral() << std::endl;
    ProtonWeighted[j]->Scale(Unity, "width");
    ProtonUnWeighted[j]->Scale(Unity, "width");
    if (ProtonWeighted[j]->Integral() > 0) ProtonWeighted[j]->Write();
    if (ProtonUnWeighted[j]->Integral() > 0) ProtonUnWeighted[j]->Write();

    ProtonWeightedAllParticles[j]->Scale(Unity, "width");
    ProtonUnWeightedAllParticles[j]->Scale(Unity, "width");
    if (ProtonWeightedAllParticles[j]->Integral() > 0) ProtonWeightedAllParticles[j]->Write();
    if (ProtonUnWeightedAllParticles[j]->Integral() > 0) ProtonUnWeightedAllParticles[j]->Write();

    ProtonWeightedAngle[j]->Scale(Unity, "width");
    ProtonUnWeightedAngle[j]->Scale(Unity, "width");
    if (ProtonWeightedAngle[j]->Integral() > 0) ProtonWeightedAngle[j]->Write();
    if (ProtonUnWeightedAngle[j]->Integral() > 0) ProtonUnWeightedAngle[j]->Write();

    ProtonWeightedAngleAllParticles[j]->Scale(Unity, "width");
    ProtonUnWeightedAngleAllParticles[j]->Scale(Unity, "width");
    if (ProtonWeightedAngleAllParticles[j]->Integral() > 0) ProtonWeightedAngleAllParticles[j]->Write();
    if (ProtonUnWeightedAngleAllParticles[j]->Integral() > 0) ProtonUnWeightedAngleAllParticles[j]->Write();

    PionWeighted[j]->Scale(Unity, "width");
    PionUnWeighted[j]->Scale(Unity, "width");
    if (PionWeighted[j]->Integral() > 0) PionWeighted[j]->Write();
    if (PionUnWeighted[j]->Integral() > 0) PionUnWeighted[j]->Write();

    PionWeightedAllParticles[j]->Scale(Unity, "width");
    PionUnWeightedAllParticles[j]->Scale(Unity, "width");
    if (PionWeightedAllParticles[j]->Integral() > 0) PionWeightedAllParticles[j]->Write();
    if (PionUnWeightedAllParticles[j]->Integral() > 0) PionUnWeightedAllParticles[j]->Write();

    PionWeightedAngle[j]->Scale(Unity, "width");
    PionUnWeightedAngle[j]->Scale(Unity, "width");
    if (PionWeightedAngle[j]->Integral() > 0) PionWeightedAngle[j]->Write();
    if (PionUnWeightedAngle[j]->Integral() > 0) PionUnWeightedAngle[j]->Write();

    PionWeightedAngleAllParticles[j]->Scale(Unity, "width");
    PionUnWeightedAngleAllParticles[j]->Scale(Unity, "width");
    if (PionWeightedAngleAllParticles[j]->Integral() > 0) PionWeightedAngleAllParticles[j]->Write();
    if (PionUnWeightedAngleAllParticles[j]->Integral() > 0) PionUnWeightedAngleAllParticles[j]->Write();



    LeptonWeighted[j]->Scale(Unity, "width");
    LeptonUnWeighted[j]->Scale(Unity, "width");
    if (LeptonWeighted[j]->Integral() > 0) LeptonWeighted[j]->Write();
    if (LeptonUnWeighted[j]->Integral() > 0) LeptonUnWeighted[j]->Write();

    LeptonWeightedAllParticles[j]->Scale(Unity, "width");
    LeptonUnWeightedAllParticles[j]->Scale(Unity, "width");
    if (LeptonWeightedAllParticles[j]->Integral() > 0) LeptonWeightedAllParticles[j]->Write();
    if (LeptonUnWeightedAllParticles[j]->Integral() > 0) LeptonUnWeightedAllParticles[j]->Write();

    LeptonWeightedAngle[j]->Scale(Unity, "width");
    LeptonUnWeightedAngle[j]->Scale(Unity, "width");
    if (LeptonWeightedAngle[j]->Integral() > 0) LeptonWeightedAngle[j]->Write();
    if (LeptonUnWeightedAngle[j]->Integral() > 0) LeptonUnWeightedAngle[j]->Write();

    LeptonWeightedAngleAllParticles[j]->Scale(Unity, "width");
    LeptonUnWeightedAngleAllParticles[j]->Scale(Unity, "width");
    if (LeptonWeightedAngleAllParticles[j]->Integral() > 0) LeptonWeightedAngleAllParticles[j]->Write();
    if (LeptonUnWeightedAngleAllParticles[j]->Integral() > 0) LeptonUnWeightedAngleAllParticles[j]->Write();

    NucleonMultiplicityWeighted[j]->Scale(Unity, "width");
    NucleonMultiplicityUnWeighted[j]->Scale(Unity, "width");
    if (NucleonMultiplicityWeighted[j]->Integral() > 0) NucleonMultiplicityWeighted[j]->Write();
    if (NucleonMultiplicityUnWeighted[j]->Integral() > 0) NucleonMultiplicityUnWeighted[j]->Write();

    ProtonMultiplicityWeighted[j]->Scale(Unity, "width");
    ProtonMultiplicityUnWeighted[j]->Scale(Unity, "width");
    if (ProtonMultiplicityWeighted[j]->Integral() > 0) ProtonMultiplicityWeighted[j]->Write();
    if (ProtonMultiplicityUnWeighted[j]->Integral() > 0) ProtonMultiplicityUnWeighted[j]->Write();

    NeutronMultiplicityWeighted[j]->Scale(Unity, "width");
    NeutronMultiplicityUnWeighted[j]->Scale(Unity, "width");
    if (NeutronMultiplicityWeighted[j]->Integral() > 0) NeutronMultiplicityWeighted[j]->Write();
    if (NeutronMultiplicityUnWeighted[j]->Integral() > 0) NeutronMultiplicityUnWeighted[j]->Write();

    PionMultiplicityWeighted[j]->Scale(Unity, "width");
    PionMultiplicityUnWeighted[j]->Scale(Unity, "width");
    if (PionMultiplicityWeighted[j]->Integral() > 0) PionMultiplicityWeighted[j]->Write();
    if (PionMultiplicityUnWeighted[j]->Integral() > 0) PionMultiplicityUnWeighted[j]->Write();



  }

  Weights->Write();


  return 0;
}
