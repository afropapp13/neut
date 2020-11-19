#include <iostream>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <cstdlib>
#include <climits>
#include <cerrno>

#include "TFile.h"
#include "TTree.h"

#include "NEUTCIncludes.hxx"
#include "nuclscat_output_blocks.hxx"

extern "C" {
  void _gfortran_set_args(int argc, char const * argv[]);
}

//Intermediate global variables are defined here.
namespace {
  //Maximum interaction vertex radius.
  Float_t RMax;
  //PDG of particle to 'track'.
  Int_t IdNucl;
  //Total momentum of initial particle.
  Float_t PTot;

  TFile* OutputFile = NULL;
  TTree* OutputTree = NULL;

  char const * OutputFileName;

  enum STR2INT_ERROR { STR2INT_SUCCESS, STR2INT_OVERFLOW, STR2INT_UNDERFLOW,
		       STR2INT_INCONVERTIBLE };
  //check for error in length of string
  STR2INT_ERROR str2int (int &i, char const *s, int base = 0) {
    char *end;
    long  l;
    errno = 0;
    l = strtol(s, &end, base);
    if ((errno == ERANGE && l == LONG_MAX) || l > INT_MAX) {
      return STR2INT_OVERFLOW;
    }
    if ((errno == ERANGE && l == LONG_MIN) || l < INT_MIN) {
      return STR2INT_UNDERFLOW;
    }
    if (*s == '\0' || *end != '\0') {
      return STR2INT_INCONVERTIBLE;
    }
    i = l;
    return STR2INT_SUCCESS;
  }
}
//Throw until you get a random position within RMax
void SetRandPosition(Float_t (&PosV)[3]){
  Float_t rad = 0;
  do {
    PosV[1] = -RMax + rlu_(0)*RMax*2;
    PosV[2] = -RMax + rlu_(0)*RMax*2;

    rad = sqrt(PosV[1]*PosV[1]
	       + PosV[2]*PosV[2]);
  } while (RMax < rad); // Keep trying until you get one that is okay.
  PosV[0] = -sqrt(RMax*RMax - rad*rad) + 0.0001;
}
//Get the following arguments: card file, output file name, pdg 2212 or 2112
int HandleCLIArgs(int argc, char const * argv[]){
  _gfortran_set_args(argc,argv);
  if(argc < 4){
    std::cerr << "[ERROR]: Fewer than 3 CLI args recieved." << std::endl;
    std::cerr << "I want ./nuclscat.exe CARD_FILE.card OUTPUT_FILE.root PDG_CODE(2112, 2212)" << std::endl;
    return 1;
  }
  OutputFileName = argv[2];

  try{

    if(str2int(IdNucl, argv[3]) != STR2INT_SUCCESS){
      throw std::invalid_argument(
				  (std::string("Expected proton (2212) or neutron (2112) PDG code, "
					       "but found: ") + argv[3]).c_str());
    }
    if((IdNucl != 2212) && (IdNucl != 2112)){
      throw std::invalid_argument(
				  (std::string("Expected proton (2212) or neutron (2112) PDG code, "
					       "but found: ") + argv[3]).c_str());
    }
  } catch(std::exception const &e){
    std::cerr << "[ERROR]: " << e.what() << std::endl;
    return 2;
  }

  if(argc > 4){
    std::cerr << "[WARN]: Ignoring " << (argc-3) << " CLI arguments:"
	      << std::endl;
    for(int i = 4; i < argc; ++i){
      std::cerr << "\t[Arg #"<< i << "]: " << argv[i] << std::endl;
    }
  }
  return 0;
}

int main(int argc, char const * argv[]){
  if(HandleCLIArgs(argc,argv)){
    return 1;
  }

  necard_();
  std::cout << "[INFO]: Read NECard." << std::endl;
  necardev_();
  std::cout << "[INFO]: Read NECardEv." << std::endl;
  DebugVals();

  if(!nevccard->nectecvt){
    std::cout << "[INFO]: Generated 0 events." << std::endl;
    return 0;
  }

  std::cout << "Trying to track particle: " << IdNucl
	    << std::endl;

  if(!(OutputFile = TFile::Open(OutputFileName, "RECREATE"))){
    return 1;
  }
  OutputTree = new TTree("h1","Nucleon Scattering Vector");

  // if( !Add_NEUTEventCommon_Branches(OutputTree) ||
  //     !Add_VectorCommon_Branches(OutputTree) ||
  //     !Add_FSIHistoryCommon_Branches(OutputTree) ||
  //     !Add_NucleusTargetCommon_Branches(OutputTree) ||
  //     !Add_NucleonFSIHistoryCommon_Branches(OutputTree) ){
  //   std::cerr << "[ERROR]: Failed to add output branches." << std::endl;
  //   return 2;
  // }
  Float_t Abspvc[100];
  Int_t Ibound = 1;
  ///Add direct common block 'add branches here'
  //#include "vcworkC.h"

  // OutputTree->SetBranchAddress("Npvc",&vcwork->nvc);
  // OutputTree->SetBranchAddress("Ipvc",&vcwork->ipvc);
  // OutputTree->SetBranchAddress("Ichvc",&vcwork->icrnvc);
  // OutputTree->SetBranchAddress("Iorgvc", &vcwork->iorgvc);
  // OutputTree->SetBranchAddress("Iflvc", &vcwork->iflgvc);
  // OutputTree->SetBranchAddress("Abspvc", Abspvc);//
  // OutputTree->SetBranchAddress("Pvc", &vcwork->pvc);

  // OutputTree->SetBranchAddress("Nvert", &fsihist->nvert);
  // OutputTree->SetBranchAddress("Posvert", &fsihist->posvert);
  // OutputTree->SetBranchAddress("Iflgvert", &fsihist->iflgvert);
  // OutputTree->SetBranchAddress("Nvcvert", &fsihist->nvcvert);
  // OutputTree->SetBranchAddress("Dirvert", &fsihist->dirvert);
  // OutputTree->SetBranchAddress("Abspvert", &fsihist->abspvert);
  // OutputTree->SetBranchAddress("Abstpvert", &fsihist->abstpvert);
  // OutputTree->SetBranchAddress("Ipvert", &fsihist->ipvert);
  // OutputTree->SetBranchAddress("Iverti", &fsihist->iverti);
  // OutputTree->SetBranchAddress("Ivertf", &fsihist->ivertf);
  // OutputTree->SetBranchAddress("Fsiprob", &fsihist->fsiprob);

  // OutputTree->SetBranchAddress("Numbndn", &neuttarget->numbndn);
  // OutputTree->SetBranchAddress("Numbndp", &neuttarget->numbndp);
  // OutputTree->SetBranchAddress("Numfrep", &neuttarget->numfrep);
  // OutputTree->SetBranchAddress("Numatom", &neuttarget->numatom);
  // OutputTree->SetBranchAddress("Ibound", Ibound);//

  // OutputTree->SetBranchAddress("Nfnvert", &nucleonfsihist->nfnvert);
  // OutputTree->SetBranchAddress("Nfiflag", &nucleonfsihist->nfiflag);
  // OutputTree->SetBranchAddress("Nfx", &nucleonfsihist->nfx);
  // OutputTree->SetBranchAddress("Nfy", &nucleonfsihist->nfy);
  // OutputTree->SetBranchAddress("Nfz", &nucleonfsihist->nfz);
  // OutputTree->SetBranchAddress("Nfpx", &nucleonfsihist->nfpx);
  // OutputTree->SetBranchAddress("Nfpy", &nucleonfsihist->nfpy);
  // OutputTree->SetBranchAddress("Nfpz", &nucleonfsihist->nfpz);
  // OutputTree->SetBranchAddress("Nfe", &nucleonfsihist->nfe);
  // OutputTree->SetBranchAddress("Nffirststep", &nucleonfsihist->nffirststep);
  // OutputTree->SetBranchAddress("Nfnstep", &nucleonfsihist->nfnstep);
  // OutputTree->SetBranchAddress("Nfecms2", &nucleonfsihist->nfecms2);
  // OutputTree->SetBranchAddress("Nfptot", &nucleonfsihist->nfptot);

  OutputTree->Branch("Npvc",&vcwork->nvc,"Npvc/I"); //number of particles
  OutputTree->Branch("Ipvc",vcwork->ipvc,"Ipvc[100]/I"); //identity of particle
  OutputTree->Branch("Ichvc",vcwork->icrnvc,"Ichvc[100]/I");
  OutputTree->Branch("Iorgvc", vcwork->iorgvc,"Iorgvc[100]/I");//identity of parent particle
  OutputTree->Branch("Iflvc", vcwork->iflgvc, "Iflvc[100]/I");
  OutputTree->Branch("Abspvc", Abspvc, "Abspvc[100]/F");//absolute momentum of particle
  OutputTree->Branch("Pvc", vcwork->pvc, "Pvc[100][3]/F");//XYZ momentum of particle
  OutputTree->Branch("Posvc", vcwork->posvc, "Posvc[3]/F");//XYZ position of particle

  OutputTree->Branch("Nvert", &fsihist->nvert, "Nvert/I");//number of vertices
  OutputTree->Branch("Posvert", fsihist->posvert, "Posvert[100][3]/F");//XYZ position of vertex
  OutputTree->Branch("Iflgvert", fsihist->iflgvert, "Iflgvert[100]/I");
  OutputTree->Branch("Nvcvert", fsihist->nvcvert, "Nvcvert/I");
  OutputTree->Branch("Dirvert", fsihist->dirvert, "Dirvert[300][3]/F");
  OutputTree->Branch("Abspvert", fsihist->abspvert, "Abspvert[300]/F");
  OutputTree->Branch("Abstpvert", fsihist->abstpvert, "Abstpvert[300]/F");
  OutputTree->Branch("Ipvert", fsihist->ipvert, "Ipvert[300]/I");
  OutputTree->Branch("Iverti", fsihist->iverti, "Iverti[300]/I");
  OutputTree->Branch("Ivertf", fsihist->ivertf, "Ivertf[300]/I");
  OutputTree->Branch("Fsiprob", &fsihist->fsiprob, "Fsiprob/F");

  OutputTree->Branch("Numbndn", &neuttarget->numbndn, "Numbndn/I");
  OutputTree->Branch("Numbndp", &neuttarget->numbndp, "Numbndp/I");
  OutputTree->Branch("Numfrep", &neuttarget->numfrep, "Numfrep/I");
  OutputTree->Branch("Numatom", &neuttarget->numatom, "Numatom/I");
  OutputTree->Branch("Ibound", &Ibound, "Ibound/I");//

  OutputTree->Branch("Nfnvert", &nucleonfsihist->nfnvert, "Nfnvert/I");
  OutputTree->Branch("Nfiflag", nucleonfsihist->nfiflag, "Nfiflag[200]/I");
  OutputTree->Branch("Nfx", nucleonfsihist->nfx, "Nfx[200]/F");
  OutputTree->Branch("Nfy", nucleonfsihist->nfy, "Nfy[200]/F");
  OutputTree->Branch("Nfz", nucleonfsihist->nfz, "Nfz[200]/F");
  OutputTree->Branch("Nfpx", nucleonfsihist->nfpx, "Nfpx[200]/F");
  OutputTree->Branch("Nfpy", nucleonfsihist->nfpy, "Nfpy[200]/F");
  OutputTree->Branch("Nfpz", nucleonfsihist->nfpz, "Nfpz[200]/F");
  OutputTree->Branch("Nfe", nucleonfsihist->nfe, "Nfe[200]/F");
  OutputTree->Branch("Nffirststep", nucleonfsihist->nffirststep, "Nffirststep[200]/I");
  OutputTree->Branch("Nfnstep", &nucleonfsihist->nfnstep, "Nfnstep/I");
  OutputTree->Branch("Nfecms2", nucleonfsihist->nfecms2, "Nfecms2[2000]/F");
  OutputTree->Branch("Nfptot", nucleonfsihist->nfptot, "Nfptot[2000]/F");

  if((nevccard->mpvevct != 1) && (nevccard->mpvevct != 2)){
    std::cerr << "[ERROR] In card file, MPV: " << nevccard->mpvevct
	      << ", should be 1 or 2." << std::endl;
    return 1;
  }

  RMax = (eftarget->c * 2);

  int number = -1;
  //Generate (nevccard->nectecvt) events.
  std::cout << "[INFO]: Generating:  " << nevccard->nectecvt << " events."
	    << std::endl;
  for(Int_t NEv = 0; NEv < nevccard->nectecvt; ++NEv){
    if ( NEv % 1000 == 0) std::cout << NEv << "/" << nevccard->nectecvt << std::endl;

    // /src/skmcsvc/vcclcm.F
    // Clears VCWORK and VCVRTX common block arrays
    vcclcm_();
    // /src/nuceff/efclfsi.F
    // Clears FSIHIST common block arrays
    efclfsi_();


    //memset stomps output blocks
    Clear_CommonBlock(NEUTEventCommon);
    Clear_CommonBlock(VectorCommon);
    Clear_CommonBlock(FSIHistoryCommon);
    Clear_CommonBlock(NucleusTargetCommon);
    Clear_CommonBlock(NucleonFSIHistoryCommon);

    //nrintr_ sets nucleonfsihist->nfnvert and nucleonfsihist->nfnstep to 0
    // which efectively clears the nucleon fsi history common block.

    SetRandPosition(VectorCommon.PosV);

    PTot = 0;

    //(1) Use mono energetic  at nevccard->pvevct[0], or
    //(2) Throw momentum between nevccard->pvevct[0] and nevccard->pvevct[1].
    if(nevccard->mpvevct == 1){
      PTot = nevccard->pvevct[0];
    } else if (nevccard->mpvevct == 2){
      PTot = (nevccard->pvevct[1] - nevccard->pvevct[0])*rlu_(0)
        + nevccard->pvevct[0];
    }
  
    vcwork->nvc = 4;

    //Setting up the 2nd particle in the VC work blocks
    vcwork->iorgvc[1] = -1;
    vcwork->iflgvc[1] = 0;
    vcwork->icrnvc[1] = 0;
    vcwork->timvc[1] = 0;

    vcwork->pvc[1][0] = PTot;
    vcwork->pvc[1][1] = 0;
    vcwork->pvc[1][2] = 0;


    //Setting up the 4th particle in the VC work block
    //This effectively fires a nucleon particle gun from VectorCommon.PosV
    //as vcwork->icrnvc[3] = 1 tells nrintr_ to 'chase' and simulate nucleon
    //FSI on this particle.
    vcwork->ipvc[3] = IdNucl;
    mcmass_( &(vcwork->ipvc[3]), &(vcwork->amasvc[3]) );

    vcwork->pvc[3][0] = PTot;
    vcwork->pvc[3][1] = 0;
    vcwork->pvc[3][2] = 0;

    posinnuc->ibound = 1;
    posinnuc->posnuc[3][0] = VectorCommon.PosV[0];
    posinnuc->posnuc[3][1] = VectorCommon.PosV[1];
    posinnuc->posnuc[3][2] = VectorCommon.PosV[2];

    vcwork->iorgvc[3] = -1;
    vcwork->iflgvc[3] = 0;
    vcwork->icrnvc[3] = 1; // chase this particle
    vcwork->timvc[3] = 0;
    //std::cout << "prenrintr = "<< vcwork->nvc << "           ";
    //Simulate nucleon FSI interactions
    // /src/nuccorspl/nrintr.F
    nrintr_();
    //std::cout << "postnrintr = "<< vcwork->nvc << std::endl;
    //Fill NE work output block
    NEUTEventCommon.IEvent = NEv;
    number = vcwork->ipvc[5];
    Copy_NEUTEventCommon_F2C();
  
    //Fill VC work output block
    Copy_VectorCommon_F2C();

    //src/nuceff/evpiprob.F
    evpiprob_();
    //Fill FSI History Output block
    Copy_FSIHistoryCommon_F2C();

    //Nuclear Target Output block
    Copy_NucleusTargetCommon_F2C();

    //Nucleon FSI History output block
    Copy_NucleonFSIHistoryCommon_F2C();
  
    for(Int_t n1 = 3; (n1 < vcwork->nvc) && (n1 < 100); ++n1){

      Abspvc[n1]=sqrt(vcwork->pvc[n1][0]*vcwork->pvc[n1][0]
		      + vcwork->pvc[n1][1]*vcwork->pvc[n1][1]
		      + vcwork->pvc[n1][2]*vcwork->pvc[n1][2]);
    }
    OutputTree->Fill();

    bool cat = true;
    //    if(vcwork->nvc==8)
      {
	//	std::cout << "Npvc == " << vcwork->nvc << std::endl;
	cat = true;
	for(int d = 0; d<vcwork->nvc; d++)
	  {
	    if(vcwork->ipvc[d]==211||vcwork->ipvc[d]==111||vcwork->ipvc[d]==-211)
	      {
		cat =false;
	      }
	  }
	    if(cat == true)
	      {
	for(int g = 0; g<vcwork->nvc; g++)
	  //for(int z = 3; z<4; z++)
	  { 
	    //if(vcwork->icrnvc[z] != 0 || vcwork->ipvc[z] < 2112|| vcwork->ipvc[z] > 2212 || vcwork->iorgvc[z] != -1 || vcwork->iflgvc[z] != 7)
	      
		//if(vcwork->ipvc[z]!=2212&&vcwork->ipvc[z]!=2112 || vcwork->icrnvc[z]!=1||vcwork->iorgvc[z]!=4||vcwork->iflgvc[z]!=0)
		// if(vcwork->icrnvc[z]!=0||vcwork->iorgvc[z]!=-1||vcwork->iflgvc[z]!=7)
		//if(vcwork->ipvc[4]==111||vcwork->ipvc[4]==211||vcwork->ipvc[4]==-111)
		//		if(vcwork->iflgvc[z] != -1 && vcwork->iflgvc[z] != 0  &&vcwork->iflgvc[z] != 3 && vcwork->iflgvc[z] != 4 && vcwork->iflgvc[z] != 7 &&vcwork->iflgvc[z] != 8 &&vcwork->iflgvc[z] != 9)
		//if(vcwork->ipvc[z]!=211 && vcwork->ipvc[z]!=111 && vcwork->ipvc[z]!=-211 && vcwork->ipvc[z]!=-111)
	    //	    if(vcwork->iflgvc[g]!=0&&vcwork->iflgvc[g]!= 7)
	    if(vcwork->iflgvc[g]==9)
	    {
		      //if(vcwork->nvc==)	

		    
		 
		
		  
		    for(int z = 0; z<vcwork->nvc; z++)
		      {
			
			std::cout << "nvc = " << vcwork->nvc << std::endl;
			std::cout <<"zth particle = " << z <<" Ipvc == " << vcwork->ipvc[z] << " Ichvc == " << vcwork->icrnvc[z] << " Iorgvc == "<< vcwork->iorgvc[z] << " Iflvc == " << vcwork->iflgvc[z]<< " AbsMom == " << Abspvc[z] << std::endl;
		      }
		    
		    
		    
		//std::cout << std::endl << std::endl;
		    
		    
		    std::cout << std::endl << std::endl;
		  }		   
	    
	  }
      }
      }    
      
  }
  
  

  OutputTree->Write();
  OutputFile->Close();

  return 0;
}
