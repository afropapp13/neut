#include <iostream>
#include <NeutRootHandlers.h>

#include "neutvect.h"
#include "neutvtx.h"

#include "neworkC.h"
#include "vcworkC.h"
#include "vcvrtxC.h"
#include "neutparamsC.h"
#include "posinnucC.h"
#include "neutmodelC.h"
#include "nieves1p1hC.h"
#include "necardC.h"
#include "nrcardC.h"
#include "nefillverC.h"
#include "neutcrsC.h"
#include "posinnucC.h"
#include "fsihistC.h"
#include "nucleonfsihistC.h"
#include "fsinucleonreweightC.h"
#include "neutcrsC.h"
#include "nrintC.h"

extern "C"
{
  int neutfillvect_(char * /* filename */, char * /* treename */,
					char * /* Branch name */, int , int, int);

  int neutfillvect(char * /* filename */, char * /* treename */,
				   char * /* Branch name */);

  int neutfillvtx_(char * /* filename */, char * /* treename */,
					char * /* Branch name */, int , int, int);

  int neutfillvtx(char * /* filename */, char * /* treename */,
				   char * /* Branch name */);


  int neutfill_(char * /* filename */, char * /* treename */,
				char * /* Vector Branch name */,
				char * /* Vertex Branch name */,
				int, int, int, int);

  int neutfill(char * /* filename */, char * /* treename */,
			   char * /* Vector Branch name */,
			   char * /* Vertex Branch name */);

  void nefillver_(void);
};


int
neutfill_(char *filename,      char *treename, 
		  char *vecbranchname, char *vtxbranchname,
		  int len_fname,int len_tname, 
		  int len_vecbname, int len_vtxbname)
{
  NeutRootHandlers nrh;  

  nrh.nulltermstr(filename,      len_fname);
  nrh.nulltermstr(treename,      len_tname);
  nrh.nulltermstr(vecbranchname, len_vecbname);
  nrh.nulltermstr(vtxbranchname, len_vtxbname);

  return neutfill(filename, treename, vecbranchname, vtxbranchname);
}

int
neutfill(char *filename,      char *treename, 
		  char *vecbranchname, char *vtxbranchname)
{
  NeutRootHandlers nrh;  
  static TTree   *t   = NULL;

  int ret;

  ret = neutfillvect(filename, treename, vecbranchname);
  if (ret != 0){
	return -1;
  }
  ret = neutfillvtx(filename, treename,  vtxbranchname);  
  if (ret != 0){
	return -2;
  }
  
  if (t == NULL){
	t = nrh.attachtree(filename, treename);
	if (t == NULL){
	  return -3;
	}
  }

  t->Fill();
  return 0;

}


int
neutfillvect_(char *filename, char *treename, char *branchname,
			  int len_fname,int len_tname, int len_bname)
{
  NeutRootHandlers nrh;  

  nrh.nulltermstr(filename,   len_fname);
  nrh.nulltermstr(treename,   len_tname);
  nrh.nulltermstr(branchname, len_bname);

  return neutfillvect(filename, treename, branchname);
}

int
neutfillvect(char *filename, char *treename, char *branchname)
{

  NeutRootHandlers nrh;

  static TTree   *t   = NULL;
  static NeutVect *nv = NULL;

  static Int_t event_no= 1;

  NeutPart pinfo;
  
  NeutFsiVert fsivinfo;
  NeutFsiPart fsipinfo;

  NeutNucFsiVert nucfsivinfo;
  NeutNucFsiStep nucfsisinfo;  

  int i,j;

  nefillver_();

  if (t == NULL){
	t = nrh.maketree(filename, treename, "Neut Tree");
	if (t == NULL) return -1;

	nv = new NeutVect();
	t->Branch(branchname,"NeutVect",&nv,1024,99);
  }
  
  nv->SetNpart(0);
  nv->SetNfsiPart(0);
  nv->SetNfsiVert(0);
  nv->SetNnucFsiVert(0);
  nv->SetNnucFsiStep(0);
  
  /****************************************************/
  nv->EventNo = event_no++;
  
  nv->TargetA  = neuttarget_.numatom;
  nv->TargetZ  = neuttarget_.numbndp;
  nv->TargetH  = neuttarget_.numfrep;
  nv->VNuclIni = nenupr_.vnuini;
  nv->VNuclFin = nenupr_.vnufin;
  nv->PFSurf   = nenupr_.pfsurf;
  nv->PFMax    = nenupr_.pfmax;
  nv->Ibound   = posinnuc_.ibound;
  
  nv->FluxID   = 0;
  
  /****************************************************/
  // interaction model and some parameters
  nv->QEModel  = nemdls_.mdlqe;
  nv->SPIModel = nemdls_.mdlspi;

  nv->COHModel = nemdls_.mdlcoh;
  nv->DISModel = nemdls_.mdldis;

  nv->QEVForm  = nemdls_.mdlqeaf;
  nv->QEMA     = nemdls_.xmaqe;

  nv->SPIForm  = neut1pi_.neiff;
  nv->RESForm  = 0;
  
  nv->SPINRType= neut1pi_.nenrtype;
  nv->RESNRType= 0;

  nv->SPICA5I  = neut1pi_.rneca5i;
  nv->SPIBGScale= neut1pi_.rnebgscl;

  nv->SPIMA    = nemdls_.xmaspi;
  nv->SPIMV    = nemdls_.xmvspi;

  nv->RESMA   = nemdls_.xmares;
  nv->RESMV   = nemdls_.xmvres;

  nv->QEMV     = nemdls_.xmvqe;

  nv->KAPPA    = nemdls_.kapp;

  nv->COHMA    = nemdls_.xmacoh;
  nv->COHR0    = nemdls_.rad0nu;
  nv->COHA1err = nemdls_.fa1coh;
  nv->COHb1err = nemdls_.fb1coh;

  // neut verion
  nv->COREVer  = neutversion_.corev;
  nv->NUCEVer  = neutversion_.nucev;
  nv->NUCCVer  = neutversion_.nuccv;
  // neut card
  nv->FrmFlg  = neutcard_.nefrmflg;
  nv->PauFlg  = neutcard_.nepauflg;
  nv->NefO16  = neutcard_.nenefo16;
  nv->ModFlg  = neutcard_.nemodflg;
  nv->SelMod  = neutcard_.neselmod;
  nv->FormLen = nenupr_.iformlen;
  nv->IPilessDcy = neutpiless_.ipilessdcy;
  nv->RPilessDcy = neutpiless_.rpilessdcy;
  nv->NucScat = nucres_.nucrescat;
  nv->NucFac  = nucres_.xnucfact;

  nv->NuceffKinVersion = nuceffver_.nefkinver;

  nv->NuceffFactorPIQE = neffpr_.fefqe;	  
  nv->NuceffFactorPIInel = neffpr_.fefinel;  
  nv->NuceffFactorPIAbs = neffpr_.fefabs;  
  nv->NuceffFactorPIQEH = neffpr_.fefqeh;
  nv->NuceffFactorPICX = neffpr_.fefcx;  
  nv->NuceffFactorPICXH = neffpr_.fefcxh;  
  nv->NuceffFactorPICoh = neffpr_.fefcoh;
  nv->NuceffFactorPIQEHKin = neffpr_.fefqehf; 	  
  nv->NuceffFactorPIQELKin = neffpr_.fefcohf; 
  nv->NuceffFactorPICXKin = neffpr_.fefcxhf; 
  nv->NuceffFactorPIAll = neffpr_.fefall;
  
  nv->NrintNucleonCascadeProb = nrint_.pcascprob;

  /****************************************************/

  nv->NVQERFG          = nievesqepar_.nvqerfg;
  nv->NVQEBind     	   = nievesqepar_.nvqebind;
  nv->NVQERPA		   = nievesqepar_.nvqerpa;
  nv->XNVRPAFP0in      = nievesqepar_.xnvrpafp0in;
  nv->XNVRPAPF0ex      = nievesqepar_.xnvrpapf0ex;
  nv->XNVRPAFstar      = nievesqepar_.xnvrpafstar;
  nv->XNVRPAF	       = nievesqepar_.xnvrpaf;
  nv->XNVRPAPILambda   = nievesqepar_.xnvrpapilambda;
  nv->XNVRPACR0        = nievesqepar_.xnvrpacr0;
  nv->XNVRPARHOLambda  = nievesqepar_.xnvrparholambda;
  nv->XNVRPAGp         = nievesqepar_.xnvrpagp;
  nv->XNVRPAXMPI       = nievesqepar_.xnvrpaxmpi;
  nv->XNVRPAXMRHO      = nievesqepar_.xnvrpaxmrho;
  nv->XNVRPAIrel       = nievesqepar_.xnvrpairel;
  nv->FFTYPE           = nievesqepar_.fftype;
  nv->NVBINDFermiCor   = nievesqepar_.nvbindfermicor;

  /****************************************************/
  nv->Mode = nework_.modene;
  nv->Totcrs = neutcrscom_.totcrsne;

  nv->CrsEnergy = neutcrscom_.crsenergy;
  
  for ( i = 0 ; i < 8 ; i++ ){
	nv->DifCrsNE[i] = neutcrscom_.difcrsne[i];
  }
  nv->Crsx   = neutcrscom_.crsx;
  nv->Crsy   = neutcrscom_.crsy;
  nv->Crsz   = neutcrscom_.crsz;
  nv->Crsphi = neutcrscom_.crsphi;
  nv->Crsq2  = neutcrscom_.crsq2;

  Int_t np = vcwork_.nvc;

  nv->SetNpart(np);

  nv->SetNprimary(nework_.numne);

  for ( i = 0 ; i < np ; i++ ){
	Double_t mass;

	nv->SetVertexID(i, vcwork_.ivtivc[i]);

	nv->SetParentIdx(i, vcwork_.iorgvc[i]);
	  
	pinfo.fPID = vcwork_.ipvc[i];
	mass       = vcwork_.amasvc[i];
	pinfo.fMass = mass;
	pinfo.fIsAlive = vcwork_.icrnvc[i];
	pinfo.fStatus  = vcwork_.iflgvc[i];
	
	pinfo.fP.SetXYZM(vcwork_.pvc[i][0],
					 vcwork_.pvc[i][1],
					 vcwork_.pvc[i][2],
					 mass);

	pinfo.fPosIni.SetXYZT(posinnuc_.posnuc[i][0],
						  posinnuc_.posnuc[i][1],
						  posinnuc_.posnuc[i][2],
						  0.);
	
	pinfo.fPosFin.SetXYZT(posinnuc_.posnuc[i][0],
						  posinnuc_.posnuc[i][1],
						  posinnuc_.posnuc[i][2],						  
						  0.);
	  
	nv->SetPartInfo(i, pinfo);
	  
  }
  

  // FSI Particles
  Int_t nfsip = fsihist_.nvcvert;
  
  if (nfsip) {
    nv->SetNfsiPart(nfsip);

    for ( i = 0 ; i < nfsip ; i++ ){
    
      fsipinfo.fPID = fsihist_.ipvert[i];
      fsipinfo.fDir.SetXYZT(fsihist_.dirvert[i][0],
			    fsihist_.dirvert[i][1],
			    fsihist_.dirvert[i][2],
			    0.);
      fsipinfo.fMomLab = fsihist_.abspvert[i];
      fsipinfo.fMomNuc = fsihist_.abstpvert[i];
      fsipinfo.fVertStart = fsihist_.iverti[i];
      fsipinfo.fVertEnd = fsihist_.ivertf[i];
  
      nv->SetFsiPartInfo(i, fsipinfo);
    }
  }
  // Generate dummy object for no FSI block (no pion) events
  // to protect from memory leak when reading from the tree
  // (You'll need to take care of this in downstream routines)
  else {
    nv->SetNfsiPart(1);
    
    fsipinfo.fPID = -1;
    fsipinfo.fDir.SetXYZT(0, 0, 0, 0);
    fsipinfo.fMomLab = -1;
    fsipinfo.fMomNuc = -1;
    fsipinfo.fVertStart = -1;
    fsipinfo.fVertEnd = -1;
    
    nv->SetFsiPartInfo(0, fsipinfo);
  }


  // FSI Vertices
  Int_t nfsiv = fsihist_.nvert;
  
  if (nfsiv) {
    nv->SetNfsiVert(nfsiv);

    for ( i = 0 ; i < nfsiv ; i++ ){
      
      fsivinfo.fPos.SetXYZT(fsihist_.posvert[i][0],
			    fsihist_.posvert[i][1],
			    fsihist_.posvert[i][2],
			    0.);

      fsivinfo.fVertID = fsihist_.iflgvert[i];
  
      nv->SetFsiVertInfo(i, fsivinfo);

    }

    nv->Fsiprob = fsihist_.fsiprob;
  } 
  // Generate dummy object for no FSI block (no pion) events
  // to protect from memory leak when reading from the tree
  // (You'll need to take care of this in downstream routines)
  else {
    nv->SetNfsiVert(1);

    fsivinfo.fPos.SetXYZT(0,0,0,0);

    fsivinfo.fVertID = -99999;
  
    nv->SetFsiVertInfo(0, fsivinfo);

    nv->Fsiprob = -1;
  }

  // Nucleon FSI Verticies
  Int_t nnucfsivert = nucleonfsihist_.nfnvert;
  
  if (nnucfsivert) {
    nv->SetNnucFsiVert(nnucfsivert);

    for ( i = 0 ; i < nnucfsivert ; i++ ){
    
      nucfsivinfo.fVertFlag      = nucleonfsihist_.nfiflag[i];
      nucfsivinfo.fVertFirstStep = nucleonfsihist_.nffirststep[i];

      nucfsivinfo.fPos.SetXYZT(nucleonfsihist_.nfx[i],
							   nucleonfsihist_.nfy[i],
							   nucleonfsihist_.nfz[i],
							   0.);
      nucfsivinfo.fMom.SetPxPyPzE(nucleonfsihist_.nfpx[i],
								  nucleonfsihist_.nfpy[i],
								  nucleonfsihist_.nfpz[i],
								  nucleonfsihist_.nfe[i]);

      nv->SetNucFsiVertInfo(i, nucfsivinfo);
    }
  }
  // Generate dummy object for no FSI block (no pion) events
  // to protect from memory leak when reading from the tree
  // (You'll need to take care of this in downstream routines)
  else {
    nv->SetNnucFsiVert(1);
	nucfsivinfo.fVertFlag      = -1;
	nucfsivinfo.fVertFirstStep = -1;

	nucfsivinfo.fPos.SetXYZT(0.,0.,0.,0.);
	nucfsivinfo.fMom.SetPxPyPzE(0.,0.,0.,0.);
								
    nv->SetNucFsiVertInfo(0, nucfsivinfo);
  }

  // Nucleon FSI Step
  Int_t nnucfsis = nucleonfsihist_.nfnstep;
  
  if ( nnucfsis ) {
    nv->SetNnucFsiStep(nnucfsis);

    for ( i = 0 ; i < nnucfsis ; i++ ){
      
      //          nucfsisinfo.fECMS2 = nucleonfsihist_.nfecms2[i];
      nucfsisinfo.fProb  = nucleonfsihist_.nfptot[i];
      nucfsisinfo.fVertFlagStep = nucleonfsihist_.nfiflagstep[i];
      nucfsisinfo.fVertFsiRhon = nucleonfsihist_.nfirhon[i];
      nucfsisinfo.fStepPel = nucleonfsihist_.nfipel[i];
      nucfsisinfo.fStepPsp = nucleonfsihist_.nfipsp[i];
      nucfsisinfo.fStepPdp = nucleonfsihist_.nfipdp[i];

      //std::cout << "fvertflag step in neutfill.cc = " << nucfsisinfo.fVertFlagStep << std::endl;
      //std::cout << "Step number in neutfill.cc = " << i << std::endl;
      /*      nucfsisinfo.fPosStep.SetXYZT(nucleonfsihist_.nfxstep[i],
			       nucleonfsihist_.nfystep[i],
			       nucleonfsihist_.nfzstep[i],
			       0.);
      std::cout << "fvert fposstep = " << nucfsisinfo.fPosStep.X() << std::endl;
      nucfsisinfo.fMomStep.SetPxPyPzE(nucleonfsihist_.nfpxstep[i],
				  nucleonfsihist_.nfpystep[i],
				  nucleonfsihist_.nfpzstep[i],

      */
      nv->SetNucFsiStepInfo(i, nucfsisinfo);

    }

  } 
  // Generate dummy object for no FSI block (no pion) events
  // to protect from memory leak when reading from the tree
  // (You'll need to take care of this in downstream routines)
  else {
    nv->SetNnucFsiStep(1);

    //    	nucfsisinfo.fECMS2 =  0.;
	nucfsisinfo.fProb  = -1.;
  
    nv->SetNucFsiStepInfo(0, nucfsisinfo);

  }

  if (!(neutcard_.quiet)){
	nv->Dump();
  }

  return 0;
}

int
neutfillvtx_(char *filename, char *treename, char *branchname,
			  int len_fname,int len_tname, int len_bname)
{
  NeutRootHandlers nrh;  

  nrh.nulltermstr(filename,   len_fname);
  nrh.nulltermstr(treename,   len_tname);
  nrh.nulltermstr(branchname, len_bname);

  return neutfillvtx(filename, treename, branchname);
}

int
neutfillvtx(char *filename, char *treename, char *branchname)
{

  NeutRootHandlers nrh;

  static TTree   *t  = NULL;
  static NeutVtx *nv = NULL;
  TLorentzVector tlv;

  static Int_t event_no= 1;
  int nvtx;

  int i;

  if (t == NULL){
	t = nrh.maketree(filename, treename, "Neut Tree");
	if (t == NULL) return -1;

	nv = new NeutVtx();
	t->Branch(branchname,"NeutVtx",&nv,1024,99);
  }

  nv->SetNvtx(0);
  
  /****************************************************/
  nv->EventNo = event_no++;
  
  nvtx  = vcvrtx_.nvtxvc;
  nv->SetNvtx(nvtx);

  for ( i = 0 ; i < nvtx ; i++ ){
	tlv.SetXYZT(vcvrtx_.pvtxvc[i][0],
				vcvrtx_.pvtxvc[i][1],
				vcvrtx_.pvtxvc[i][2],
				vcvrtx_.timvvc[i]);

	nv->SetPos(i, tlv);
	  
  }
  if (!(neutcard_.quiet)){
	nv->Dump();
  }

  return 0;
}
