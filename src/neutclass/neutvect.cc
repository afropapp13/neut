#include "neutpart.h"
#include "neutfsipart.h"
#include "neutfsivert.h"
#include "neutnucfsivert.h"
#include "neutnucfsistep.h"
#include "neutvect.h"

#include <iostream>

NeutVect::NeutVect(Int_t np, Int_t nfsip, Int_t nfsiv)
{

  ClearVars();

  fPartInfo = NULL;
  fFsiPartInfo = NULL;
  fFsiVertInfo = NULL;
  fNucFsiVertInfo = NULL;
  fNucFsiStepInfo = NULL;

  fZeroVect.SetXYZM(0.,0.,0.,0.);

  SetNpart(np);
  SetNfsiPart(nfsip);
  SetNfsiVert(nfsiv);

}

void
NeutVect::ClearVars()
{
  // TO DO: Implement for all variables

  EventNo  =  0;

  TargetA  =  0;
  TargetZ  =  0;
  TargetH  =  0;
  Ibound   =  -1;

  VNuclIni =  0;
  VNuclFin =  0;

  PFSurf   =  0;
  PFMax    =  0;
  
  FluxID   = -1;

  ////////////////////////////////////////////
  Mode     =  0;

  Totcrs   =  0;

  Fsiprob  =  1;

  fVertexID.Set(0);
  fParentIdx.Set(0);

  {
	Int_t i;
	for(i = 0 ; i < 16 ; i++){
	  fRandomSeed[i] = 0;
	}
  }

}

void
NeutVect::SetNprimary(Int_t nprimary)
{
  fNprimary = nprimary;
}

void
NeutVect::SetNpart(Int_t npart)
{
  fNpart = npart;

  if (fNpart > 0){
	if (fPartInfo == NULL){
	  fPartInfo = new TObjArray(fNpart);
	  fPartInfo->SetOwner(kTRUE);
	}else{
	  fPartInfo->Expand(fNpart);
	}
	
	fVertexID.Set(fNpart);
	
	fParentIdx.Set(fNpart);
	
	/*
	if (fPosIni == NULL){
	  fPosIni = new TLorentzVector[fNpart];
	}
	if (fPosFin == NULL){
	  fPosFin = new TLorentzVector[fNpart];
	}
	*/
  }else{
	delete fPartInfo;
	fPartInfo = NULL;
	ClearVars();
  }
}


void
NeutVect::SetNfsiPart(Int_t npart)
{
  fNfsiPart = npart;

  if (fNfsiPart > 0){
	if (fFsiPartInfo == NULL){
	  fFsiPartInfo = new TObjArray(fNfsiPart);
	  fFsiPartInfo->SetOwner(kTRUE);
	}else{
	  fFsiPartInfo->Expand(fNfsiPart);
	}

	
  }else{
	delete fFsiPartInfo;
	fFsiPartInfo = NULL;
	//ClearVars();
  }
}

void
NeutVect::SetNfsiVert(Int_t nvert)
{
  fNfsiVert = nvert;

  if (fNfsiVert > 0){
	if (fFsiVertInfo == NULL){
	  fFsiVertInfo = new TObjArray(fNfsiVert);
	  fFsiVertInfo->SetOwner(kTRUE);
	}else{
	  fFsiVertInfo->Expand(fNfsiVert);
	}

	
  }else{
	delete fFsiVertInfo;
	fFsiVertInfo = NULL;
	//ClearVars();
  }
}

void
NeutVect::SetVertexID(Int_t idx , Int_t vtxid)
{
  fVertexID.AddAt(vtxid,idx);
};

void
NeutVect::SetVertexID(Int_t npart , Int_t *vtxid_array)
{
  fVertexID.Set(npart, vtxid_array);
};

void
NeutVect::SetParentIdx(Int_t np , Int_t idx)
{
  fParentIdx.AddAt(idx,np);
};

void
NeutVect::SetParentIdx(Int_t npart , Int_t *idx_array)
{
  fParentIdx.Set(npart, idx_array);
};

void
NeutVect::SetPartInfo(int idx, NeutPart PInfo)
{

  NeutPart *PInfo_p;

  if (fNpart<=idx){
  	SetNpart(idx);
  }

  PInfo_p = new NeutPart();

  *PInfo_p = PInfo;
  
  fPartInfo->AddAt(PInfo_p,idx);

}

void
NeutVect::SetPartInfo(Int_t npart, NeutPart *PInfo_array)
{

  Int_t i;

  for (i = 0 ; i < npart ; i++){
	SetPartInfo(i,PInfo_array[i]);
  }

};

void
NeutVect::SetFsiPartInfo(int idx, NeutFsiPart PInfo)
{

  NeutFsiPart *PInfo_p;

  if (fNfsiPart<=idx){
  	SetNfsiPart(idx);
  }

  PInfo_p = new NeutFsiPart();

  *PInfo_p = PInfo;
  
  fFsiPartInfo->AddAt(PInfo_p,idx);

}

void
NeutVect::SetFsiPartInfo(Int_t npart, NeutFsiPart *PInfo_array)
{

  Int_t i;

  for (i = 0 ; i < npart ; i++){
	SetFsiPartInfo(i,PInfo_array[i]);
  }

};


void
NeutVect::SetFsiVertInfo(int idx, NeutFsiVert VInfo)
{

  NeutFsiVert *VInfo_p;

  if (fNfsiVert<=idx){
  	SetNfsiVert(idx);
  }

  VInfo_p = new NeutFsiVert();

  *VInfo_p = VInfo;
  
  fFsiVertInfo->AddAt(VInfo_p,idx);

}

void
NeutVect::SetFsiVertInfo(Int_t nvert, NeutFsiVert *VInfo_array)
{

  Int_t i;

  for (i = 0 ; i < nvert ; i++){
	SetFsiVertInfo(i,VInfo_array[i]);
  }

};

void
NeutVect::SetNnucFsiVert(Int_t nnucvert)
{
  fNnucFsiVert = nnucvert;

  if (fNnucFsiVert > 0){
	if (fNucFsiVertInfo == NULL){
	  fNucFsiVertInfo = new TObjArray(fNnucFsiVert);
	  fNucFsiVertInfo->SetOwner(kTRUE);
	}else{
	  fNucFsiVertInfo->Expand(fNnucFsiVert);
	}

	
  }else{
	delete fNucFsiVertInfo;
	fNucFsiVertInfo = NULL;
	//ClearVars();
  }
}

void
NeutVect::SetNucFsiVertInfo(int idx, NeutNucFsiVert VInfo)
{

  NeutNucFsiVert *VInfo_p;

  if (fNnucFsiVert<=idx){
  	SetNnucFsiVert(idx);
  }

  VInfo_p = new NeutNucFsiVert();

  *VInfo_p = VInfo;
  
  fNucFsiVertInfo->AddAt(VInfo_p,idx);

}

void
NeutVect::SetNucFsiVertInfo(Int_t nvert, NeutNucFsiVert *VInfo_array)
{

  Int_t i;

  for (i = 0 ; i < nvert ; i++){
	SetNucFsiVertInfo(i,VInfo_array[i]);
  }

};

void
NeutVect::SetNnucFsiStep(Int_t nnucstep)
{
  fNnucFsiStep = nnucstep;

  if (fNnucFsiStep > 0){
	if (fNucFsiStepInfo == NULL){
	  fNucFsiStepInfo = new TObjArray(fNnucFsiStep);
	  fNucFsiStepInfo->SetOwner(kTRUE);
	}else{
	  fNucFsiStepInfo->Expand(fNnucFsiStep);
	}

	
  }else{
	delete fNucFsiStepInfo;
	fNucFsiStepInfo = NULL;
	//ClearVars();
  }
}

void
NeutVect::SetNucFsiStepInfo(int idx, NeutNucFsiStep VInfo)
{

  NeutNucFsiStep *VInfo_p;

  if (fNnucFsiStep<=idx){
  	SetNnucFsiStep(idx);
  }

  VInfo_p = new NeutNucFsiStep();

  *VInfo_p = VInfo;
  
  fNucFsiStepInfo->AddAt(VInfo_p,idx);

}

void
NeutVect::SetNucFsiStepInfo(Int_t nvert, NeutNucFsiStep *VInfo_array)
{

  Int_t i;

  for (i = 0 ; i < nvert ; i++){
	SetNucFsiStepInfo(i,VInfo_array[i]);
  }

};

NeutVect::~NeutVect()
{

  if (fPartInfo) delete fPartInfo;
  if (fFsiVertInfo) delete fFsiVertInfo;
  if (fFsiPartInfo) delete fFsiPartInfo;

}

void
NeutVect::Dump()
{

  int i, np;

  std::cout << "Event #        :" << EventNo  << "\n";
  
  std::cout << "Target A       :" << TargetA  << "\n";
  std::cout << "Target Z       :" << TargetZ  << "\n";
  std::cout << "Target Free H  :" << TargetH  << "\n";
  std::cout << "Bound nucleon  :" << Ibound  << "\n";
  std::cout << "VNuclIni       :" << VNuclIni << "\n";
  std::cout << "VNuclFin       :" << VNuclFin << "\n";
  std::cout << "PF Surface     :" << PFSurf   << "\n";
  std::cout << "PF Maximum     :" << PFMax    << "\n";
  std::cout << "Flux ID        :" << FluxID   << "\n";

  std::cout << "QE Model       :" << QEModel   << "\n";
  std::cout << "SPI Model      :" << SPIModel  << "\n";
  std::cout << "COH Model      :" << COHModel  << "\n";
  std::cout << "DIS Model      :" << DISModel  << "\n";

  std::cout << "QE Vector F/F  :" << QEVForm   << "\n";
  std::cout << "QE MA          :" << QEMA      << "\n";
  std::cout << "SPI MA         :" << SPIMA     << "\n";
  std::cout << "QE MV          :" << QEMV      << "\n";
  std::cout << "SPI MV         :" << SPIMV     << "\n";

  std::cout << "KAPPA          :" << KAPPA     << "\n";
  std::cout << "COH MA         :" << COHMA     << "\n";
  std::cout << "COH R0         :" << COHR0     << "\n";
  std::cout << "COH A1 err     :" << COHA1err  << "\n";
  std::cout << "COH b1 err     :" << COHb1err  << "\n";

  std::cout << "neutcore ver.  :" << COREVer   << "\n";
  std::cout << "nuceff ver.    :" << NUCEVer   << "\n";
  std::cout << "neccorspl ver. :" << NUCCVer   << "\n";

  std::cout << "Fermi motion   :" << FrmFlg   << "\n";
  std::cout << "Pauli blocking :" << PauFlg   << "\n";
  std::cout << "Nuceff in O16  :" << NefO16   << "\n";
  std::cout << "Flag on mode   :" << ModFlg   << "\n";
  std::cout << "Formation zone :" << FormLen  << "\n";
  std::cout << "Pi-less Decay  :" << IPilessDcy  << "\n";
  std::cout << "Pi-less Decay R:" << RPilessDcy  << "\n";
  std::cout << "nucleon FSI    :" << NucScat  << "\n";
  std::cout << "FSI Xsec factor:" << NucFac   << "\n";

  std::cout << "FSI LowE QE MFP Scaling          :" << NuceffFactorPIQE     << "\n";
  std::cout << "FSI Hadron Production MFP Scaling:" << NuceffFactorPIInel	  << "\n";
  std::cout << "FSI Absorption MFP Scaling       :" << NuceffFactorPIAbs 	  << "\n";
  std::cout << "FSI HighE QE MFP Scaling         :" << NuceffFactorPIQEH 	  << "\n";
  std::cout << "FSI LowE CX MFP Scaling          :" << NuceffFactorPICX 	  << "\n";
  std::cout << "FSI HighE CX MFP Scaling         :" << NuceffFactorPICXH    << "\n";
  std::cout << "FSI Forward Scatter MFP Scaling  :" << NuceffFactorPICoh	  << "\n";
  std::cout << "FSI HighE QE Kinematic Factor    :" << NuceffFactorPIQEHKin << "\n";
  std::cout << "FSI LowE QE Kinematic Factor     :" << NuceffFactorPIQELKin << "\n";
  std::cout << "FSI HighE CX Mix Factor          :" << NuceffFactorPICXKin  << "\n";
  std::cout << "FSI Total Mean Free Path Factor  :" << NuceffFactorPIAll  << "\n";

  
  std::cout << "Intr. mode     :" << Mode   << "\n";
  std::cout << "Total Xsec     :" << Totcrs   << " x 10^(-38) cm^2\n";

  std::cout << "# of particles =" << Npart() << "\n";	
  std::cout << "# of primary particles =" << Nprimary() << "\n";	

  np = Npart();

  for ( i = 0 ; i < np ; i++ ){
	std::cout << "i=" << i << "\n";
	std::cout << "Vertex         =" << VertexID(i) << "\n";		
	std::cout << "Parent Index   =" << ParentIdx(i) << "\n";

	std::cout << "Particle Code  = " << (PartInfo(i))->fPID   << "\n";
	std::cout << "Particle Mass  = " << (PartInfo(i))->fMass   << "\n";
	std::cout << "Particle Mom.  =(" << (PartInfo(i))->fP.Px() << ","
			  << (PartInfo(i))->fP.Py() << "," 
			  << (PartInfo(i))->fP.Pz() << ","
			  << (PartInfo(i))->fP.E()  << ")"
			  << "\n";
	std::cout << "Particle Flag  = " << (PartInfo(i))->fIsAlive << "\n";
	std::cout << "Particle Stat. = " << (PartInfo(i))->fStatus  << "\n";
	std::cout << "Particle Pos(1)=(" << (PartInfo(i))->fPosIni.X() << ","
			  << (PartInfo(i))->fPosIni.Y() << "," 
			  << (PartInfo(i))->fPosIni.Z() << ","
			  << (PartInfo(i))->fPosIni.T()  << ")"
			  << "\n";
	std::cout << "Particle Pos(2)=(" << (PartInfo(i))->fPosFin.X() << ","
			  << (PartInfo(i))->fPosFin.Y() << "," 
			  << (PartInfo(i))->fPosFin.Z() << ","
			  << (PartInfo(i))->fPosFin.T()  << ")"
			  << "\n";
  }

  np = NfsiPart();

  std::cout << "\n # of FSI particles =" << np << "\n";

  for ( i = 0 ; i < np ; i++ ){
	std::cout << "i=" << i << "\n";

	std::cout << "Particle Code  = " << (FsiPartInfo(i))->fPID   << "\n";
	std::cout << "Particle Dir.  =(" << (FsiPartInfo(i))->fDir.Px() << ","
			  << (FsiPartInfo(i))->fDir.Py() << "," 
			  << (FsiPartInfo(i))->fDir.Pz() << ")"
			  << "\n";
  	std::cout << "Momentum (Lab)        = " << (FsiPartInfo(i))->fMomLab   << "\n"; 
	std::cout << "Momentum (Nuc. Rest)  = " << (FsiPartInfo(i))->fMomNuc   << "\n";
	std::cout << "Starting->Ending Vertex  = " << (FsiPartInfo(i))->fVertStart << "->" << (FsiPartInfo(i))->fVertEnd  << "\n";

  }



  np = NfsiVert();

  std::cout << "\n # of FSI vertices =" << np << "\n";

  for ( i = 0 ; i < np ; i++ ){
    std::cout << "Vertex # " << i+1 << ": ";
    std::cout << "Type=" << (FsiVertInfo(i))->fVertID << ", Pos=(";
	std::cout << (FsiVertInfo(i))->fPos.Px() << ","
		  << (FsiVertInfo(i))->fPos.Py() << "," 
		  << (FsiVertInfo(i))->fPos.Pz() << ")"
		  << "\n";
  }

  std::cout << "FSI probability :" << Fsiprob   << "\n";
  
}


ClassImp(NeutVect)
