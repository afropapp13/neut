void
chkreadneutvect()
{

  Int_t i, j;
  Int_t nevents;

  gSystem->Load("neutfsipart.so");
  gSystem->Load("neutfsivert.so");
  gSystem->Load("neutnucfsipart.so");
  gSystem->Load("neutnucfsivert.so");
  gSystem->Load("neutvtx.so");
  gSystem->Load("neutpart.so");
  gSystem->Load("neutvect.so");
  
  TTree *tnv;
  NeutVect *nv;

  TFile f("neutvect.root");
  tnv = (TTree *)(f.Get("neuttree"));
  nv = new NeutVect();

  tnv->SetBranchAddress("vectorbranch",&nv);

  nevents = tnv->GetEntries();

  for ( j = 0 ; j < nevents ; j++ ){

	tnv->GetEntry(j);
	cout << "Event #        :" << nv->EventNo << "\n";
	cout << "Target A       :" << nv->TargetA << "\n";
	cout << "Target Z       :" << nv->TargetA << "\n";
	cout << "VNuclIni       :" << nv->VNuclIni << "\n";
	cout << "VNuclFin       :" << nv->VNuclFin << "\n";
	cout << "PF Surface     :" << nv->PFSurf   << "\n";
	cout << "PF Maximum     :" << nv->PFMax    << "\n";
	cout << "Flux ID        :" << nv->FluxID   << "\n";

	cout << "Intr. mode     :" << nv->Mode   << "\n";

	cout << "# of particles =" << nv->Npart() << "\n";	

	for ( i = 0 ; i < nv->Npart() ; i++ ){
	  cout << "i=" << i << "\n";
	  cout << "Vertex         =" << nv->VertexID(i) << "\n";
	  cout << "Parent Index   =" << nv->ParentIdx(i) << "\n";

	  cout << "Particle Code  = " << (nv->PartInfo(i))->fPID   << "\n";
	  cout << "Particle Mass  = " << (nv->PartInfo(i))->fMass   << "\n";
	  cout << "Particle Mom.  =(" << (nv->PartInfo(i))->fP.Px() << ","
	  	   << (nv->PartInfo(i))->fP.Py() << "," 
		   << (nv->PartInfo(i))->fP.Pz() << ","
		   << (nv->PartInfo(i))->fP.E()  << ")"
		   << "\n";
	  cout << "Particle Flag  = " << (nv->PartInfo(i))->fIsAlive << "\n";
	  cout << "Particle Stat. = " << (nv->PartInfo(i))->fStatus  << "\n";
	  cout << "Particle Pos(1)=(" << (nv->PartInfo(i))->fPosIni.X() << ","
	  	   << (nv->PartInfo(i))->fPosIni.Y() << "," 
		   << (nv->PartInfo(i))->fPosIni.Z() << ","
		   << (nv->PartInfo(i))->fPosIni.T()  << ")"
		   << "\n";
	  cout << "Particle Pos(2)=(" << (nv->PartInfo(i))->fPosFin.Px() << ","
	  	   << (nv->PartInfo(i))->fPosFin.Y() << "," 
		   << (nv->PartInfo(i))->fPosFin.Z() << ","
		   << (nv->PartInfo(i))->fPosFin.T()  << ")"
		   << "\n";
	}
  }
}
  
  
  
