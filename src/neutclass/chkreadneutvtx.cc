void
chkreadneutvtx()
{

  Int_t i, j;
  Int_t nevents;

  gSystem->Load("neutvtx.so");
  
  TTree  *tnv;
  NeutVtx *nv;

  TFile f("neutvtx.root");
  tnv = (TTree *)(f.Get("NeutVrtex"));

  nv = new NeutVtx();
  tnv->SetBranchAddress("Vertex",&nv);

  nevents = tnv->GetEntries();

  for ( j = 0 ; j < nevents ; j++ ){

	tnv->GetEntry(j);

	std::cout << "Event #" << nv->EventNo << "\n";

	cout << "# of vertexes =" << nv->Nvtx() << "\n";	
	for (i = 0 ; i < nv->Nvtx() ; i++){
	  cout << "i=" << i << "\n";
	  
	  cout << "Vertex Pos(1)=(" << (nv->Pos(i))->X() << ","
		   << (nv->Pos(i))->Y() << "," 
		   << (nv->Pos(i))->Z() << ","
		   << (nv->Pos(i))->T()  << ")"
		   << "\n";
	}
  }
}
  
  
  
