void
chkfillneutvtx()
{

  Int_t i, j;

  gSystem->Load("neutvtx.so");
  
  NeutVtx *nv;

  nv = new NeutVtx();

  TFile f("neutvtx.root","recreate");
  TTree tarr("NeutVrtex", "Neut Vertex info.");

  tarr.Branch("Vertex","NeutVtx",&nv,1024,99);

  for ( j = 0 ; j < 100 ; j++ ){

	nv->EventNo = j;

	std::cout << "Event #" << EventNo << "\n";

	nv->SetNvtx(0);

	Int_t nvtx = (Int_t)(gRandom->Rndm()*30+1.0);

	nv->SetNvtx(nvtx);

	cout << "# of vertexes =" << nv->Nvtx() << "\n";	

	for ( i = 0 ; i < nvtx ; i++ ){
	  cout << "i=" << i << "\n";

	  TLorentzVector pinfo;

	  pinfo.SetXYZT(gRandom->Rndm()*200-100.,
					gRandom->Rndm()*200-100.,
					gRandom->Rndm()*200-100.,
					0.);
	  
	  nv->SetPos(i, pinfo);
	  
	  cout << "Vertex Pos(1)=(" << (nv->Pos(i))->X() << ","
	  	   << (nv->Pos(i))->Y() << "," 
		   << (nv->Pos(i))->Z() << ","
		   << (nv->Pos(i))->T()  << ")"
		   << "\n";
	}
	tarr.Fill();  
  }
  f.Write();

}
  
  
  
