void
chkfillneutvect()
{

  Int_t i, j;

  gSystem->Load("neutpart.so");
  gSystem->Load("neutvect.so");
  
  NeutVect *nv;

  nv = new NeutVect();

  TFile f("neutvect.root","recreate");
  TTree tarr("NeutVector", "Neut Vector info.");

  tarr.Branch("Vector","NeutVect",&nv,1024,99);

  NeutPart pinfo;

  for ( j = 0 ; j < 100 ; j++ ){

	nv->SetNpart(0);

	/****************************************************/
	nv->EventNo = j;
	cout << "Event #        :" << nv->EventNo << "\n";

	if (gRandom->Rndm()>0.11){
	  nv->TargetA  =   8;
	  nv->TargetZ  =   8;
	  nv->VNuclIni = -27.;
	  nv->VNuclFin =   0.;
	  nv->PFSurf   = 236.;
	  nv->PFMax    = 236.;
	}else{
	  nv->TargetA = 1;
	  nv->TargetZ = 1;
	  nv->VNuclIni =   0.;
	  nv->VNuclFin =   0.;
	  nv->PFSurf   =   0.;
	  nv->PFMax    =   0.;
	}
	nv->FluxID     = 0;

	cout << "Target A       :" << nv->TargetA << "\n";
	cout << "Target Z       :" << nv->TargetA << "\n";
	cout << "VNuclIni       :" << nv->VNuclIni << "\n";
	cout << "VNuclFin       :" << nv->VNuclFin << "\n";
	cout << "PF Surface     :" << nv->PFSurf   << "\n";
	cout << "PF Maximum     :" << nv->PFMax    << "\n";
	cout << "Flux ID        :" << nv->FluxID   << "\n";

	/****************************************************/
	nv->Mode = (Int_t)(gRandom->Rndm()*30-15.);

	cout << "Intr. mode     :" << nv->Mode   << "\n";

	Int_t np = (Int_t)(gRandom->Rndm()*30+1.0);

	nv->SetNpart(np);

	cout << "# of particles =" << nv->Npart() << "\n";	

	for ( i = 0 ; i < np ; i++ ){
	  cout << "i=" << i << "\n";
	  Double_t mass;

	  nv->SetVertexID(i,1);

	  cout << "Vertex         =" << nv->VertexID(i) << "\n";

	  if (i < 2){
		nv->SetParentIdx(i,-1);
	  }else{
		nv->SetParentIdx(i,1);
	  }

	  cout << "Parent Index   =" << nv->ParentIdx(i) << "\n";
	  
	  pinfo.fPID = (Int_t)(gRandom->Rndm()*1000);
	  mass = (gRandom->Rndm()*1000);
	  pinfo.fMass = mass;
	  if (i < 2){
		pinfo.fIsAlive = false;
		pinfo.fStatus  = -1;
	  }else{
		pinfo.fIsAlive = true;
		pinfo.fStatus  = 0;
	  }

	  pinfo.fP.SetXYZM(gRandom->Rndm()*1000,
					   gRandom->Rndm()*1000,
					   gRandom->Rndm()*1000,
					   mass);

	  pinfo.fPosIni.SetXYZT(gRandom->Rndm()*2-1.,
							gRandom->Rndm()*2-1.,
							gRandom->Rndm()*2-1.,
							0.);
	  
	  pinfo.fPosFin.SetXYZT(gRandom->Rndm()*2-1.,
							gRandom->Rndm()*2-1.,
							gRandom->Rndm()*2-1.,
							0.);
	  
	  nv->SetPartInfo(i, pinfo);
	  
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
	  cout << "Particle Pos(2)=(" << (nv->PartInfo(i))->fPosFin.X() << ","
	  	   << (nv->PartInfo(i))->fPosFin.Y() << "," 
		   << (nv->PartInfo(i))->fPosFin.Z() << ","
		   << (nv->PartInfo(i))->fPosFin.T()  << ")"
		   << "\n";
	}
	tarr.Fill();  
  }
  f.Write();

}
  
  
  
