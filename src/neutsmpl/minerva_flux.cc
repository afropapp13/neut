{
  Double_t numu[17]     = {0.310,0.409,0.504,0.526,0.423,0.253,0.137,0.081,
						   0.055,0.043,0.036,0.031,0.027,0.024,0.021,0.019,
						   0.017};
  Double_t antinumu[17] = {0.281,0.368,0.444,0.448,0.349,0.205,0.106,0.061,
						   0.038,0.029,0.022,0.018,0.016,0.013,0.012,0.010,
						   0.009};

  TFile f("Minerva_flux.root","RECREATE");

  TH1D *minervaflux_tbl[2];

  minervaflux_tbl[0] = new TH1D("minerva_flux_numu","Minerva numu flux",
								17,1.5,10.0);
  minervaflux_tbl[1] = new TH1D("minerva_flux_numubar","Minerva numu bar flux",
								17,1.5,10.0);
  Int_t i;
  for ( i = 1 ; i <= 17 ; i++ ){
	minervaflux_tbl[0]->SetBinContent(i,numu[i-1]);
	minervaflux_tbl[1]->SetBinContent(i,antinumu[i-1]);
  }

  TSpline3 *minervaflux_spline[2];

  for ( i = 0 ; i < 2 ; i++ ){
	minervaflux_spline[i] = new TSpline3(minervaflux_tbl[i]);
  }

  TH1D *minervaflux_fine[2];
  minervaflux_fine[0] = new TH1D("minerva_flux_fine_numu",
								"Minerva numu fine flux",85,1.5,10.0);
  minervaflux_fine[1] = new TH1D("minerva_flux_fine_numubar",
								"Minerva numu bar fine flux",85,1.5,10.0);

  Double_t energy;
  for ( i = 1 ; i <= 85 ; i++ ){
	energy = 1.45 + ((Double_t)(i)) * 0.1;
	minervaflux_fine[0]->SetBinContent(i,minervaflux_spline[0]->Eval(energy));
	minervaflux_fine[1]->SetBinContent(i,minervaflux_spline[1]->Eval(energy));
  }

  f->Write();
  f->Close();
  
}

