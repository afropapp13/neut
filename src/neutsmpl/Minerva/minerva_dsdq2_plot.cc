minerva_dsdq2_plot()
{
  
//                  0.0,    -0.025, -0.05, -0.10,-0.20,-0.40,-0.80,-1.2,-2.0};
  Double_t q2_exptbl[8] 
	= { 0.0125, 0.0375, 0.075, 0.15, 0.30, 0.60, 1.0, 1.6};
  Double_t q2_experr[8]
	= { 0.0125, 0.0125, 0.025, 0.05, 0.10, 0.20, 0.2, 0.4};
  
  
  Double_t dsdq2_numu[8]
	= {0.761,1.146,1.343,1.490,1.063,0.582,0.242,0.097};
  Double_t dsdq2_err1_numu[8]
	= {0.035,0.047,0.034,0.028,0.019,0.013,0.014,0.008};
  Double_t dsdq2_err2_numu[8] 
	= {0.097,0.137,0.156,0.170,0.120,0.074,0.053,0.024};

  Double_t dsdq2_numubar[8]
	= {0.813,1.061,1.185,1.096,0.777,0.340,0.123,0.041};
  Double_t dsdq2_err1_numubar[8]
	= {0.035,0.045,0.033,0.024,0.016,0.009,0.009,0.004};
  Double_t dsdq2_err2_numubar[8]
	= {0.102,0.134,0.150,0.135,0.101,0.050,0.024,0.010};
  
  Double_t dsdq2_err_numu[8];
  Double_t dsdq2_err_numubar[8];

  Int_t i;
  for ( i = 0 ; i < 8 ; i++ ){
	dsdq2_err_numu[i]    = dsdq2_err1_numu[i]    + dsdq2_err2_numu[i];
	dsdq2_err_numubar[i] = dsdq2_err1_numubar[i] + dsdq2_err2_numubar[i];
  }



  TGraphErrors *data_numu, *data_numubar;

  data_numu    = new TGraphErrors(8,q2_exptbl,dsdq2_numu,
								  q2_experr,dsdq2_err2_numu);
							   
  data_numubar = new TGraphErrors(8,q2_exptbl,dsdq2_numubar,
								  q2_experr,dsdq2_err_numubar);


  TFile *f_numu_sf, *f_numubar_sf;
  TFile *f_numu_rfg,*f_numubar_rfg;

  f_numu_sf     = new TFile("sf/hists_neut_minerva_numu.root");
  f_numubar_sf  = new TFile("sf/hists_neut_minerva_numubar.root");

  TH1D *h_q2_numu_sf,  *h_q2_numubar_sf;
  TH1D *h_q2_ccqe_numu_sf,  *h_q2_ccqe_numubar_sf;

  h_q2_numu_sf     = (TH1D *)f_numu_sf->Get("q2rec_all");
  h_q2_numubar_sf  = (TH1D *)f_numubar_sf->Get("q2rec_all");
  h_q2_ccqe_numu_sf     = (TH1D *)f_numu_sf->Get("q2rec_ccqe");
  h_q2_ccqe_numubar_sf  = (TH1D *)f_numubar_sf->Get("q2rec_ccqe");

  TCanvas *c_numu,*c_numubar;
  c_numu    = new TCanvas("q2_numu",   "q^2 for numu");
  c_numubar = new TCanvas("q2_numubar","q^2 for numubar");

  c_numu->SetTitle("q2 for numu");
  c_numu->cd();

  data_numu->SetLineColor(kBlack);
  data_numu->SetLineWidth(2);
  data_numu->SetMarkerColor(kBlack);
  data_numu->SetDrawOption("ap");

  h_q2_numu_sf->SetTitle("q2 for numu");
  h_q2_numu_sf->GetXaxis()->SetTitle("q^2_{rec} (GeV/c)^2");
  h_q2_numu_sf->SetAxisRange(0., 3.);
  h_q2_numu_sf->SetLineColor(kBlue);
  h_q2_numu_sf->SetFillColor(kBlue);
  h_q2_numu_sf->SetLineWidth(3);
  h_q2_numu_sf->Draw("");

  h_q2_ccqe_numu_sf->SetTitle("q2 for numu");
  h_q2_ccqe_numu_sf->SetAxisRange(0., 3.);
  h_q2_ccqe_numu_sf->SetLineColor(kGreen);
  h_q2_ccqe_numu_sf->SetFillColor(kGreen);
  h_q2_ccqe_numu_sf->SetLineWidth(3);
  h_q2_ccqe_numu_sf->Draw("SAME");

  data_numu->SetLineColor(kRed);
  data_numu->SetLineWidth(2);
  data_numu->SetMarkerColor(kRed);
  data_numu->SetDrawOption("p");
  data_numu->Draw("p same");

  c_numu->Print("numu_q2.pdf");

  c_numubar->SetTitle("q2 for numubar");
  c_numubar->cd();

  h_q2_numubar_sf->SetTitle("q2 for numubar");
  h_q2_numubar_sf->GetXaxis()->SetTitle("q^2_{rec} (GeV/c)^2");
  h_q2_numubar_sf->SetAxisRange(0., 3.);
  h_q2_numubar_sf->SetLineColor(kBlue);
  h_q2_numubar_sf->SetFillColor(kBlue);
  h_q2_numubar_sf->SetLineWidth(3);
  h_q2_numubar_sf->Draw("");

  h_q2_ccqe_numubar_sf->SetTitle("q2 for numubar");
  h_q2_ccqe_numubar_sf->SetAxisRange(0., 3.);
  h_q2_ccqe_numubar_sf->SetLineColor(kGreen);
  h_q2_ccqe_numubar_sf->SetFillColor(kGreen);
  h_q2_ccqe_numubar_sf->SetLineWidth(3);
  h_q2_ccqe_numubar_sf->Draw("SAME");

  data_numubar->SetLineColor(kRed);
  data_numubar->SetLineWidth(2);
  data_numubar->SetMarkerColor(kRed);
  data_numubar->SetDrawOption("p");

  data_numubar->Draw("p same");

  c_numubar->Print("numubar_q2.pdf");

}

