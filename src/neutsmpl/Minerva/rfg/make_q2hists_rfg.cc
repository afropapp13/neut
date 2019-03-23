#include <iostream>

using namespace std;

void
fill_q2hists(char *, char *, Double_t );


void
make_q2hists_rfg()
{

  Int_t i, j;
  Int_t nevents;

  gSystem->Load("../../../neutclass/neutfsipart.so");
  gSystem->Load("../../../neutclass/neutfsivert.so");
  gSystem->Load("../../../neutclass/neutvtx.so");
  gSystem->Load("../../../neutclass/neutpart.so");
  gSystem->Load("../../../neutclass/neutvect.so");
  
  TTree *tnv;
  NeutVect *nv;

  /* numu */
  /* flux : neutrinos / cm^2 / 10^-8 */
  /* cross-section   cm^2 / nucleon / 10^-38 */

  Double_t num_numu, rate_numu;

  Double_t num_numubar, rate_numubar;

  Double_t numnu, rate;

  TFile *f_hists_numu,*f_hists_numubar;

  TH1D *h_flux_numu,    *h_rate_numu;
  TH1D *h_flux_numubar, *h_rate_numubar;

  f_hists_numu    
	= new TFile("minerva_flux_fine_numu_o.root");
  f_hists_numubar 
	= new TFile("minerva_flux_fine_numubar_o.root");

  h_flux_numu    = (TH1D *)f_hists_numu->Get("flux_numu");
  h_rate_numu    = (TH1D *)f_hists_numu->Get("evtrt_numu");

  h_flux_numubar = (TH1D *)f_hists_numubar->Get("flux_numub");
  h_rate_numubar = (TH1D *)f_hists_numubar->Get("evtrt_numub");

  num_numu      = h_flux_numu->GetSumOfWeights();
  rate_numu     = h_rate_numu->GetSumOfWeights();

  num_numubar   = h_flux_numubar->GetSumOfWeights();
  rate_numubar  = h_rate_numubar->GetSumOfWeights();

  rate_numu    = rate_numu / num_numu / 100000.;
  rate_numu    = rate_numu * 13./6. ;

  rate_numubar = rate_numubar / num_numubar / 100000.;
  rate_numubar = rate_numubar * 13./7. ;

  fill_q2hists("neut_minerva_numu_rfg.root",
			   "hists_neut_minerva_numu_rfg.root",
			   rate_numu);

  fill_q2hists("neut_minerva_numubar_rfg.root",
			   "hists_neut_minerva_numubar_rfg.root",
			   rate_numubar);
}

void
fill_q2hists(char *in_fname, char *out_fname, Double_t rate)
{
  Int_t j;
  Int_t nevents;

  TFile *f;
  f = new TFile(in_fname,"READ");
  if ( f == NULL ){
	cout << "Failed to open " << in_fname << endl;
	return;
  }
  tnv = (TTree *)(f->Get("neuttree"));
  nv = new NeutVect();

  tnv->SetBranchAddress("vectorbranch",&nv);

  nevents = tnv->GetEntries();


  TFile *f2;
  f2 = new TFile(out_fname,"RECREATE");
  if ( f2 == NULL ){
	cout << "Failed to (re)create " << out_fname << endl;
	return;
  }

  TH1D *h_p_mu_ccqe, *h_angle_mu_ccqe, *h_q2_ccqe, *h_enurec_ccqe, *h_q2rec_ccqe;
  TH1D *h_p_mu_mec, *h_angle_mu_mec, *h_q2_mec *h_enurec_mec, *h_q2rec_mec;
  TH1D *h_p_mu_all, *h_angle_mu_all, *h_q2_all *h_enurec_all, *h_q2rec_all;

  h_p_mu_ccqe     = new TH1D("pmu_ccqe","muon momentum ( CCQE )", 200,0.,10.);
  h_angle_mu_ccqe = new TH1D("anglemu_ccqe","muon angle ( CCQE )", 18,0.,180.);
  h_q2_ccqe       = new TH1D("q2_ccqe","q2 ( CCQE )",             200,0.,10.);
  h_enurec_ccqe   = new TH1D("enurec_ccqe","enu rec( CCQE )",      18,0.,180.);
  h_q2rec_ccqe    = new TH1D("q2rec_ccqe","q2 rec( CCQE )",       200,0.,10.);

  h_p_mu_mec     = new TH1D("pmu_mec","muon momentum ( MEC ) ",  200,0.,10.);
  h_angle_mu_mec = new TH1D("anglemu_mec","muon angle ( MEC )",   18,0.,180.);
  h_q2_mec       = new TH1D("q2_mec","q2 ( MEC )",               200,0.,10.);
  h_enurec_mec   = new TH1D("enurec_mec","enu rec( MEC )",        18,0.,180.);
  h_q2rec_mec    = new TH1D("q2rec_mec","q2 rec( MEC )",         200,0.,10.);

  h_p_mu_all     = new TH1D("pmu_all","muon momentum ( ALL )",  200,0.,10.);
  h_angle_mu_all = new TH1D("anglemu_all","muon angle ( ALL )",  18,0.,180.);
  h_q2_all       = new TH1D("q2_all","q2 ( ALL )",              200,0.,10.);
  h_enurec_all   = new TH1D("enurec_all","enu rec( ALL )",       18,0.,180.);
  h_q2rec_all    = new TH1D("q2rec_all","q2 rec( ALL )",        200,0.,10.);


  Double_t emu, pmu, angle_mu,q2, cos_mu;
  Double_t xmn, xmmu, eb;
  Double_t enu_rec, q2_rec;

  Int_t    loc;

  for ( j = 0 ; j < nevents ; j++ ){
	loc = 0;
	tnv->GetEntry(j);
	//	cout << "Intr. mode     :" << nv->Mode   << "\n";

	if (abs(nv->Mode) == 1){
	  loc = 2;
	}
	if (abs(nv->Mode) == 2){
	  loc = 3;
	}

	if (loc > 0){
	  //	  cout << "Parent Index   =" << nv->ParentIdx(i) << "\n";
	  //	  cout << "Particle Code  = " << (nv->PartInfo(i))->fPID   << "\n";
	  //	  cout << "Particle Mass  = " << (nv->PartInfo(i))->fMass   << "\n";
	  //	  cout << "Particle Mom.  =(" << (nv->PartInfo(i))->fP.Px() << ","
	  //		   << (nv->PartInfo(i))->fP.Py() << "," 
	  //		   << (nv->PartInfo(i))->fP.Pz() << ","
	  //		   << (nv->PartInfo(i))->fP.E()  << ")"
	  //		   << "\n";

	  emu = (nv->PartInfo(loc))->fP.E()/1000.;
	  pmu = (nv->PartInfo(loc))->fP.P()/1000.;
	  xmn  = 0.939;
	  xmmu = 0.106;
	  eb   =-0.025;

	  angle_mu = (nv->PartInfo(0))->fP.Angle((nv->PartInfo(loc))->fP.Vect());
	  angle_mu = angle_mu / 3.1415926535 * 180.;

	  q2 = 
		((nv->PartInfo(0))->fP - (nv->PartInfo(loc))->fP) *
		((nv->PartInfo(0))->fP - (nv->PartInfo(loc))->fP);
	  q2 = q2 * (-1.e-6);

	  cos_mu 
		= cos((nv->PartInfo(0))->fP.Angle((nv->PartInfo(loc))->fP.Vect()));
	  enu_rec
		= ( (xmn + eb )*emu - ( 2*xmn*eb+eb*eb+xmmu*xmmu) / 2 )
		/ ( xmn + eb - emu + pmu*cos_mu );
	  q2_rec = (2*enu_rec*(emu - pmu*cos_mu)-xmmu*xmmu);
	  
	  h_p_mu_all->Fill(pmu);
	  h_angle_mu_all->Fill(angle_mu);
	  h_q2_all->Fill(q2,rate/0.050);

	  h_enurec_all->Fill(enu_rec);
	  h_q2rec_all->Fill(q2_rec,rate/0.050);

	  if (abs(nv->Mode) == 1){
		h_p_mu_ccqe->Fill(pmu);
		h_angle_mu_ccqe->Fill(angle_mu);
		h_q2_ccqe->Fill(q2,rate/0.050);
		h_enurec_ccqe->Fill(enu_rec);
		h_q2rec_ccqe->Fill(q2_rec,rate/0.050);

		/*
		cout << "E_nu,   q2, angle    =(" 
			 << (nv->PartInfo(0))->fP.E()*1.e-3 << ","
			 << q2                     << "," 
			 << angle_mu               << ")"
			 << endl;
		
		cout << "enurec, q2rec, cos_mu =(" 
			 << enu_rec                << ","
			 << q2_rec                 << ")" 
			 << cos_mu                 << ")"
			 << endl;
		cout << "-" << endl;
		*/
		
	  }

	  if (abs(nv->Mode) == 2){
		h_p_mu_mec->Fill(pmu);
		h_angle_mu_mec->Fill(angle_mu);
		h_q2_mec->Fill(q2,rate/0.050);
		h_enurec_mec->Fill(enu_rec);
		h_q2rec_mec->Fill(q2_rec,rate/0.050);
	  }

	  /*
	  cout << "E_nu,   q2, angle =(" 
		   << (nv->PartInfo(i))->fP.E()*1.e-3 << ","
		   << q2                     << "," 
		   << angle_mu               << ")"
		   << endl;

	  cout << "enurec, q2rec     =(" 
		   << enu_rec                << ","
		   << q2rec                  << ")" 
		   << endl;
	  */
	}
  }
  f2->Write();
  f2->Close();

}
    
