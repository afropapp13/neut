 //STL
#include <vector>
#include <iostream>
#include <string>
#include <sstream>

#include <fstream>
#include <set>
#include <cstdlib>


//ROOT
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TMath.h"
#include "TAxis.h"

void ResetStyle(TLegend *obj, Double_t mar, const Double_t ff)
{
  if(mar>0){
    obj->SetMargin(mar);
  }
  obj->SetFillStyle(-1);
  obj->SetBorderSize(-1);

}

void checknscat(){
  gStyle->SetOptStat(kFALSE);


  //  TFile file("../checknscat2/neut536/inputs3/pc_nom_test.root"); 
  TFile file("Outputnuc360nrprd300.root");

  TTree *h1 = (TTree*)file.Get("h1");

  Int_t           Npvc;
  Int_t           Ipvc[100];   //[Npvc] 
  Int_t           Iorgvc[100];   //[Npvc]
  Int_t           Ichvc[100];   //[Npvc]
  Int_t           Iflvc[100];   //[Npvc]
  Float_t         Abspvc[100];   //[Npvc]
  Float_t         Pvc[100][3];   //[Npvc]
  Float_t         Posvc[3];
  Int_t           Nvert;
  Int_t           Nfnvert;
  Float_t         Posvert[150][3];   //[Nvert]
  Int_t           Iflgvert[150];   //[Nvert]
  Int_t           Nvcvert;
  Float_t         Dirvert[900][3];   //[Nvcvert]
  Float_t         Abspvert[900];   //[Nvcvert]
  Float_t         Abstpvert[900];   //[Nvcvert]
  Int_t           Ipvert[900];   //[Nvcvert]
  Int_t           Iverti[900];   //[Nvcvert]
  Int_t           Ivertf[900];   //[Nvcvert]
  Float_t         Fsiprob;
  Int_t           Numbndn;
  Int_t           Numbndp;
  Int_t           Numfrep;
  Int_t           Numatom;

  TBranch        *b_Npvc;   //!
  TBranch        *b_Ipvc;   //!
  TBranch        *b_Iorgvc;   //!
  TBranch        *b_Ichvc;   //!
  TBranch        *b_Iflvc;   //!
  TBranch        *b_Abspvc;   //!
  TBranch        *b_Pvc;   //!
  TBranch        *b_Posvc;   //!
  TBranch        *b_Nvert;   //!
  TBranch        *b_Nfnvert;
  TBranch        *b_Posvert;   //!
  TBranch        *b_Iflgvert;   //!
  TBranch        *b_Nvcvert;   //!
  TBranch        *b_Dirvert;   //!
  TBranch        *b_Abspvert;   //!
  TBranch        *b_Abstpvert;   //!
  TBranch        *b_Ipvert;   //!
  TBranch        *b_Iverti;   //!
  TBranch        *b_Ivertf;   //!
  TBranch        *b_Fsiprob;   //!
  TBranch        *b_Numbndn;   //!
  TBranch        *b_Numbndp;   //!
  TBranch        *b_Numfrep;   //!
  TBranch        *b_Numatom;   //!

  h1->SetBranchAddress("Npvc", &Npvc, &b_Npvc);
  h1->SetBranchAddress("Ipvc", Ipvc, &b_Ipvc);
  h1->SetBranchAddress("Iorgvc", Iorgvc, &b_Iorgvc);
  h1->SetBranchAddress("Ichvc", Ichvc, &b_Ichvc);
  h1->SetBranchAddress("Iflvc", Iflvc, &b_Iflvc);
  h1->SetBranchAddress("Abspvc", Abspvc, &b_Abspvc);
  h1->SetBranchAddress("Pvc", Pvc, &b_Pvc);
  h1->SetBranchAddress("Posvc", Posvc, &b_Posvc);
  h1->SetBranchAddress("Nvert", &Nvert, &b_Nvert);
  h1->SetBranchAddress("Nfnvert", &Nfnvert, &b_Nfnvert);
  h1->SetBranchAddress("Posvert", Posvert, &b_Posvert);
  h1->SetBranchAddress("Iflgvert", Iflgvert, &b_Iflgvert);
  h1->SetBranchAddress("Nvcvert", &Nvcvert, &b_Nvcvert);
  h1->SetBranchAddress("Dirvert", Dirvert, &b_Dirvert);
  h1->SetBranchAddress("Abspvert", Abspvert, &b_Abspvert);
  h1->SetBranchAddress("Abstpvert", Abstpvert, &b_Abstpvert);
  h1->SetBranchAddress("Ipvert", Ipvert, &b_Ipvert);
  h1->SetBranchAddress("Iverti", Iverti, &b_Iverti);
  h1->SetBranchAddress("Ivertf", Ivertf, &b_Ivertf);
  h1->SetBranchAddress("Fsiprob", &Fsiprob, &b_Fsiprob);
  h1->SetBranchAddress("Numbndn", &Numbndn, &b_Numbndn);
  h1->SetBranchAddress("Numbndp", &Numbndp, &b_Numbndp);
  h1->SetBranchAddress("Numfrep", &Numfrep, &b_Numfrep);
  h1->SetBranchAddress("Numatom", &Numatom, &b_Numatom);
  //basicTreeChain->SetBranchAddress("BasicHeader", &BasicHeader);

  int nentries = h1->GetEntries();
  
  int max_bin = 3000;
  int n_bin = max_bin/100; //max_bin/10
  
  
  //int max_bin = 300;
  //int n_bin = max_bin/300;
  int a_nBins = 90;
  float a_xMin = 0;
  float a_xMax = 180; // degrees
  double dtheta = (a_xMax-a_xMin)/a_nBins * TMath::Pi()/180.;

  TH1F* all = new TH1F("all", "", n_bin, 0., max_bin);
  TH1F* tot = new TH1F("tot", "", n_bin, 0., max_bin);
  TH1F* number_particles = new TH1F("number_particles", "", 30, 0., 30);
  TH1F* final_tot = new TH1F("final_tot", "", n_bin, 0., max_bin);
  TH1F* final_tot_p = new TH1F("final_tot_p", "", n_bin, 0., max_bin);
  TH1F* final_tot_n = new TH1F("final_tot_n", "", n_bin, 0., max_bin);
  TH1F* final_tot_pim = new TH1F("final_tot_pim", "", n_bin, 0., max_bin);
  TH1F* final_tot_pi0 = new TH1F("final_tot_pi0", "", n_bin, 0., max_bin);
  TH1F* final_tot_pip = new TH1F("final_tot_pip", "", n_bin, 0., max_bin);
  TH1F* Nvertices = new TH1F("Nvertices", "", 40, 0., 40);
  TH1F* reac = new TH1F("reac", "", n_bin, 0., max_bin);
  TH1F* survial = new TH1F("survial", "   ", n_bin, 0., max_bin);
  TH1F* QuasiElas = new TH1F("QuasiElas", "   ", n_bin, 0., max_bin);
  TH1F* Elastic = new TH1F("Elastic", "   ", n_bin, 0., max_bin);
  TH1F* Spiprod = new TH1F("Spiprod", "   ", n_bin, 0., max_bin);
  TH1F* Dpiprod = new TH1F("Dpiprod", "   ", n_bin, 0., max_bin);
  TH1F* Morepiprod = new TH1F("Morepiprod", "   ", n_bin, 0., max_bin);
  TH1F* InitialMom = new TH1F("InitialMom", "   ", n_bin, 0., max_bin);

  TH1F* diffxsec = new TH1F("diffxsec", "   ", a_nBins, a_xMin, a_xMax);
  TH1F* multiplicity = new TH1F("multiplicity", "   ", 20, 0, 20);
  TH2F* outphasespc = new TH2F("outphasespc", "   ", n_bin/5, 0., max_bin, 100, -1, 1.);
  int NuclSurv = 0;
  int NuclDie = 0;
  int PiSurv = 0;
  int PiDie = 0;
  int counttot = 0;
  int countelas = 0;
  int countinelas = 0;
  int countothers = 0;
  int countdiff = 0;
  double outgoingp = 0.;
  float weight_ang;
  double momdiff = 0.;

  for (int i = 0; i < nentries; i++){
    Nvertices->Fill(Nfnvert);    
    //    std::cout << "NVERT =  " << Nvert << std::endl;
    if(!(i % 100000)){std::cout << "checked " << i << "/" << nentries << std::endl;}

	  h1->GetEntry(i);
	  // std::cout << i << "   " << Ipvc[4] << std::endl;


	  //std::cout << "***" << std::endl;
	  //std::cout << mode << std::endl;

	  //std::cout << Npvc << std::endl;
	  //std::cout << Nvert << std::endl;
//for (int z = 0; z < 100; ++z) {
//std::cout << z << std::endl;
//std::cout << "    " << Ipvc[z] << std::endl;
//std::cout << "    " << Iorgvc[z] << std::endl;
//std::cout << "    " << Ichvc[z] << std::endl;
//std::cout << "    " << Iflvc[z] << std::endl;
//std::cout << "    " << Abspvc[z] << std::endl;
//for (int z2 = 0; z2 < 3; ++z2) {
//std::cout << "    " << Pvc[z][z2] << std::endl;
//}
//}

//for (int z2 = 0; z2 < 3; ++z2) {
//std::cout << Posvc[z2] << std::endl;
//}
//std::cout << Posvert[150][3] << std::endl;
//std::cout << Iflgvert[150] << std::endl;
//std::cout << Nvcvert << std::endl;
//std::cout << Dirvert[900][3] << std::endl;
//std::cout << Abspvert[900] << std::endl;
//std::cout << Abstpvert[900] << std::endl;
//std::cout << Ipvert[900] << std::endl;
//std::cout << Iverti[900] << std::endl;
//std::cout << Ivertf[900] << std::endl;
//std::cout << Fsiprob << std::endl;
//std::cout << Numbndn << std::endl;
//std::cout << Numbndp << std::endl;
//std::cout << Numfrep << std::endl;
//std::cout << Numatom << std::endl;

	  Float_t evt_weight;
	  evt_weight = 1;

	  all->Fill(Abspvc[3],evt_weight); //Absolute momentum of 3rd particle
	  //find outgoing momentum
	  //they live at different places in the array!!
	  // int ipart=6;
	  // while(Ipvc[ipart]!=2212) ipart++;
	  // outgoingp = Abspvc[ipart];
	  
	  //for(int x =0; x<Npvc;x++)
	  //{std::cout << "particle: " << x << " total momentum: "<< Abspvc[x] << std::endl;
	  //} 
	  
	  if (Npvc == 5) {
	    
	    //	    std::cout << "particle 3 = " << Ipvc[3] << "particle 4 = " << Ipvc[4] << " particle 5 = " << Ipvc[5] << std::endl <<"particle 5 momentum = " << Abspvc[5]<<"particle 1 momentum = " << Abspvc[1]<<"particle 2 momentum = " << Abspvc[2]<< "particle 3 momentum = " << Abspvc[3] << "particle 4 momentum = " << Abspvc[4] << std::endl;
	    //std::cout << "  particle 3 X mom = " << Pvc[3][0] << "  particle 4 X mom = " << Pvc[4][0] << "  particle 3 Y mom = " << Pvc[3][1] << "  particle 4 Y mom = " << Pvc[4][1] << "  particle 3 Z mom = " << Pvc[3][2] << "  particle 4 Z mom = " << Pvc[4][2] << std::endl 
	    //	    std::cout << "  particle 3 Ichvc " << Ichvc[3] << "  particle 4 Ichvc " << Ichvc[4]  << "  particle 3 Iorgvc " << Iorgvc[3] << "  particle 4 Iorgvc " << Iorgvc[4]<< "  particle 3 Iflvc " << Iflvc[3]<< "  particle 4 Iflvc " << Iflvc[4] << std::endl << std::endl;






	  }
	  //	  Nvertices->Fill(Nvert);  
	  number_particles->Fill(Npvc);
	  // survial / Normal
	  if (Npvc== 4) {  //if there are 4 particles in total, (ingoing and outgoing)
		  // if (Iflvc[3]==0) {

		  survial->Fill(Abspvc[3],evt_weight);


	  } else {
	    //Nvertices->Fill(Nvert);
	    tot->Fill(Abspvc[3],evt_weight);
	    InitialMom->Fill(Abspvc[3],evt_weight);
	    counttot++;
	    //	    final_tot->Fill(Abspvc[4],evt_weight);
	   
	    



		  //  counttot++;

		  // if (Abspvc[3]==Abspvc[4]) std::cout << "ELASTIC!" << std::endl;
		  int particle=-1;
		  int PID=-1;
		  int nnucl=0;
		  int nType[5];
		  nType[0]=nType[1]=nType[2]=nType[3]=nType[4]=0;

		  float cost, px[3], pabs;
		  // Determine final state
		  // std::cout <<"------------------------------------------" << std::endl;
		  //number_particles->Fill(Npvc);
		  for (int ipvc=3; ipvc<Npvc; ipvc++) { //for particles beyond the initial particles
		    
		    
		    if (Ichvc[ipvc]!=1) continue;//check if it is detected
					  int cat;
					  cat++;	  
			  // Check the particle type
			  if (Ipvc[ipvc] == 2212) {particle=0;}//p
			  else if (Ipvc[ipvc] == 2112) {particle=1;}//n
			  else if (Ipvc[ipvc] == -211) {particle=2;}//pim
			  else if (Ipvc[ipvc] == 111) {particle=3;}//pi0
			  else if (Ipvc[ipvc] == 211) {particle=4;}//pip
			  else continue;

			  px[0] = Pvc[ipvc][0];
			  px[1] = Pvc[ipvc][1];
			  px[2] = Pvc[ipvc][2];
			  pabs = Abspvc[ipvc]; //save momentum values
			  cost = px[0]/pabs; //define cos t

			  PID=ipvc;  

			  nType[particle]++;      //count each particle type
			  nnucl++;
			  
			  if (ipvc > 3){ //if npvc>ipvc
			    
			    final_tot->Fill(Abspvc[ipvc],evt_weight);                                                                                                     
		
			    if (Ipvc[ipvc] == 2212) {final_tot_p->Fill(Abspvc[ipvc],evt_weight);}//p
			    else if (Ipvc[ipvc] == 2112) {final_tot_n->Fill(Abspvc[ipvc],evt_weight);}//n
			    else if (Ipvc[ipvc] == -211) {final_tot_pim->Fill(Abspvc[ipvc],evt_weight);}//pim
			    else if (Ipvc[ipvc] == 111) {final_tot_pi0->Fill(Abspvc[ipvc],evt_weight);}//pi0
			    else if (Ipvc[ipvc] == 211) {final_tot_pip->Fill(Abspvc[ipvc],evt_weight);}//pip
			    //final_tot->Fill(Abspvc[ipvc],evt_weight);
			  


			  }

		  }

		  if ((nType[0]+nType[1])!=1) multiplicity->Fill(nType[0]+nType[1]-1); //if more than one nucleon, fill multiplicity
		  float theta = acos(cost)*180/TMath::Pi(); //define theta from cos t
		  float sint = sqrt(1-cost*cost); //define sin t from cos t
		  if (sint==0) sint=0.000001; //remove 0 errors
		  weight_ang =  1/(sint*(2*TMath::Pi()*dtheta)); //provide weighted angle from dtheta (dependent on number of bins)

		  outphasespc->Fill(pabs,cost,evt_weight); //fill with absolute momentum and cos t to define the phase space
		  // Quasielasscat Scattering

		  //		  if (Npvc==6){
		
		  if(Npvc == 6){
		  QuasiElas->Fill(Abspvc[3],evt_weight);
			  // elasscatout->Fill(Abspvc[PID],evt_weight);//outgoing p for final nucleon
			  countelas++;

		  }
		  else if (Npvc>6){
			  reac->Fill(Abspvc[3],evt_weight);

		  }
		  // Hadronic Production
		  if (nType[2]+nType[3]+nType[4]==1){
			  Spiprod->Fill(Abspvc[3],evt_weight);        
		  }

		  else if (nType[2]+nType[3]+nType[4]==2){
			  Dpiprod->Fill(Abspvc[3],evt_weight);        
		  }         
		  else if (nType[2]+nType[3]+nType[4]>2){
			  Morepiprod->Fill(Abspvc[3],evt_weight);
		  }


	  }

	  }//end of event loop

	  TAxis *xaxis = tot->GetXaxis();
	  Int_t nbins  = xaxis->GetNbins();
	  Float_t xmin = xaxis->GetXmin();
	  Float_t xmax = xaxis->GetXmax();
	  std::cout << nbins << "  " << xmin << "  " << xmax << std::endl;

	  double eftarget = 2.355;//C
	  // double eftarget = 2.69;//O
	  double rMax = 2.5*eftarget;

	  double targetArea = TMath::Pi()*rMax*rMax*10.;  // mb
	  std::cout << "targetArea: " << targetArea << " " << nentries << std::endl;
	  std::cout << "counttot: " << counttot << std::endl;
	  std::cout << "countelas: " << countelas << std::endl;
	  std::cout << "countothers: " << countothers << std::endl;
	  std::cout <<"cat: " << cat      << std::endl;
	  TH1F *totxsec = new TH1F("totxsec","",nbins,xmin,xmax);
	  TH1F *final_totxsec = new TH1F("final_totxsec","",nbins,xmin,xmax);

	  TH1F *quasielaxsec = new TH1F("quasielaxsec","",nbins,xmin,xmax);
	  TH1F *elaxsec = new TH1F("elaxsec","",nbins,xmin,xmax);

	  TH1F *spixsec = new TH1F("spixsec","",nbins,xmin,xmax);
	  TH1F *dpixsec = new TH1F("dpixsec","",nbins,xmin,xmax);
	  TH1F *morepixsec = new TH1F("morepixsec","",nbins,xmin,xmax);
	  TH1F *initmomxsec = new TH1F("initmomxsec","",nbins,xmin,xmax);
	  
	  TH1F *reaxsec = new TH1F("reaxsec","",nbins,xmin,xmax);
	  TH1F *pdiffxsec = new TH1F("pdiffxsec","",a_nBins,a_xMin,a_xMax);

	  for (Int_t bin=0;bin<=nbins;bin++) {

		  if (all->GetBinContent(bin)==0) continue;

		  totxsec->SetBinContent(bin,(tot->GetBinContent(bin))*targetArea/(all->GetBinContent(bin)));
		 totxsec->SetBinError(bin,(tot->GetBinError(bin))*targetArea/(all->GetBinContent(bin)));

		  //final_totxsec->SetBinContent(bin,(final_tot->GetBinContent(bin))*targetArea/(all->GetBinContent(bin)));

		  quasielaxsec->SetBinContent(bin,(QuasiElas->GetBinContent(bin))*targetArea/(all->GetBinContent(bin)));
		  quasielaxsec->SetBinError(bin,(QuasiElas->GetBinError(bin))*targetArea/(all->GetBinContent(bin)));

		  elaxsec->SetBinContent(bin,(Elastic->GetBinContent(bin))*targetArea/(all->GetBinContent(bin)));
		  elaxsec->SetBinError(bin,(Elastic->GetBinError(bin))*targetArea/(all->GetBinContent(bin)));


		  reaxsec->SetBinContent(bin,(reac->GetBinContent(bin))*targetArea/(all->GetBinContent(bin)));
		  reaxsec->SetBinError(bin,(reac->GetBinError(bin))*targetArea/(all->GetBinContent(bin)));


		  spixsec->SetBinContent(bin,(Spiprod->GetBinContent(bin))*targetArea/(all->GetBinContent(bin)));
		  spixsec->SetBinError(bin,(Spiprod->GetBinError(bin))*targetArea/(all->GetBinContent(bin)));

		  dpixsec->SetBinContent(bin,(Dpiprod->GetBinContent(bin))*targetArea/(all->GetBinContent(bin)));
		  dpixsec->SetBinError(bin,(Dpiprod->GetBinError(bin))*targetArea/(all->GetBinContent(bin)));

		  morepixsec->SetBinContent(bin,(Morepiprod->GetBinContent(bin))*targetArea/(all->GetBinContent(bin)));
		  morepixsec->SetBinError(bin,(Morepiprod->GetBinError(bin))*targetArea/(all->GetBinContent(bin)));

		  initmomxsec->SetBinContent(bin,(InitialMom->GetBinContent(bin))*targetArea/(all->GetBinContent(bin)));

	  }




	  totxsec->SetEntries(tot->GetEntries());
	  //final_totxsec->SetEntries(final_tot->GetEntries());

	  quasielaxsec->SetEntries(QuasiElas->GetEntries());
	  elaxsec->SetEntries(Elastic->GetEntries());

	  reaxsec->SetEntries(reac->GetEntries());
	  spixsec->SetEntries(Spiprod->GetEntries());
	  dpixsec->SetEntries(Dpiprod->GetEntries());
	  morepixsec->SetEntries(Morepiprod->GetEntries());
	  initmomxsec->SetEntries(InitialMom->GetEntries());
	  pdiffxsec->SetEntries(diffxsec->GetEntries());




	  TFile* mmmm = new TFile("rootplots/OUTPUTNVERTCurrent.root","RECREATE");
	  number_particles->Write();
	  totxsec->Write();
	  quasielaxsec->Write();
	  spixsec->Write();
	  dpixsec->Write();
	  morepixsec->Write();
	  reaxsec->Write();
	  Nvertices->Write();
	  final_tot->Write();
	  final_tot_p->Write();
	  final_tot_n->Write();
	  final_tot_pim->Write();
	  final_tot_pi0->Write();
	  final_tot_pip->Write();

	  
	  final_tot_p->SetEntries(final_tot_p->GetEntries());
	  //quasielaxsec->Write();
	  //spixsec->Write();
	  //dpixsec->Write();
	  //morepixsec->Write();
	  //initmomxsec->Write();
	  //reaxsec->Write();

	  mmmm->Close();

	  TCanvas *c1 = new TCanvas("c1", "different scales hists",600,400);

	  TLegend* leg = new TLegend(0.15, 0.7, 0.55, 0.9);


	  totxsec->GetXaxis()->SetTitle("Momentum (MeV/c)");
	  totxsec->GetYaxis()->SetTitle("#sigma (mb)");
	  totxsec->SetAxisRange(0., 600.,"Y");
	  //	  totxsec->SetBinError(bin,(tot->GetBinError(bin))*targetArea/(all->GetBinContent(bin)));

	  totxsec->SetLineColor(kBlue);
	  totxsec->Draw();


	  quasielaxsec->SetLineColor(kRed);
	  quasielaxsec->Draw("SAME");
	  spixsec->SetLineColor(kCyan+1);
	  spixsec->Draw("SAME");    
	  dpixsec->SetLineColor(kGreen+1);
	  dpixsec->Draw("SAME");
	  morepixsec->SetLineColor(kBlack);
	  morepixsec->Draw("SAME");
	  //initmomxsec->SetLineColor(kMagenta);
	  //initmomxsec->Draw("SAME");
	  reaxsec->SetLineColor(kOrange+1);
	  reaxsec->Draw("SAME");
	  //"chist"
	  /*
	  double mom_totP[] = {525.2,658.5,830.9,965.5,
			       628.4,684.3,772.4,817.2,879.8,881.3,914.6,943.0,1011.0,1027.3,1031.3,1040.8,1046.1,1086.1,1092.7,1098.0,1113.7,1155.2,1159.1,1166.8,
			       1545.8};
	  double xsec_totP[] = {353.,296.,292.,285.,
				305.,290.,281.,284.,286.,286.,293.,295.,299.,303.,304.,299.,309.,312.,318.,314.,315.,324.,324.,322.,
				365.};
	  double err_totP[] = {7.,3.,6.,14.,
			       9.,12.,4.,6.,6.,3.,5.,4.,5.,7.,7.,8.,4.,8.,5.,4.,5.,11.,10.,6.,
			       3.};


	  TGraphErrors *ProCdata = new TGraphErrors(sizeof(mom_totP)/sizeof(double),mom_totP,xsec_totP,0,err_totP);
	  ProCdata->Draw("P");
	  ResetStyle(leg, 0, 1);
	  */
	  leg->AddEntry(totxsec,"total");
	  leg->AddEntry(quasielaxsec,"quasi-elastic");
	  leg->AddEntry(spixsec,"1 #pi production");
	  leg->AddEntry(dpixsec,"2 #pi production");
	  leg->AddEntry(morepixsec,"more #pi production");
	  //leg->AddEntry(initmomxsec,"initial momentum");
	  leg->AddEntry(reaxsec,"reactive");
	  leg->Draw("SAME");


	  c1->Update();
	  c1->SaveAs("FinalOutputNVERTCurrent.root");
	  //	  c1->Update();
	  
	  //	  TCanvas *c2 = new TCanvas("c2", "Momentum Output",600,400);
	  //c2->Update();
	  //c2->SaveAs("interactionChannel_testcase2.root");
	  //c2->Update();



  }

