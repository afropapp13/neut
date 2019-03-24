#include <TApplication.h>
#include <TRint.h>
#include <TROOT.h> 
#include <TCanvas.h> 
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TFile.h>
#include "N1p1h.h"
#include "timer.h"
#include <math.h>

double Random(void);

int main(int argc, char **argv) {

 char inputtensor[256]; 

 sprintf(inputtensor,"C_FULL_old.had"); 

 double Enu = 1.; 
 double tmu = 0.69; 
 double xcos = 0.9; 
 int nt; 
 int inu = 1;
 double tpmin; 
 //double xCrossSect; 

 // int idnucl[4];
 
 double pm,cos;
 double crosx; 
 double r,dp,tp;
 int Nuc; 
 double Pnu[4],Pt[4],Pl[4],Pp[4]; 
 // TRint *theApp = new TRint("ROOT example", &argc, argv);

 TFile *fluxFile = new TFile("../neutsmpl/skflux_dump_11a_511.root");

 TH1D *flux;
 
 fluxFile->GetObject("t2k_skflux_numu",flux); 
 
 if( argc < 3 ) exit(0); 

 bool  bind = false; 
  
 N1p1h n1p1h(argv[1],12,bind); 

 n1p1h.InitializeNucleus(16,true); 
 
 char name[256];

 sprintf(name,"%s",argv[2]); 

 TFile hfile(name,"RECREATE","1p1h Nieves test tree"); 
 
 TTree *tree = new TTree("N1p1h","n1p1h");
 
 tree->Branch("enu",&Enu,"enu/D"); 
 tree->Branch("pm",&pm,"pm/D");
 tree->Branch("cos",&xcos,"cos/D");
 tree->Branch("r",&r,"r/D");
 tree->Branch("dp",&dp,"dp/D");
 tree->Branch("tp",&tp,"tp/D");
 tree->Branch("crosx",&crosx,"crosx/D");
 tree->Branch("ntries",&nt,"ntries/I"); 
 tree->Branch("tpmin",&tpmin,"tpmin/D"); 
 
 tree->Branch("Pnu",Pnu,"Pnu[4]/D");
 tree->Branch("Pl",Pl,"Pl[4]/D");
 tree->Branch("Pt",Pt,"Pt[4]/D");
 tree->Branch("Pp",Pp,"Pp[4]/D");
 
 tree->Branch("Nuc",&Nuc,"Nuc/I");
 
 // ht.SetRFG(); 

 double Emax = n1p1h.GetEmax(); 

 TH1D electron("electron","  ",1000,0.,Emax); 
 TH1D positron("positron","  ",1000,0.,Emax); 
 TH1D muon("muon","  ",1000,0.,Emax); 
 TH1D amuon("amuon","  ",1000,0.,Emax); 
 TH1D tau("tau","  ",1000,0.,Emax); 
 TH1D atau("atau","  ",1000,0.,Emax); 
 TH2D tmucos("tmucos"," ",100,0.,3.,100,-1.,1.); 
 TH1F tries("tries"," ",1000,0,10e+5);
 TH1F RejEnu("RejEnu"," ",1000,0,Emax);
 
 TH1D electronO("electronO","  ",1000,0.,Emax); 
 TH1D positronO("positronO","  ",1000,0.,Emax); 
 TH1D muonO("muonO","  ",1000,0.,Emax); 
 TH1D amuonO("amuonO","  ",1000,0.,Emax); 

 TH1D electronmax("maxelectron","  ",1000,0.,Emax); 
 TH1D positronmax("maxpositron","  ",1000,0.,Emax); 
 TH1D muonmax("maxmuon","  ",1000,0.,Emax); 
 TH1D amuonmax("maxamuon","  ",1000,0.,Emax); 
 TH1D taumax("maxtau","  ",1000,0.,Emax); 
 TH1D ataumax("maxatau","  ",1000,0.,Emax);

 TH1D muonCO("muonCO","  ",1000,0.,Emax); 
 TH1D amuonCO("amuonCO","  ",1000,0.,Emax); 

 TH1D elCO("elCO","  ",1000,0.,Emax); 
 TH1D posCO("posCO","  ",1000,0.,Emax); 

 double pn,vc; 
 
 for( int i = 0. ; i < 1000 ; i++ ) {
   
   double E = electron.GetBinCenter(i+1);

   double xel = n1p1h.IntegralCrossSection(12,12,E);
   double xpos = n1p1h.IntegralCrossSection(-12,12,E); 
   
   double xmu = n1p1h.IntegralCrossSection(14,12,E);
   double xamu = n1p1h.IntegralCrossSection(-14,12,E); 

   double xelO = n1p1h.IntegralCrossSection(12,16,E);
   double xposO = n1p1h.IntegralCrossSection(-12,16,E); 
   
   double xmuO = n1p1h.IntegralCrossSection(14,16,E);
   double xamuO = n1p1h.IntegralCrossSection(-14,16,E); 

   double xtau = n1p1h.IntegralCrossSection(16,12,E);
   double xatau = n1p1h.IntegralCrossSection(-16,12,E); 

   electron.SetBinContent(i+1,xel); 
   positron.SetBinContent(i+1,xpos); 
   muon.SetBinContent(i+1,xmu); 
   amuon.SetBinContent(i+1,xamu); 
   tau.SetBinContent(i+1,xtau); 
   atau.SetBinContent(i+1,xatau); 

   electronO.SetBinContent(i+1,xelO); 
   positronO.SetBinContent(i+1,xposO); 
   muonO.SetBinContent(i+1,xmuO); 
   amuonO.SetBinContent(i+1,xamuO); 

   double rat =  (xmuO-xmu)/xmu;

   muonCO.SetBinContent(i+1,rat);

   rat =  (xamuO-xamu)/xamu;
   amuonCO.SetBinContent(i+1,rat);

   rat =  (xelO-xel)/xel;
   elCO.SetBinContent(i+1,rat);

   rat =  (xposO-xpos)/xpos;
   posCO.SetBinContent(i+1,rat);


   double mxel = n1p1h.CrossSectionmax(12,12,E);
   double mxpos = n1p1h.CrossSectionmax(-12,12,E); 
   
   double mxmu = n1p1h.CrossSectionmax(14,12,E);
   double mxamu = n1p1h.CrossSectionmax(-14,12,E); 

   double mxtau = n1p1h.CrossSectionmax(16,12,E);
   double mxatau = n1p1h.CrossSectionmax(-16,12,E); 

   electronmax.SetBinContent(i+1,mxel); 
   positronmax.SetBinContent(i+1,mxpos); 
   muonmax.SetBinContent(i+1,mxmu); 
   amuonmax.SetBinContent(i+1,mxamu); 
   taumax.SetBinContent(i+1,mxtau); 
   ataumax.SetBinContent(i+1,mxatau); 

 }

  //  ht.SetDebug(true); 
 
 // std::cout << xCrossSect <<  "  " << ht.DoubleDifferential(14,12,Enu,tmu,xcos) << std::endl;  

 electron.SetLineColor(2);
 positron.SetLineColor(2);
 positron.SetLineStyle(2); 
 muon.SetLineColor(3);
 amuon.SetLineColor(3);
 amuon.SetLineStyle(2); 
 tau.SetLineColor(4);
 atau.SetLineColor(4);
 atau.SetLineStyle(2); 

 electronmax.SetLineColor(2);
 positronmax.SetLineColor(2);
 positronmax.SetLineStyle(2); 
 muonmax.SetLineColor(3);
 amuonmax.SetLineColor(3);
 amuonmax.SetLineStyle(2); 
 taumax.SetLineColor(4);
 ataumax.SetLineColor(4);
 ataumax.SetLineStyle(2); 

 TCanvas *c1 = new TCanvas("c1"); 
 c1->cd(); 
 electron.Draw("L C"); 
 positron.Draw("L C same");
 muon.Draw("L C same"); 
 amuon.Draw("L  C same"); 
 tau.Draw("L C same"); 
 atau.Draw("L C same"); 
 c1->Update(); 

 TCanvas *c2 = new TCanvas("c1"); 
 c2->cd(); 
 electronmax.Draw("L C"); 
 positronmax.Draw("L C same");
 muonmax.Draw("L C same"); 
 amuonmax.Draw("L  C same"); 
 taumax.Draw("L C same"); 
 ataumax.Draw("L C same"); 
 c2->Update(); 
 
 Emax = n1p1h.GetMaxEnergy(); 
 
 std::cout << " ..................................... " << std::endl; 
 std::cout << " Get events " << std::endl; 
 std::cout << " ..................................... " << std::endl;  

 timer t; 
 int ntries = 0; 

 double tm = 0.; 

 int evtot = 10000000;

 TH1D *fluxO = (TH1D*)flux->Clone("fluxO");
 TH1D *fluxC = (TH1D*)flux->Clone("fluxC"); 
 
 for(int i = 0; i < flux->GetNbinsX(); i++ ) {
   double E = flux->GetBinCenter(i+1);
   double fval = flux->GetBinContent(i+1); 
   fluxO->SetBinContent(i+1,muonO.Interpolate(E)*fval);
   fluxC->SetBinContent(i+1,muon.Interpolate(E)*fval); 
 }

 
 for(int i = 0; i < evtot; i++ ) {
  double TLepton,coslepton; 
  //  Enu = (Emax-0.2)*Random()+0.2; 

  Nuc = 12; 
  if( Random() < 0.5 ) Nuc = 16; 

  if( Nuc == 12 ) 
    Enu = fluxC->GetRandom(); 

  if( Nuc == 16 )
    Enu = fluxO->GetRandom(); 
  
  if( Enu < 0.15 ) { i--; continue; }
  
  double pnu[4],p[4][4];
  int idpart[4],parent[4];
  double R;

  pnu[1] = pnu[2] = 0; 
  pnu[0] = pnu[3] = Enu; 
    
  t.restart(NULL);

  n1p1h.GenerateVectors(14,Nuc,pnu,p,idpart,parent,R,crosx); 

  crosx *= (double) Nuc/2.;
  
  tm += t.elapsed_time(); 

  ntries++;

  pm = sqrt(p[2][1]*p[2][1]+p[2][2]*p[2][2]+p[2][3]*p[2][3])/1000.;
  xcos= (p[0][1]*p[2][1]+p[0][2]*p[2][2]+p[0][3]*p[2][3])/(sqrt(p[0][1]*p[0][1]+p[0][2]*p[0][2]+p[0][3]*p[0][3])*pm*1000.);
  tp = sqrt(p[1][1]*p[1][1]+p[1][2]*p[1][2]+p[1][3]*p[1][3])/1000.; 
  dp=sqrt(p[3][1]*p[3][1]+p[3][2]*p[3][2]+p[3][3]*p[3][3])/1000.;
  r=R;

  for(int i=0;i<4;i++){
    Pnu[i]=p[0][i];
    Pl[i]=p[2][i];
    Pt[i]=p[1][i];
    Pp[i]=p[3][i];
  }

  if( ((ntries*100)%evtot) == 0 ) {
    std::cout << " Processed " << i*100/evtot << " % of the events " << std::endl; 
    std::cout << crosx << "  " << Enu << "  " << tmu << "  " << xcos << "   " << r << "  " << dp <<  "  " << vc << std::endl;
    std::cout << " time per event " << tm/(double)ntries << " sec" << std::endl;   
    std::cout << " Fraction of success " << (double)i/(double)ntries*100 << "%" << std::endl; 
  }

  if( crosx > 0 ) {
    nt = n1p1h.GetNumberofItLastEvent(); 
    tpmin = n1p1h.GetMinTargetMomentum(); 
    tree->Fill();
    tries.Fill(nt);   
  }
  else {
    RejEnu.Fill(Enu); 
    i--;
  } 
 }

 tree->Write();
 electron.Write();
 positron.Write(); 
 muon.Write();
 amuon.Write();
 electronO.Write();
 positronO.Write(); 
 muonO.Write();
 amuonO.Write();
 tau.Write();
 atau.Write();
 muonCO.Write();
 amuonCO.Write();
 elCO.Write();
 posCO.Write();

 electronmax.Write();
 positronmax.Write(); 
 muonmax.Write();
 amuonmax.Write();
 taumax.Write();
 ataumax.Write(); 

 tries.Write();
 RejEnu.Write();

 hfile.Close(); 

}
