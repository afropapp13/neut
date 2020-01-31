#include <TROOT.h> 
#include <TRint.h>
#include <TCanvas.h> 
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TFile.h>

#include "HT2p2h.h"
#include <math.h>

#if 1 
extern"C" {

  // CC2p2h data blocks

  // void readhadtensor_(char *name);

  //  void diffcrosssection_(int *inu,double *Enu,double *Tmu,double *xcos,double *xCrossSect);

  struct nieves2p2hpar_common nieves2p2hpar_;
  
  float rlu_(int *dummy){
	
#define RNDMX  ((double) 0x7fffffff)
    double a = random()/RNDMX;
	
    return a; 
	
  }
  
}
#endif 

int main(int argc, char **argv) {

 char inputtensor[256]; 

 sprintf(inputtensor,"C_FULL_old.had"); 

 double Enu = 1.; 
 double tmu = 0.69; 
 double xcos = 0.9; 
 int inu = 1;
 double xCrossSect; 

 int idnucl[4];
 
 double pm[4],q[4],qi1[4],qi2[4],qf1[4],qf2[4];

 double r;
 
 double val = 0.0 ;
 
 bool applybind = true;

 // nieves2p2hpar_.nv2p2hqval = 1; /* qval is set to 0 */
 nieves2p2hpar_.nv2p2hqval = 2; /* r-dependent qval */

 int hmode = 0; 
 
 if( argc > 1 ) {

   for(int ii = 0; ii < argc; ii++ ) {

     if( !strncmp(argv[ii],"-b2b",4) ) {
       std::cout << " Back to back kinematics " << std::endl;
       hmode = 1;
     } else if( !strncmp(argv[ii],"-f2f",4) ) {
       std::cout << " Forward to forward kinematics " << std::endl;
       hmode = 2; 
     } else if ( !strncmp(argv[ii],"-nbind",6) ) {
       applybind = false; 
     } else if ( !strncmp(argv[ii],"-hap0",5) ) {
       ii++;
       val = atof(argv[ii]);  
     }
   }
 }

 argc = 1; 

 TRint *theApp = new TRint("ROOT example", &argc, argv);
 
 HT2p2h ht("./HT2p2h/Nieves",applybind,false); 

 ht.SetHFSparameters(val,0.);
 ht.SetHadronKinMode(hmode);
 
 char name[256];

 sprintf(name,"HT2p2h_V%f_M%d.root",val,hmode); 
 
 if( argc > 2 ) 
   if( atoi(argv[2]) )
       sprintf(name,"HT2p2h_%f_NoBind.root",val); 

 
 TFile hfile(name,"RECREATE","2p2h HT test tree"); 
 
 TTree *tree = new TTree("HT2p2h","2p2h");
 
 tree->Branch("enu",&Enu,"enu/D"); 
 tree->Branch("pm",pm,"em/D:pxm/D:pym/D:pzm/D");
 tree->Branch("r",&r,"r/D");
 tree->Branch("pi1",qi1,"ei1/D:pxi1/D:pyi1/D:pzi1/D");
 tree->Branch("pi2",qi2,"ei2/D:pxi2/D:pyi2/D:pzi2/D");
 tree->Branch("pf1",qf1,"ef1/D:pxf1/D:pyf1/D:pzf1/D");
 tree->Branch("pf2",qf2,"ef2/D:pxf2/D:pyf2/D:pzf2/D");
 tree->Branch("idn",idnucl,"idi1/I:idi2/I:idf1/I:idf2/I");
 
 // ht.SetRFG(); 
 
 TH1D electron("electron","  ",1000,0.,30.); 
 TH1D positron("positron","  ",1000,0.,30.); 
 TH1D muon("muon","  ",1000,0.,30.); 
 TH1D amuon("amuon","  ",1000,0.,30.); 
 TH1D tau("tau","  ",1000,0.,30.); 
 TH1D atau("atau","  ",1000,0.,30.); 
 TH2D tmucos("tmucos"," ",100,0.,3.,100,-1.,1.); 
 TH1F tries("tries"," ",100,-0.5,99.5); 

 TH2D initial("initial","   ",100,0.,.300,100,0.,.300);
 TH2D final("final","   ",100,0.,2.000,100,0.,2.00); 

 TH1D htries("htries"," ",10007,-6.5,10000.5); 

 TH1D EBalance("EBalance","",1000,-0.0001,0.0001);
 TH1D pxBalance("pxBalance","",1000,-0.0001,0.0001);
 TH1D pyBalance("pyBalance","",1000,-0.0001,0.0001);
 TH1D pzBalance("pzBalance","",1000,-0.0001,0.0001);
 TH1D RH("RH","",1000,0.,10.); 


 
 // std::cout << " Cross Section test " << ht.DoubleDifferential(14,12,1.00164 , 0.584074 , 0.652624 ) << std::endl;
 
 for( int i = 0. ; i < 1000 ; i++ ) {
   
   double E = electron.GetBinCenter(i+1);

   double xel = ht.IntegralCrossSection(12,12,E);
   double xpos = ht.IntegralCrossSection(-12,12,E); 
   
   double xmu = ht.IntegralCrossSection(14,12,E);
   double xamu = ht.IntegralCrossSection(-14,12,E); 

   double xtau = ht.IntegralCrossSection(16,12,E);
   double xatau = ht.IntegralCrossSection(-16,12,E); 
   
   electron.SetBinContent(i+1,xel); 
   positron.SetBinContent(i+1,xpos); 
   muon.SetBinContent(i+1,xmu); 
   amuon.SetBinContent(i+1,xamu); 
   tau.SetBinContent(i+1,xtau); 
   atau.SetBinContent(i+1,xatau); 
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

 TCanvas *c1 = new TCanvas("c1"); 
 c1->cd(); 
 electron.Draw("L C"); 
 positron.Draw("L C same");
 muon.Draw("L C same"); 
 amuon.Draw("L  C same"); 
 tau.Draw("L C same"); 
 atau.Draw("L C same"); 
 c1->Update(); 

 for(int i = 0; i < 100000; i++ ) {
   double TLepton,coslepton; 
   Enu = 2.8*ht.Random()+0.2; 


   double pnu[4];
   double p[6][4];
   int idpart[6];
   int parent[6];
  

   pnu[0]=pnu[3]=Enu;pnu[1]=pnu[2]=0.;

   int hadrontries = ht.GenerateVectors(14,12,pnu,p,idpart,parent,r);

   if( hadrontries < 0 ) { i--; continue;}
   
   RH.Fill(r);
   
   for(int i = 0; i < 4; i++ ) {
     pm[i] = p[3][i];
     qi1[i] = p[1][i]; qi2[i] = p[2][i];qf1[i] = p[4][i]; qf2[i] = p[5][i];
     q[i] = pnu[i]-pm[i];
   }

   tree->Fill(); 

   if( q[0]+qi1[0]+qi2[0]-qf1[0]-qf2[0] > 0.04 ) std::cout << q[0] << " " << q[1] << "  " << q[2] << "  " << q[3] << std::endl;
   
   EBalance.Fill(q[0]+qi1[0]+qi2[0]-qf1[0]-qf2[0]);
   pxBalance.Fill(q[1]+qi1[1]+qi2[1]-qf1[1]-qf2[1]);
   pyBalance.Fill(q[2]+qi1[2]+qi2[2]-qf1[2]-qf2[2]);
   pzBalance.Fill(q[3]+qi1[3]+qi2[3]-qf1[3]-qf2[3]); 

   
 }
 
 tree->Write();
 electron.Write();
 positron.Write(); 
 muon.Write();
 amuon.Write();
 tau.Write();
 atau.Write(); 

 hfile.Close(); 

 
 
 TCanvas *c2 = new TCanvas("c2"); 
 c2->cd(); 
 tmucos.Draw("colz"); 
 c2->Update();


 TCanvas *c3 = new TCanvas("c3"); 
 c3->cd(); 
 tries.Draw(); 
 c3->Update();


TCanvas *c4 = new TCanvas("c4");
 c4->Divide(2,2); 
 c4->cd(1); 
 initial.Draw("colz");
 c4->cd(2); 
 final.Draw("colz");
 c4->cd(3);
 htries.Draw();
 c4->cd(4);
 RH.Draw(); 
 c4->Update();

 TCanvas *c5 = new TCanvas("c5");
 c5->Divide(2,2); 
 c5->cd(1); 
 EBalance.Draw();
 c5->cd(2); 
 pxBalance.Draw();
 c5->cd(3);
 pyBalance.Draw();
 c5->cd(4);
 pzBalance.Draw();
 c5->Update();
 
 theApp->Run(); 
 
}



