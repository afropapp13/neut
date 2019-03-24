#include <TApplication.h>
#include <TRint.h>
#include <TROOT.h> 
#include <TCanvas.h> 
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TFile.h>
#include <TMinuit.h> 
#include "N1p1h.h"
#include <math.h>
#include "f2cdeclarations.h"

double Random(void);

double Enu = 1.; 

void xsfunction(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

  int NRC = 10; int NPC = 10;

  double TLep = par[0];
  double xcos = par[1];
  double R = par[2];
  double dP = par[3];
  double vc; 
  
  double xsec = sigthbruno_(&Enu,&TLep,&xcos,&NRC,&NPC,&R,&dP,&vc);

  std::cout << xsec << std::endl; 
  
  if( xsec <= 0 )
    f = 1e+33; 
  else
    f = log(xsec);

  return; 
}



int main(int argc, char **argv) {
 int inu = 1;
 double xCrossSect; 

 int idnucl[4];
 
 // TRint *theApp = new TRint("ROOT example", &argc, argv);

 double val = 0.0 ;
 
 if( argc > 1 )
   val = atof(argv[1]); 

 bool applybind = true; 
  
 N1p1h n1p1h(applybind); 

 TMinuit *ptMinuit = new TMinuit(4);  //initialize TMinuit with a maximum of 5 params
  //
  //  select verbose level:
  //    default :     (58 lines in this test)
  //    -1 : minimum  ( 4 lines in this test)
  //     0 : low      (31 lines)
  //     1 : medium   (61 lines)
  //     2 : high     (89 lines)
  //     3 : maximum (199 lines in this test)
  //
  ptMinuit->SetPrintLevel();
  // set the user function that calculates chi_square (the value to minimize)
  ptMinuit->SetFCN(xsfunction);

  Double_t arglist[10];
  Int_t ierflg = 0;

  arglist[0] = 1;
  ptMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  Enu = 1.;
  double Tlep = 1.-0.1056;
  double xcos = 1.;
  double R = 3.;
  double dp = 0.250;
  
  // Set starting values and step sizes for parameters
  static Double_t vstart[4] = {Tlep, xcos, R , dp};
  static Double_t step[4] = {Tlep*1. , xcos*.1 ,R*.1 ,dp*.1};
  ptMinuit->mnparm(0, "tl", vstart[0], step[0], 0,0,ierflg);
  ptMinuit->mnparm(1, "cos", vstart[1], step[1], 0,0,ierflg);
  ptMinuit->mnparm(2, "R", vstart[2], step[2], 0,0,ierflg);
  ptMinuit->mnparm(3, "fp", vstart[3], step[3], 0,0,ierflg);

  ptMinuit->mnexcm("SCAN", arglist ,0,ierflg);
  
  // Now ready for minimization step
  arglist[0] = 500;
  arglist[1] = 1.;
  ptMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

}
