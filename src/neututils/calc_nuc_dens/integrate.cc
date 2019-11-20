//////////////////////////////////////////////////////
//
// Purpose: Make table of integrals of nuclear density
//
//////////////////////////////////////////////////////

#include <iostream>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>

#include "efpionC.h"
#include "necardC.h"

#include "neutnuclei.h"

using namespace std;

extern "C" {
float efabrho_(double *r);
void nesettarg_();
}

double efabrho(double *x, double *par) {

  double r = x[0];
  double f = efabrho_(&r);
  return f;
}

double efabrho_rsq(double *x, double *par) {

  double r = x[0];
  double f =
      r * r * efabrho_(&r); // *4*pi (cancels when we normalize in analyze.cc)
  return f;
}

int main(int argc, char *argv[]) {

  TFile *outputfile = new TFile("nucldens.root", "RECREATE");

  for (int iNuc = 0; iNuc < nNuclei; iNuc++) {

    cout << "Proton # = " << nBoundProtons[iNuc] << endl;
    ;

    neuttarget_.numbndp = nBoundProtons[iNuc];
    nesettarg_();

    TF1 fn_efabrho(Form("efabrho_%dprotons", nBoundProtons[iNuc]), efabrho, 0,
                   max_r, 0);
    TF1 fn_efabrho_rsq("efabrho_rsq", efabrho_rsq, 0, max_r, 0);
    fn_efabrho.Write();

    TGraph g_efabrho_int;
    g_efabrho_int.SetName(Form("efabrho_%dprotons_int", nBoundProtons[iNuc]));

    int nPoints = 0;
    double r_last = 0;
    double integral = 0;
    for (double r = 0; r <= 2.5 * eftarget_.c; r += r_step) {

      integral += fn_efabrho_rsq.Integral(r_last, r);
      r_last = r;

      g_efabrho_int.SetPoint(nPoints, r, integral);
      nPoints++;

    }

    g_efabrho_int.Write();
  }
}
