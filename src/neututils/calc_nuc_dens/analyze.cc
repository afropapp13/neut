#include "neutnuclei.h"

void analyze() {
  
  TFile *infile = new TFile("nucldens.root");
  TFile *outfile = new TFile("nucldens_norm.root","RECREATE");

  TGraph *g_efabrho_int[nNuclei];
  
  // Load graphs
  for (int iNuc=0; iNuc<nNuclei; iNuc++) {
    g_efabrho_int[iNuc] = new TGraph(*(TGraph*)infile->Get(Form("efabrho_%dprotons_int",nBoundProtons[iNuc])));
  }

  // Choose scaling factor

  // Final Coulomb factor for use in eftrace.F (see coulomb_factor.xls)
  Double_t scale_factor = 1.439964;

  // Charge (e)
  //Double_t scale_factor = 1;

  // Normalize
  for (int iNuc=0; iNuc<nNuclei; iNuc++) {
    
    int nPoints = g_efabrho_int[iNuc]->GetN();
    
    double x,max;
    g_efabrho_int[iNuc]->GetPoint(nPoints-1,x,max);
    
    for (int ipoint=0; ipoint<nPoints; ipoint++) {
      double y;
      
      g_efabrho_int[iNuc]->GetPoint(ipoint,x,y);

      // Normalized
      //g_efabrho_int[iNuc]->SetPoint(ipoint,x,y/max);
      
      g_efabrho_int[iNuc]->SetPoint(ipoint,x,y/max*nBoundProtons[iNuc]*scale_factor);

    }
    
    g_efabrho_int[iNuc]->Write();
  }
  

  // Write out to text file
  FILE * pFile;
  pFile = fopen ("nucdens_int_coulomb.dat","w");
  fprintf (pFile, "     ");

  for (int iNuc=0; iNuc<nNuclei; iNuc++) {
    fprintf (pFile, "%11.d ",nBoundProtons[iNuc]);
  }
  fprintf (pFile, "\n");

  int i=0;
  for (float r=0; r<=max_r; r+=r_step) {

    fprintf (pFile, "%4.1f ",r);

    for (int iNuc=0; iNuc<nNuclei; iNuc++) {

      Double_t x,y;
      if (g_efabrho_int[iNuc]->GetPoint(i,x,y)<0)
	y = nBoundProtons[iNuc]*scale_factor;

      fprintf (pFile, "%11e ",y);

    }

    fprintf (pFile, "\n");

    i++;
  }
  /**/
  fclose (pFile);	

  /*
  // Draw
  TCanvas *c = new TCanvas(1);
  for (int iNuc=nNuclei-1; iNuc>=0; iNuc--) {
    g_efabrho_int[iNuc]->SetLineColor(iNuc+1);
    
    if (iNuc==nNuclei-1)
      g_efabrho_int[iNuc]->Draw("AC");
    else
      g_efabrho_int[iNuc]->Draw("C same");
  }
  /**/
}
