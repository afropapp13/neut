//Konosuke Iwamoto
//June 8, 2014

//DESCRIPTION:

//Script to read in muon & neutrino kinematics in text file, attach photon under certain probability

// 3/31/2014 17:43 EST - Completed radiative correction script using interpolation and integrals of 2D distributions and splines that are already made; no need to run over text file multiple times.
// 6/8/2014 17:43 EST - Completed radiative correction script using interpolation for probability density of photon as well.
// 7/8/2014  14:58 EST - Turned off FV constraint
// 7/30/2014 15:37 EST - nue rad CCQE added
// 10/29/2014 16:56 JST - Modified for NEUT
// 11/11/2014 16:32 JST - Completed the modification for NEUT

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>

#include <TROOT.h>
#include <TString.h>
#include <TMath.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TSpline.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TRandom.h>
#include <TF1.h>
#include <TArrow.h>
#include "necardC.h"
#include "vcworkC.h"
#include "neworkC.h"
#include "neutfilepathC.h"
#include "posinnucC.h"

class SplineEval
{
public:
    SplineEval(TSpline3* s)
    {
        spline=s;
    }
    
    Double_t operator()(Double_t* x, Double_t* /*p*/)
    {
        if ((spline->Eval(x[0])) >0)
        {
            return (spline->Eval(x[0]));
        }
        else
        {
            return 0;
        }
    }
    
private:
    TSpline3* spline;
};

Double_t ProbCalCode(Double_t EP, Double_t CosTheta, Double_t Enu1, Double_t Enu2, Double_t ProbCal, TH2D *Hist, TH2D *HistA);
TF1* plotSP(TGraph *GRF, TSpline3 *Spl,Double_t *EGammaArray, Double_t *pArray);

extern "C"
{
  void radcorr();
  void radcorr_();
}

void radcorr_() // Fortran wrapper
{
  switch ( neutradcorr_.iradcorr ){
  case 1:
	radcorr();
	break;
  case 0:
	break;
	
  default:
	std::cerr << "Undefined iradcorr ( neut card ) value was set"
			  << neutradcorr_.iradcorr
			  << std::endl;
	exit(1);
  }
  return;
}

void radcorr() //Apply radiative corrections on muon samples
{
    
    //RESETS RANDOMIZATION
  //    gRandom->SetSeed();
    //////////////////////
  
////////////////////////////////////
    
///////////////VARIABLES TO READ INPUT FILE/////////////////////

    
    //Store necessary variables to calculate probability and photon energy
    
   Double_t pabsl = sqrt(vcwork_.pvc[2][0]*vcwork_.pvc[2][0]+vcwork_.pvc[2][1]*vcwork_.pvc[2][1]+vcwork_.pvc[2][2]*vcwork_.pvc[2][2]);
   Double_t pabsn = sqrt(vcwork_.pvc[0][0]*vcwork_.pvc[0][0]+vcwork_.pvc[0][1]*vcwork_.pvc[0][1]+vcwork_.pvc[0][2]*vcwork_.pvc[0][2]);
    
   Double_t dirl[3];
   Double_t dirn[3];
    for (int i =0;i<3;i++)
    {
   dirl[i] = vcwork_.pvc[2][i]/pabsl;
   dirn[i] = vcwork_.pvc[0][i]/pabsn;
    }
    
   Double_t CosTheta = (dirl[0]*dirn[0])+(dirl[1]*dirn[1])+(dirl[2]*dirn[2]);
   Double_t Ep1 = (pabsn)/1000.;//Converts MeV to GeV
    
    
    // variables used to store the calculated values
    TF1 *exmple;
    TGraph *GEXMSet=0;
    TSpline3 *SpEXMSet=0;
    
    Double_t exmplArray[24];
    
    Double_t ProbCal = 0.;
    Double_t ProbRand = 0.;
    Double_t ProbMax = 0.;
    
    
    Double_t Ephoton = 0.;
    
    
    Double_t Pil = 0.;
    Double_t Pfl = 0.;
    Double_t pf[3] = {0.,0.,0.};
    Double_t pp[3] = {0.,0.,0.};



    //Maximum and minimum of Integral
    const Double_t min = 0.03;
    const Double_t max = 3.5;
    
        const Double_t mm = 105.658366/1000.; //muon mass in GeV/c2
        const Double_t me = 0.510998/1000.;   //electron mass in GeV/c2
///////////////////////////////////////////////////////////////
    
    
    
    
//Opens necessary distributions and splines from root files for the first run
    static int flag;

    //MUON CCQE
   static TH2D *Hist[7];
   static TH2D *HistA[7];
   static TH2D *HistInt[6];
   static TH2D *HistIntA[6];
    
   static TH2D *HistEP[7][24];
   static TH2D *HistEPA[7][24];
    
   static TH2D *HistEPInt[6][24];
   static TH2D *HistEPIntA[6][24];
    
    ///
    
    //ELECTRON CCQE
   static TH2D *HistE[7];
   static TH2D *HistEA[7];
   static TH2D *HistEInt[6];
   static TH2D *HistEIntA[6];
    
   static TH2D *HistEEP[7][24];
   static TH2D *HistEEPA[7][24];
    
   static TH2D *HistEEPInt[6][24];
   static TH2D *HistEEPIntA[6][24];
    




    
    Double_t EnuList[37]= {0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.15, 1.35, 1.55, 1.75, 1.95, 2.1, 2.5, 2.9, 3.3, 3.7, 4.1, 4.75, 5.5, 6.25, 7., 7.75, 8.5, 9.25, 10.};
    
    Double_t CtList[28]={-0.95, -0.85, -0.75, -0.65, -0.55, -0.45, -0.35, -0.25, -0.15, -0.05, 0.05, 0.15, 0.25,0.35, 0.45, 0.55, 0.65, 0.75, 0.81, 0.83, 0.85, 0.87, 0.89, 0.91, 0.93, 0.95, 0.97, 0.99};
    
    
    Double_t EGList[24]={0.015,0.025,0.0375,0.0525,0.07,0.09,0.1125,0.1375,0.1625,0.1875,0.225,0.275,0.325,0.375,0.45,0.55,0.7,0.9,1.15,1.45,1.8,2.25,2.75,3.5};
    
    Double_t EPVal;
    
    Int_t INTEnuA[5] = {0,13,18,24,30};
    Int_t INTEnuB[5] = {12,17,22,28,36};
    Int_t INTHist[5] = {0,1,2,4,6};
    
    Int_t INTEnuIntA[6] = {12,17,22,23,28,29};
    Int_t INTEnuIntB[6] = {13,18,23,24,29,30};
    Int_t INTHistInt[6] = {0,1,2,3,4,5};
    
    TString numbering[7]={"0","1","2","3","4","5","6"};
    TString NMB[50]={"_1","_2","_3","_4","_5","_6","_7","_8","_9","_10","_11","_12","_13","_14","_15","_16","_17","_18","_19","_20","_21","_22","_23","_24","_25","_26","_27","_28","_29","_30","_31","_32","_33","_34","_35","_36","_37","_38","_39","_40","_41","_42","_43","_44","_45","_46","_47","_48","_49","_50"};
 





    ///
    
    if (flag == 0) //READ MODE (reads in the root files)
{
    TFile *MCfile = new TFile("RadnumuCCQE.root");
    TFile *MCfileE= new TFile("RadnueCCQE.root");
    TString mode[2] = {"M","E"};
    
   
    for (int b=0;b<7;b++)
    {
        Hist[b]= new TH2D(*(TH2D*)MCfile->Get("Prob"+mode[0]+numbering[b]));
        HistA[b]= new TH2D(*(TH2D*)MCfile->Get("ProbA"+mode[0]+numbering[b]));
        
        HistE[b]= new TH2D(*(TH2D*)MCfileE->Get("Prob"+mode[1]+numbering[b]));
        HistEA[b]= new TH2D(*(TH2D*)MCfileE->Get("ProbA"+mode[1]+numbering[b]));
        
    }
    
    
    for (int c=0;c<6;c++)
    {
        HistInt[c]= new TH2D(*(TH2D*)MCfile->Get("ProbInt"+mode[0]+numbering[c]));
        HistIntA[c]= new TH2D(*(TH2D*)MCfile->Get("ProbIntA"+mode[0]+numbering[c]));
        
        HistEInt[c]= new TH2D(*(TH2D*)MCfileE->Get("ProbInt"+mode[1]+numbering[c]));
        HistEIntA[c]= new TH2D(*(TH2D*)MCfileE->Get("ProbIntA"+mode[1]+numbering[c]));
        
    }
    
    
    
    
    for (int x = 0;x<7;x++)
    {
        for (int y =0;y<24;y++)
        {
            HistEP[x][y] = new TH2D(*(TH2D*)MCfile->Get("EP"+mode[0]+numbering[x]+NMB[y]));
            HistEPA[x][y] = new TH2D(*(TH2D*)MCfile->Get("EPA"+mode[0]+numbering[x]+NMB[y]));
            
            HistEEP[x][y] = new TH2D(*(TH2D*)MCfileE->Get("EP"+mode[1]+numbering[x]+NMB[y]));
            HistEEPA[x][y] = new TH2D(*(TH2D*)MCfileE->Get("EPA"+mode[1]+numbering[x]+NMB[y]));
            
            
            if (x<6)
            {
                HistEPInt[x][y] = new TH2D(*(TH2D*)MCfile->Get("EPInt"+mode[0]+numbering[x]+NMB[y]));
                HistEPIntA[x][y] = new TH2D(*(TH2D*)MCfile->Get("EPIntA"+mode[0]+numbering[x]+NMB[y]));
                HistEEPInt[x][y] = new TH2D(*(TH2D*)MCfileE->Get("EPInt"+mode[1]+numbering[x]+NMB[y]));
                HistEEPIntA[x][y] = new TH2D(*(TH2D*)MCfileE->Get("EPIntA"+mode[1]+numbering[x]+NMB[y]));
            }
            
            
        }
    }
    

    
    
    flag = 1;  //TURNS OFF THE READ MODE    
    
}


///////////////////////////////////////////////////////////////
    
    ///SET THE MAXIMUM PROBABILITY

    ProbMax = HistA[6]->Interpolate(10,0.99);

       
//////////////
    
    
    //Gets probability of generating photons using the distributions in ProbDist.root
    

    
    if (vcwork_.ipvc[0] == 16) return;
    if (vcwork_.icrnvc[2] == 0) return;
    
        
        if (vcwork_.ipvc[0] == 14)  //NUMU
        {
            
            if (Ep1 < EnuList[0])
            {
                ProbCal = 0;
            }
            
            for (int A=0;A<5;A++)
            {
                ProbCal = ProbCalCode(Ep1,CosTheta,EnuList[INTEnuA[A]],EnuList[INTEnuB[A]],ProbCal,Hist[INTHist[A]],HistA[INTHist[A]]);
            }
            
            for (int B=0;B<6;B++)
                
            {
                ProbCal = ProbCalCode(Ep1,CosTheta,EnuList[INTEnuIntA[B]],EnuList[INTEnuIntB[B]],ProbCal,HistInt[INTHistInt[B]],HistIntA[INTHistInt[B]]);
            }
            
            for (int C=0;C<24;C++)
            {
                for (int A=0;A<5;A++)
                {
                    EPVal = ProbCalCode(Ep1,CosTheta,EnuList[INTEnuA[A]],EnuList[INTEnuB[A]],ProbCal,HistEP[INTHist[A]][C],HistEPA[INTHist[A]][C]);
                }
                
                for (int B=0;B<6;B++)
                    
                {
                    EPVal = ProbCalCode(Ep1,CosTheta,EnuList[INTEnuIntA[B]],EnuList[INTEnuIntB[B]],ProbCal,HistEPInt[INTHistInt[B]][C],HistEPIntA[INTHistInt[B]][C]);
                }
                
                exmplArray[C] = EPVal;
                
            }
            
            
            
            
        }
        
        
        if (vcwork_.ipvc[0] == 12) //NUE
        {
            
            
            if (Ep1 < EnuList[0])
            {
                ProbCal = 0;
            }
            
            for (int A=0;A<5;A++)
            {
                ProbCal = ProbCalCode(Ep1,CosTheta,EnuList[INTEnuA[A]],EnuList[INTEnuB[A]],ProbCal,HistE[INTHist[A]],HistEA[INTHist[A]]);
            }
            
            for (int B=0;B<6;B++)
                
            {
                ProbCal = ProbCalCode(Ep1,CosTheta,EnuList[INTEnuIntA[B]],EnuList[INTEnuIntB[B]],ProbCal,HistEInt[INTHistInt[B]],HistEIntA[INTHistInt[B]]);
            }
            
            
            
            
            
            for (int C=0;C<24;C++)
            {
                for (int A=0;A<5;A++)
                {
                    EPVal = ProbCalCode(Ep1,CosTheta,EnuList[INTEnuA[A]],EnuList[INTEnuB[A]],ProbCal,HistEEP[INTHist[A]][C],HistEEPA[INTHist[A]][C]);
                }
                
                for (int B=0;B<6;B++)
                    
                {
                    EPVal = ProbCalCode(Ep1,CosTheta,EnuList[INTEnuIntA[B]],EnuList[INTEnuIntB[B]],ProbCal,HistEEPInt[INTHistInt[B]][C],HistEEPIntA[INTHistInt[B]][C]);
                }
                
                exmplArray[C] = EPVal;
                
            }
            
          
            
           
         


        }

	exmple = plotSP(GEXMSet, SpEXMSet, EGList, exmplArray);
	Pil = pabsl/1000.;

  if (vcwork_.ipvc[0] == 14)  //NUMU
        {
	  Ephoton = exmple->GetRandom(min,sqrt(Pil*Pil+mm*mm)-mm);
	  Pfl = sqrt((sqrt(Pil*Pil+mm*mm)-Ephoton)*(sqrt(Pil*Pil+mm*mm)-Ephoton)-mm*mm) ;
	}

  else if (vcwork_.ipvc[0] == 12) //NUE
        {
          Ephoton = exmple->GetRandom(min,sqrt(Pil*Pil+me*me)-me);
	  Pfl = sqrt((sqrt(Pil*Pil+me*me)-Ephoton)*(sqrt(Pil*Pil+me*me)-Ephoton)-me*me);
	}
	    

	   for (int i=0;i<3;i++)
            {
	      pf[i] = 1000.*Pfl*dirl[i];
	      pp[i] = 1000.*Ephoton*dirl[i];
            }
                 
        
        //Generates random value uniform in range [0,1], which is used for the criteria to attach photon or not
        ProbRand = gRandom->Uniform(0,1);
       

       ////////////
        
        
	     if (ProbCal>ProbRand) //ATTACH PHOTON     
       {
            
            //ADDING PHOTON INFORMATION

            vcwork_.nvc++;
            vcwork_.ipvc[vcwork_.nvc-1] = 22;
            vcwork_.amasvc[vcwork_.nvc-1] = 0.;
            vcwork_.iorgvc[vcwork_.nvc-1] = 3;
            vcwork_.iflgvc[vcwork_.nvc-1] = 8;
            vcwork_.icrnvc[vcwork_.nvc-1] = 1;
            vcwork_.timvc[vcwork_.nvc-1]  = vcwork_.timvc[2];
            vcwork_.ivtivc[vcwork_.nvc-1]   = vcwork_.ivtivc[2];
            vcwork_.ivtfvc[vcwork_.nvc-1]   = vcwork_.ivtfvc[2];

            for (int i=0;i<3;i++)
            {
                vcwork_.pvc[vcwork_.nvc-1][i] = pp[i];
                vcwork_.posivc[vcwork_.nvc-1][i]=  vcwork_.posivc[2][i];
                vcwork_.posfvc[vcwork_.nvc-1][i]=  vcwork_.posfvc[2][i];
		posinnuc_.posnuc[vcwork_.nvc-1][i]=posinnuc_.posnuc[2][i];
                vcwork_.pvc[2][i] = pf[i];
            }


       }


}


Double_t ProbCalCode(Double_t EP, Double_t CosTheta, Double_t Enu1, Double_t Enu2, Double_t ProbCal, TH2D *Hist, TH2D *HistA)
//Gets probability of generating photons using the distributions in ProbDist.root
{
    
    if(EP >= Enu1 && EP < Enu2)
    {
        if(CosTheta < 0.8)
        {
            ProbCal = Hist->Interpolate(EP,CosTheta);
            
        }
        if(CosTheta > 0.8)
        {
            ProbCal = HistA->Interpolate(EP,CosTheta);
        }
    }
    
    return ProbCal;
    
}

TF1* plotSP(TGraph *GRF, TSpline3 *Spl,Double_t *EGammaArray, Double_t *pArray)

{
    GRF=new TGraph(23,EGammaArray,pArray);
    Spl= new TSpline3("",GRF);
    SplineEval* se=new SplineEval(Spl);
    TF1* f1Spline=new TF1("f1", se, 0, 3.5, 0, "SplineEval");
    return f1Spline;
}

