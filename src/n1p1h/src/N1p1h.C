#include "N1p1h.h"
#include <dirent.h>
#include <stdlib.h>
#include <math.h>
#include <TApplication.h>
#include <TRint.h>
#include <TROOT.h> 
#include <TFile.h>
#include <TTree.h>
#include "f2cdeclarations.h"
#include <TH1F.h>
#include <fstream>

#ifdef WITH_NEUT
#include "neutmodelC.h"
#include "nieves1p1h.h"
#endif

#define PNOFSET 10000

double Random(void) {
  double a =  (double)random()/(double)RAND_MAX;
  return a; 
}


bool N1p1h::ReadConfig(const char *dir){

#ifndef WITH_NEUT
  char filename[1024];

  snprintf(filename,sizeof(filename),"%s/Configuration",dir);

  std::cout << " Reading configuration from " << filename << std::endl; 

  std::ifstream infile(filename); 

  if( infile.fail() ) return false; 

  char item[256];
  double val; 


  do {
    char line[1024];

    infile.getline(line,1024);
    
    if( line[0] == '#' ) continue; 
    
    sscanf(line,"%s %lf",item,&val);

    if( strncmp(item,"MA",2) ==0 ) xmag = val;
    else if ( strncmp(item,"IRPA",4) ==0 ) irpa = val;
    else if ( strncmp(item,"RFG",3) ==0 ) RFG = val;
    else if ( strncmp(item,"BIND",4) ==0 ) ApplyBindEnergy = val;
    else if ( strncmp(item,"FP0IN",5) ==0 ) RPAfp0in = val;
    else if ( strncmp(item,"FP0EX",5) ==0 ) RPAfp0ex = val;
    else if ( strncmp(item,"FSTAR",5) ==0 ) RPAfstar = val;
    else if ( strncmp(item,"F",1) ==0 )     RPAf = val;
    else if ( strncmp(item,"PILAMBDA",8) ==0 )   RPApilambda = val;
    else if ( strncmp(item,"CR0",3) ==0 )        RPAcr0 = val;
    else if ( strncmp(item,"RHOLAMBDA",9) ==0 )  RPArholambda = val;
    else if ( strncmp(item,"GP",2) ==0 )     RPAgp = val;
    else if ( strncmp(item,"XMPI",4) ==0 )   RPAxmpi = val;
    else if ( strncmp(item,"XMRHO",5) ==0 )  RPAxmrho = val;
    else if ( strncmp(item,"IREL",4) ==0 )   RPAirel = val;
    else if ( strncmp(item,"EBIND",5) ==0 )   Ebind = val;
    //    else if ( strncmp(item,"FFTYPE",6) ==0 )   FFtype = val;
    
  } while( ! infile.eof() );

#else
  // Copy parameres from the neut card
  
  xmag = nemdls_.xmaqe;
  irpa = nievesqepar_.nvqerpa;
  if (nievesqepar_.nvqerfg == 0){
	RFG  = false;
  }else{
	RFG  = true;
  }
  if (nievesqepar_.nvqebind == 0){
	ApplyBindEnergy = false;
  }else{
	ApplyBindEnergy = true;
  }

  RPAfp0in     = nievesqepar_.xnvrpafp0in;
  RPAfp0ex     = nievesqepar_.xnvrpapf0ex;
  RPAfstar     = nievesqepar_.xnvrpafstar;
  RPAf         = nievesqepar_.xnvrpaf;
  RPApilambda  = nievesqepar_.xnvrpapilambda;
  RPAcr0       = nievesqepar_.xnvrpacr0;
  RPArholambda = nievesqepar_.xnvrparholambda;
  RPAgp        = nievesqepar_.xnvrpagp;
  RPAxmpi      = nievesqepar_.xnvrpaxmpi;
  RPAxmrho     = nievesqepar_.xnvrpaxmrho;
  RPAirel      = nievesqepar_.xnvrpairel;
  //  FFtype       = nievesqepar_.fftype;
  
#endif

  card_params[param_names[0]]   = xmag;
  if (irpa){
	card_params[param_names[1]] = 1;
  }else{
	card_params[param_names[1]] = 0;
  }
  if (RFG){
	card_params[param_names[2]] = 1;
  }else{
	card_params[param_names[2]] = 0;
  }
  if (ApplyBindEnergy){
	card_params[param_names[3]] = 1;
  }else{
	card_params[param_names[3]] = 0;
  }
  card_params[param_names[4]]   = RPAfp0in;
  card_params[param_names[5]]   = RPAfp0ex;
  card_params[param_names[6]]   = RPAf;
  card_params[param_names[7]]   = RPAfstar;
  card_params[param_names[8]]   = RPApilambda;
  card_params[param_names[9]]   = RPAcr0;
  card_params[param_names[10]]  = RPArholambda;
  card_params[param_names[11]]  = RPAgp;
  card_params[param_names[12]]  = RPAxmpi;
  card_params[param_names[13]]  = RPAxmrho;
  card_params[param_names[14]]  = RPAirel;
  card_params[param_names[15]]  = 0.;
  //  card_params[param_names[16]]  = FFtype;

  if( !ApplyBindEnergy ) 
    std::cout << " Bind energy ignored " << std::endl; 
  else
    std::cout << " Bind energy used " << std::endl;

  if(Ebind != 0)
    std::cout << "Nucleus created in an excited state (GeV): " << Ebind << std::endl;
  
  return true;
}

N1p1h::N1p1h(){
  return;
}


N1p1h::N1p1h(std::string directory, bool d){
  Initialize(directory,d);
  InitializeAllNuclei(); 
}

N1p1h::N1p1h(bool d){
  Initialize("Tables",d);
  InitializeAllNuclei(); 
}

N1p1h::N1p1h(std::string directory, int nuclei, bool d  ) {
  Initialize(directory,d);
  InitializeNucleus(nuclei,true);
}

void
N1p1h::Initialize(bool d) {
  Initialize("Tables", d);
  InitializeAllNuclei();
}

void
N1p1h::Initialize(std::string directory, bool d){

  int i;

  xmag     = 1.05; 
  irpa     = 1; //1->NUCLEAR EFFECT
  RFG      = false; 
  ApplyBindEnergy = true; // Default is always on.  
  RPAfp0in = 0.33; 
  RPAfp0ex = 0.45; 
  RPAf     = 1; 
  RPAfstar = 2.13; 
  RPApilambda = 1200.;
  RPAcr0      = 2.0; 
  RPArholambda = 2500.;
  RPAgp        =  0.63; 
  RPAxmpi = 139.57;
  RPAxmrho = 770.0; 
  RPAirel  = 1;
  //  FFtype = 1;
 
  Ebind=0;
  
  debug = d; 

  static const char param_name_const[N1P1H_CARD_PARAMS][32]=
	{"#MA:\0",         "#RPA:\0",     "#RFG:\0",         "#BIND:\0",
	 "#RPAfp0in:\0",   "#RPAfp0ex:\0","#RPAf:\0",   	 "#RPAfstar:\0",
	 "#RPApilambda:\0","#RPAcr0:\0",  "#RPArholambda:\0","#RPAgp:\0",
	 "#RPAxmpi:\0",    "#RPAxmrho:\0","#RPAirel:\0",     "#Enubin:\0"};
	 //	 "#FFtype:\0"};

  for ( i = 0 ; i < N1P1H_CARD_PARAMS; i++ ){
	snprintf(param_names[i],32,"%s",param_name_const[i]);
  }

  ReadConfig(directory.c_str()); 

  snprintf(dirtable,sizeof(dirtable),
		   "%s",directory.c_str());

  Enubin = 0.1;   // 100 MeV 
  Emax   = 10.5;   // 10 GeV

  RFG = false; 
 
  NRC= 1;//0->loop on R /  1->Random R  /  2->free nucleon
  NPC= 1;//0->dur     / 1->dur1(integral) /  2->durt(prand)
    
  // Parameters 
  
  GFermi = 1.1664e-5;   // Gfermi is in GeV^(-2))
  facconv = 0.19733*0.19733*1.e+15;

  unit = pow(10,-41.);

  int    anti     = 1;//1->neutrino;-1->antineutrinos
  double type     = 12;//12,16,40.;
  int contrlepton = 1;//1->muons,0->electron
  double rm; 

  initialization_(&xmag,&irpa,&rm,&anti,&type,&contrlepton);

  constantsinitialization_(&xmag,&irpa);

  rpaparameters_(&RPAfp0in,&RPAfp0ex,&RPAf,&RPAfstar,&RPApilambda,&RPAcr0,&RPArholambda,&RPAgp,&RPAxmpi,&RPAxmrho,&RPAirel);

  leptonmass[11] = 0.000510998910;  // lepton mass in GeV
  leptonmass[13] = 0.105658357; // muon mass in GeV
  leptonmass[15] = 1.77682; // tau mass in GeV 

  leptonmass[12] = 0.;  // neutrino mass in GeV
  leptonmass[14] = 0.; //  neutrino mass in GeV
  leptonmass[16] = 0.; //  neutrino mass in GeV 

  daughter[-12] = -11 ;
  daughter[12] = 11;
  daughter[-14] = -13;
  daughter[14] = 13;
  daughter[-16] = -15;
  daughter[16] = 15;

  idneutron = 2112;
  idproton = 2212;
  
  protonmass = 0.93827208;
  neutronmass = 0.93956542;
  
  neutrinoIdlist.push_back(-12);
  neutrinoIdlist.push_back(12);
  neutrinoIdlist.push_back(-14);
  neutrinoIdlist.push_back(14);
  neutrinoIdlist.push_back(-16);
  neutrinoIdlist.push_back(16);

  std::cout << " N1p1h creating nuclei " << std::endl; 
  
  return;
}

void N1p1h::InitializeAllNuclei(void) {
  InitializeNucleus(1,true);
  InitializeNucleus(12,true);
  InitializeNucleus(16,true);
  InitializeNucleus(27,true);
  InitializeNucleus(28,true);
  InitializeNucleus(40,true);
  InitializeNucleus(56,true);
  InitializeNucleus(208,true);
  return; 
}

void N1p1h::InitializeNucleus(int id, bool CInt) {

  if( Nuclei[id]  ) // if exist delete to create a new one. 
    delete Nuclei[id];

  if(id == 1){
    if( CInt )//force H integral calculation
      ComputeHintegrals(1);
    return;
  }
  
  Nuclei[id] = new Nucleus(id);
  
  if( CInt )   // force the integral calculation 
    ComputeIntegrals(id);
  
  return; 
}


void N1p1h::ComputeIntegrals(int nuclei){
  
  CheckNuclei(nuclei);  

  double hbarc = 197.3269602;
  
  int binmax = Emax/Enubin;
  //  double rmax = Nuclei[nuclei]->FindRMaximumDensity();
  int id,lid;
  double ml;

  double Enu,Enu2;
  double Xsect,Xsect2,MaxXsect;
  double maxcos,maxtl;
  int ntot;
  double rmax,ldpmax,vcmax,q2max,emumaxemumin;
  double qval=0.;
  double lr,ang,xcos,tmommax;
  //double pmommax; 
  double vc; 
  double q0,dq,Elloc,Plloc,coseno,q2;
  double Emumax,Emumin;	
  double Tlepmax,Tlepmin,tl;
  double El;
  double tmommin;
  double ldp,ldpsaved;
  double correctionangle;
  double xs,xs2,dd;
  double cosp;
  double xsect2d;

  if( !ReadFromTable(nuclei) ) {
    
    double pi = acos(-1.);
    
    for( unsigned int il = 0; il < neutrinoIdlist.size(); il++ ) { // Looping for all neutrinos and antineutrinos.
      id = neutrinoIdlist[il];
      lid = daughter[id];
      Precompindx pindx(id,nuclei);
      std::vector<double> integral;
      std::vector<double> maximal;
      ml = leptonmass[abs(lid)];
      double ml2 = ml*ml;

      SelectNucleusLepton(lid,nuclei); 
      
      //      binmax = 20; 
      
      std::cout << " Neutrino " << id << " Lepton  " << lid << " mass " << ml << std::endl; 
      
      for( int i = 0; i < binmax; i++ ) {  //looping over the neutrino energy
	Enu = (double)i*Enubin+Enubin*0.5;
	//Enu=1.0;
	Enu2 = Enu*Enu;
	Xsect = Xsect2 = 0.;
	MaxXsect = 0.;
	maxcos   = 0.;
	maxtl    = 0.; 
	xsect2d  = 0.;
	
	ntot = 1000000; 

	rmax = 0.;
	ldpmax = 0.;
	vcmax = 0.; 
	q2max =0.;
	emumaxemumin = 0.; 

	//radial position
	double Rmax = Nuclei[nuclei]->GetMaximumRadious();

	if(ApplyBindEnergy)
	  qval= (Nuclei[nuclei]->GetQvalueGeV(id)) + Ebind;
	else
	  qval=0;
	
	for( int  iloop = 0 ; iloop < ntot; iloop++ ){
	  //angle neutrino-lepton
	  ang = Random()*pi;
	  xcos = cos(ang);
	  
	  double intfactor4 = 1.; 
	  //double intfactor2 = 1.; 
	  
	  correctionangle = sqrt(1.-xcos*xcos)*pi;
	  
	  //intfactor2 *= correctionangle;
	  intfactor4 *= correctionangle;
	  
	  lr = Random()*Rmax;
	  
	  intfactor4 *= Rmax;
	  
	  //nucleus simulation
	  if( irpa ) 
	    vcinterpolate_(&lr,&vc); //vc is in GeV
	  else 
	    vc = 0.; 

	  tmommax = Nuclei[nuclei]->GetMaximumFermi(id,lr);
	  	    
	  //lepton kinetic energy
	  Emumax=EmuMax(Enu,xcos,ml,id,tmommax,vc,qval);

	  if( Emumax > Enu ) Emumax = Enu;
	  	  	  
	  Tlepmax = Emumax-ml;

	  if( Tlepmax < 0.) continue;
	  
	  Emumin=EmuMin(Enu,xcos,ml,id,tmommax,vc,qval);	  

	  if( Emumin < ml ) Emumin = ml; 

	  Tlepmin = Emumin-ml;

	  if( Tlepmin > Tlepmax ) continue;

	  tl= Tlepmin+(Tlepmax-Tlepmin)*Random();
	  
	  intfactor4 *= (Tlepmax-Tlepmin);
	  //	  intfactor2 *= (Tlepmax-Tlepmin);
	  
	  El = ml+tl;
	  Elloc=El-vc;

	  if( Elloc-ml <= 0.) continue;	
  
	  //end lepton kinetic energy (note that min is not used/ !over restriction)
	  
	  //local calculation !!!!!!!!!!!for target momentum determination
	  double cosang = cos(ang/2.);
	  coseno =2.*cosang*cosang-1.;
	  q0=Enu-El-qval;
	  if( q0 < 0. ) continue;
	  
	  double Plloc2 = Elloc*Elloc-ml2;

	  Plloc=sqrt(Plloc2);
	  
	  double dq2 = Enu2+Plloc2-2*Enu*Plloc*coseno;

	  dq=sqrt(dq2);
	  
	  q2=q0*q0-dq2;

	  //End local
	  
	  //Target hadron momentum !!IN GeV
	  if(id < 0.)
	    tmommin = targetminmon_(&q0,&dq,&protonmass,&neutronmass);
	  else if(id > 0.)
	    tmommin = targetminmon_(&q0,&dq,&neutronmass,&protonmass);
	  	  
	  if(tmommin > tmommax)  continue; 

	  if(std::isnan(tmommin))tmommin=0;

	  ldp = tmommin+Random()*(tmommax-tmommin);
	  
	  ldpsaved = ldp;
	  
	  intfactor4 *= (tmommax-tmommin)*1000./hbarc;
	  
	  double Ep;
	  if(id<0.){
	    Ep = sqrt(protonmass*protonmass+ldp*ldp);
	    cosp = 0.5*((Ep+q0)*(Ep+q0)-(neutronmass*neutronmass+ldp*ldp+dq*dq))/(dq*ldp);
	  }else if(id>0.){
	    Ep = sqrt(neutronmass*neutronmass+ldp*ldp);
	    cosp = 0.5*((Ep+q0)*(Ep+q0)-(protonmass*protonmass+ldp*ldp+dq*dq))/(dq*ldp);	    
	  }
	  
	  if(cosp > 1. || cosp < -1.)continue;

	  //end target mom
	  
	  // Four differential
	  
	  double xsl = FourDifferential(id,nuclei,Enu,tl,xcos,lr,ldp,vc);

	  if( xsl < 0. || std::isnan(xsl)  ) xsl=0.;

	  xs =  xsl*intfactor4;
	  xs2 =  xsl*xsl*intfactor4;

	  if( xs > MaxXsect ) { MaxXsect = xs; maxcos = xcos; maxtl = tl; rmax = lr; ldpmax = ldpsaved; vcmax = vc;q2max=q2; emumaxemumin = Emumax-Emumin;
	  }

	  Xsect   += xs;
	  Xsect2   += xs2;
	  //	  xsect2d  += dd;
	}

	Xsect   /= (double)ntot;
	Xsect2  /= (double)ntot;
	//	xsect2d /= (double)ntot; 
	
	integral.push_back(Xsect);
	maximal.push_back(MaxXsect);
	
        if( debug ) 
	std::cout << " nuclei " << nuclei << " Lepton " << lid << " Enu " << Enu << " xsect " << Xsect << "+-" << sqrt(Xsect2-Xsect*Xsect)/sqrt(ntot) << " max " << MaxXsect << " ( " << maxcos << " , " << maxtl << " , " << rmax << " , " << q2max  << " , " <<  emumaxemumin << " ) " << std::endl; 

      }
      
      
      IntCrossSection[pindx] = integral;
      MaxCrossSection[pindx] = maximal; 
      
      
    }

    RecomputeMaximum(nuclei); 

    WriteTable(nuclei);
  }

return; 
}


bool N1p1h::ReadFromTable(int nuclei) {
  
  char Table[256]; 

  snprintf(Table,sizeof(Table),
		   "%s/CrossSection_%d.dat",dirtable,nuclei); 

  std::ifstream infile(Table); 

  if( infile.fail() ) return false; 

  std::cout << " Reading data from " << Table << std::endl; 

  std::vector<double> elintegral;
  std::vector<double> elmaximal;
  std::vector<double> aelintegral;
  std::vector<double> aelmaximal;
  std::vector<double> muonintegral;
  std::vector<double> muonmaximal;
  std::vector<double> amuonintegral;
  std::vector<double> amuonmaximal;
  std::vector<double> tauintegral;
  std::vector<double> taumaximal;
  std::vector<double> atauintegral;
  std::vector<double> ataumaximal;

  double enu; 

  double eltot, elmax, aeltot, aelmax, mutot, mumax, amutot, amumax, tautot, taumax, atautot, ataumax; 

  bool all_ok;
  
  char tmpstr[1024],comment[1024],comment2[1024];
  int  i,ret;
  float tmpval,diff;

  while(1){
	while(1){
	  infile.getline(&(tmpstr[0]),sizeof(tmpstr));
	  if (infile.eof()){
		return false;
	  }
	  if (strstr(tmpstr,"#PARAMSTART")){
		break;
	  }
	}
	for( i = 0 ; i < N1P1H_CARD_PARAMS ; i++ ){
	  card_params_chk[param_names[i]] = false;
	}

	while(1){
	  infile.getline(&(tmpstr[0]),sizeof(tmpstr));
	  if (infile.eof()){
		return false;
	  }
	  if (strstr(tmpstr,"#TABLESTART")){
		break;
	  }
	  sscanf(tmpstr,"%s %f",comment,&tmpval);
	  diff = card_params[comment] - tmpval;
	  if ( card_params[comment] > 0 ){
		diff = diff / card_params[comment];
	  }
	  if (strstr(comment,"Enubin")){
		diff = 0.;
		Enubin = tmpval;
	  }
	  if ( fabs(diff)<0.01 ){
		card_params_chk[comment] = true;
	  }else{
		break;
	  }
	}
	all_ok = true;
	for( i = 0 ; i < N1P1H_CARD_PARAMS; i++ ){
	  if (card_params_chk[param_names[i]] != true){
		all_ok = false;
		break;
	  }
	}
	if (all_ok == true){
	  break;
	}
  }

  while( 1 ){
	infile.getline(&(tmpstr[0]),sizeof(tmpstr));
	if (tmpstr[0]=='#'){
	  break;
	}
	ret = sscanf(tmpstr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
				 &enu,    &eltot, &elmax,  &aeltot, &aelmax, 
				 &mutot,  &mumax, &amutot, &amumax, &tautot,
				 &taumax, &atautot,&ataumax );
	if (ret != 13){
	  break;
	}
    elintegral.push_back(eltot);
    elmaximal.push_back(elmax);
    aelintegral.push_back(aeltot);
    aelmaximal.push_back(aelmax);
    muonintegral.push_back(mutot);
    muonmaximal.push_back(mumax);
    amuonintegral.push_back(amutot);
    amuonmaximal.push_back(amumax);
    tauintegral.push_back(tautot);
    taumaximal.push_back(taumax);
    atauintegral.push_back(atautot);
    ataumaximal.push_back(ataumax);

  }

  IntCrossSection[Precompindx(12,nuclei)] = elintegral;
  MaxCrossSection[Precompindx(12,nuclei)] = elmaximal; 
  IntCrossSection[Precompindx(-12,nuclei)] = aelintegral;
  MaxCrossSection[Precompindx(-12,nuclei)] = aelmaximal; 
  IntCrossSection[Precompindx(14,nuclei)] = muonintegral;
  MaxCrossSection[Precompindx(14,nuclei)] = muonmaximal; 
  IntCrossSection[Precompindx(-14,nuclei)] = amuonintegral;
  MaxCrossSection[Precompindx(-14,nuclei)] = amuonmaximal; 
  IntCrossSection[Precompindx(16,nuclei)] = tauintegral;
  MaxCrossSection[Precompindx(16,nuclei)] = taumaximal; 
  IntCrossSection[Precompindx(-16,nuclei)] = atauintegral;
  MaxCrossSection[Precompindx(-16,nuclei)] = ataumaximal; 

  return true;
}



void N1p1h::RecomputeMaximum(int nuclei) {
  bool changed = false; 

  do{
    changed = false; 
    for( unsigned int i = 1; i < MaxCrossSection[Precompindx(12,nuclei)].size()-1 ; i++ ) {
      for( unsigned int il = 0; il < neutrinoIdlist.size(); il++ ) { // Looping for all neutrinos and antineutrinos.
	int id = neutrinoIdlist[il];
       
	if( (MaxCrossSection[Precompindx(id,nuclei)][i-1]+ MaxCrossSection[Precompindx(id,nuclei)][i+1])/2. -  MaxCrossSection[Precompindx(id,nuclei)][i] >  MaxCrossSection[Precompindx(id,nuclei)][i]*1.001 ) {
	  std::cout << MaxCrossSection[Precompindx(id,nuclei)][i] ; 
	  MaxCrossSection[Precompindx(id,nuclei)][i] =  (MaxCrossSection[Precompindx(id,nuclei)][i-1]+ MaxCrossSection[Precompindx(id,nuclei)][i+1])/2. ; 
	  changed = true;
	  if( debug ) 
	  std::cout << " >>> " <<  MaxCrossSection[Precompindx(id,nuclei)][i] << std::endl; 
	}  
      }
    }
  } while( changed ); 
    
}


bool N1p1h::WriteTable(int nuclei) {
  
  char Table[256]; 

  snprintf(Table,sizeof(Table),
		   "%s/CrossSection_%d.dat",dirtable,nuclei); 
  
  std::ofstream offile(Table,std::ofstream::app); 
  
  if( offile.fail() ) return false; 
  offile << "#PARAMSTART"     << std::endl;

  offile << "#MA: " << xmag << std::endl;


  offile << "#RPA:          " << irpa << std::endl;

  if (RFG){
  offile << "#RFG:            1" << std::endl;
  }else{
  offile << "#RFG:            0" << std::endl;
  }
  if (ApplyBindEnergy){
  offile << "#BIND:           1" << std::endl;
  }else{
  offile << "#BIND:           0" << std::endl;
  }
  offile << "#RPAfp0in:     " << RPAfp0in     << std::endl;
  offile << "#RPAfp0ex:     " << RPAfp0ex     << std::endl;
  offile << "#RPAf:         " << RPAf         << std::endl;
  offile << "#RPAfstar:     " << RPAfstar     << std::endl;
  offile << "#RPApilambda:  " << RPApilambda  << std::endl;
  offile << "#RPAcr0:       " << RPAcr0       << std::endl;
  offile << "#RPArholambda: " << RPArholambda << std::endl;
  offile << "#RPAgp:        " << RPAgp        << std::endl;
  offile << "#RPAxmpi:      " << RPAxmpi      << std::endl;
  offile << "#RPAxmrho:     " << RPAxmrho     << std::endl;
  offile << "#RPAirel:      " << RPAirel      << std::endl;
  offile << "#Enubin:       " << Enubin       << std::endl;
  offile << "#TABLESTART"     << std::endl;

  for( unsigned int i = 0; i < IntCrossSection[Precompindx(12,nuclei)].size() ; i++ ) {

    double enu = ((double)i*Enubin+Enubin*0.5);

    offile << enu << "  " 
           << IntCrossSection[Precompindx(12,nuclei)][i]  << "  " <<  MaxCrossSection[Precompindx(12,nuclei)][i] << "  " 
	   << IntCrossSection[Precompindx(-12,nuclei)][i] << "  " << MaxCrossSection[Precompindx(-12,nuclei)][i] << "  " 
	   << IntCrossSection[Precompindx(14,nuclei)][i]  << "  " << MaxCrossSection[Precompindx(14,nuclei)][i]  << "  " 
	   << IntCrossSection[Precompindx(-14,nuclei)][i] << "  " << MaxCrossSection[Precompindx(-14,nuclei)][i] << "  " 
	   << IntCrossSection[Precompindx(16,nuclei)][i]  << "  " << MaxCrossSection[Precompindx(16,nuclei)][i]  << "  " 
	   << IntCrossSection[Precompindx(-16,nuclei)][i] << "  " << MaxCrossSection[Precompindx(-16,nuclei)][i] << " \n";
   
  }
  offile << "#"     << std::endl;

  return true;
}



double N1p1h::IntegralCrossSection(int id,int nuclei,double Enu){

  CheckNuclei(nuclei);  
  
  int ibin = (int) ((Enu+Enubin*0.001)/Enubin); 

  Precompindx pindx(id,nuclei);  
  
  double a0 = IntCrossSection[pindx][ibin]; 
  double a1 = IntCrossSection[pindx][ibin+1];

  int binmax = Emax/Enubin;

  double cross_section; 

  if( ibin+1  < binmax ) 
     cross_section = (a1-a0)/Enubin*(Enu-Enubin*(double)ibin)+a0;  
  else if ( ibin >= binmax  )
     cross_section = IntCrossSection[pindx][binmax-1]; 
  else if ( ibin < binmax ) 
    cross_section = a0;
  else
    cross_section = 0.0;

  return cross_section; 
}

double N1p1h::CrossSectionmax(int id,int nuclei,double Enu){

  CheckNuclei(nuclei);  

  int ibin = (int) ((Enu+Enubin*0.001)/Enubin); 

  Precompindx pindx(id,nuclei);  
  
  double a0 = MaxCrossSection[pindx][ibin]; 
  double a1 = MaxCrossSection[pindx][ibin+1];

  double Max_cross_section; 

  int binmax = Emax/Enubin;

  if( ibin+1  < binmax ) 
    Max_cross_section = (a1-a0)/Enubin*(Enu-Enubin*(double)ibin)+a0;  
  else if ( ibin >= binmax  )
    Max_cross_section = MaxCrossSection[pindx][binmax-1];  
  else if ( ibin < binmax ) 
    Max_cross_section = a0;
  else 
    Max_cross_section = 0.0;
    
  return Max_cross_section*1.; 
}

double   N1p1h::DoubleDifferential(int id,int nuclei, double Enu,double TLep,double xcos) {

  if( nuclei == 1 ){
    int lid = daughter[id];
    SelectNucleusLepton(lid,nuclei);
    double xs=sigthfree_(&Enu,&TLep,&xcos)*unit;
    cosH=xcos;
    return xs; 
  }else {
    NPC = NRC = 0; 
    int lid = daughter[id];
    SelectNucleusLepton(lid,nuclei); 
    return sigthnacho_(&Enu,&TLep,&xcos)*unit;
  }
}


double   N1p1h::FourDifferential(int id,int nuclei, double Enu,double TLep,double xcos,double &R, double &dP, double &vc) {
 
  if( nuclei == 1 ) {
    int lid = daughter[id];
    SelectNucleusLepton(lid,nuclei); 
    return sigthfree_(&Enu,&TLep,&xcos)*unit; 
  }
  else {
    NRC = 10;NPC = 10;
    int lid = daughter[id]; 
    SelectNucleusLepton(lid,nuclei);  
    return sigthbruno_(&Enu,&TLep,&xcos,&NRC,&NPC,&R,&dP,&vc)*unit;
  }
}


void N1p1h::GenerateVectors(int id,int nuclei,double pnu[4],  double p[4][4], int idpart[4], int parent[4], double &R, double &xsec){

  double hbarc = 197.3269602;

  Precompindx pindx(id,nuclei);
  
  double Enu = pnu[0];
  
  double Enu2 = Enu*Enu; 

  double dP,vc; 

  int lid = daughter[id];

  double mass = leptonmass[abs(lid)];
  double mass2 = mass*mass; 

  double Tlepmax = Enu-mass;
  
  if( ApplyBindEnergy ) 
    Tlepmax -= ((Nuclei[nuclei]->GetQvalueGeV(id)) + Ebind);  
  
  int ibin = (int) ((Enu+Enubin*0.001)/Enubin);

  double cmax = CrossSectionmax(id,nuclei,Enu);
    
  double cx,val; 

  int    it = 0; 
  double pi = acos(-1.);

  numberofiterationslastevent = 0; 

  //radial position
  double Rmax = Nuclei[nuclei]->GetMaximumRadious();

  SelectNucleusLepton(lid,nuclei);  

  //  Nuclei[nuclei]->DumpInfo();

  double qval; 
  if(ApplyBindEnergy)
    qval= (Nuclei[nuclei]->GetQvalueGeV(id)) + Ebind;
  else
    qval=0;
  
  int r0 = 0; 
  int r1 = 0; 
  int r2 = 0; 
  int r3 = 0;
  int r4 = 0; 
 
  do {

    double intfactor4 = Rmax; 

    it++; 
    double ang = Random()*pi;
    double xcos = cos(ang); 

    intfactor4 *= sqrt(1.-xcos*xcos)*pi;

    R = Random()*Rmax;
    
    double vc = 0.;
 
    //nucleus simulation
    if( irpa ) 
      vcinterpolate_(&R,&vc); //vc is in GeV
    
    //lepton kinetic energy
    double tmommax = Nuclei[nuclei]->GetMaximumFermi(id,R);	  
    double Emumax = EmuMax(Enu,xcos,mass,id,tmommax,vc,qval);
    double Emumin = EmuMin(Enu,xcos,mass,id,tmommax,vc,qval);
    
    if( Emumax > Enu ) Emumax = Enu;
    if( Emumin < mass ) Emumin = mass; 
    
    // Emumax = Enu; Emumin=mass; 
    
    double Tlepmax = Emumax-mass;
    double Tlepmin = Emumin-mass;

    // Keep the values for debugging 

    TlepmaxKept = Tlepmax; 
    TlepminKept = Tlepmin; 
        
    if( Tlepmax < 0. || Tlepmin > Tlepmax ) { 
      cx = 0; val = 1.;
      r0++;
      continue;
    } 
    
    double Tlep = Tlepmin+(Tlepmax-Tlepmin)*Random();

    intfactor4 *= (Tlepmax-Tlepmin);
    
    double El = mass+Tlep;
    double Elloc=El-vc;

    if(Elloc-mass <= 0.){
      cx = 0; val = 1.;
      r1++;
      continue;
    }

    double Plloc2 = Elloc*Elloc-mass2;
    if(Plloc2 < 0.){
      cx = 0; val = 1.;
      r1++;
      continue;
    }    
    double Plloc = sqrt(Plloc2);
    
    // End lepton kinetic energy (note that min is not used/ !over restriction
    
    // Local calculation !!!!!!!!!!!for target momentum determination
    double cosang = cos(ang/2.);
    double coseno = 2.*cosang*cosang-1.;
    double q0     = Enu-El-qval;
    if( q0 < 0. ){
      cx = 0; val = 1.;
      continue;
    }
    
    double dq     = sqrt(Enu2+Plloc2-2.*Enu*Plloc*coseno);
    
    double tmommin = 0; 

    //Target hadron momentum !!IN GeV
    if(id < 0.) 
      tmommin = targetminmon_(&q0,&dq,&protonmass,&neutronmass);
    else if(id > 0.)
      tmommin = targetminmon_(&q0,&dq,&neutronmass,&protonmass);
    
    if(std::isnan(tmommin))tmommin=0;

    tmomminkept = tmommin; 
    
    if( tmommin >= tmommax ){
      cx = 0; val = 1;
      r2++;
      continue;
    }

    dP = tmommin+Random()*(tmommax-tmommin);

    double Ep,cosp;
    if(id<0.){
      Ep = sqrt(protonmass*protonmass+dP*dP);
      cosp = 0.5*((Ep+q0)*(Ep+q0)-(neutronmass*neutronmass+dP*dP+dq*dq))/(dq*dP);
    }else if(id>0.){
      Ep = sqrt(neutronmass*neutronmass+dP*dP);
      cosp = 0.5*((Ep+q0)*(Ep+q0)-(protonmass*protonmass+dP*dP+dq*dq))/(dq*dP);	    
    }

    if(cosp < -1. || cosp > 1.){
      cx = 0; val = 1;
      continue;
    }

    intfactor4 *= (tmommax-tmommin)*1000./hbarc;
        
    cx =  FourDifferential(id,nuclei,Enu,Tlep,xcos,R,dP,vc)*intfactor4;

    if( cx <= 0.0 ) { 
      cx = 0;
      val = 1.;
      r3++;
      continue; 
    }
    
    if( cx > cmax*1.2 ) { 
      std::cout << " Enu " << Enu << " max " << cmax << " < " << cx << "( " <<  MaxCrossSection[pindx][ibin] << " , " <<  MaxCrossSection[pindx][ibin+1] << " ) " << std::endl; 
      cmax = cx; 
    }

    val = cmax*Random(); 

  } while ( cx < val && it < 50e+4 ); 


  //  std::cout << (double)r0/(double)it << "  " << 100.*(double)r1/(double)it << "  " << 100.*(double)r2/(double)it << "  " << 100.*(double)r3/(double)it << "  " << std::endl;

  xsec = cx;

  if( xsec == 0  ) std::cout << it << std::endl; 

  double z[3];

  z[0] = pnu[1]/Enu;
  z[1] = pnu[2]/Enu;
  z[2] = pnu[3]/Enu;

  comptmomentum_(z);
  
  for(int i=0;i<4;i++){
    p[0][i]=fourvectors_.vPnu[i];
    p[1][i]=fourvectors_.vPn[i];
    p[2][i]=fourvectors_.vPmu[i];
    p[3][i]=fourvectors_.vPp[i];
  }

  parent[0] = 0;
  parent[1] = 0;
  parent[2] = 1;
  parent[3] = 2;

  idpart[0] = id;
  if( id < 0 ) {
    idpart[1] = idproton;
    idpart[3] = idneutron; 
  }
  else {
    idpart[3] = idproton;
    idpart[1] = idneutron; 
  }
  idpart[2] = daughter[id]; 
  
  numberofiterationslastevent = it; 

  return; 
}

void N1p1h::CheckNuclei(int nuclei){
  if( Nuclei[nuclei] == NULL && nuclei != 1 ) {
	
    std::cout << " N1p1h Error: nuclei " << nuclei << " not available " << std::endl;
    exit(1);
    /*
    std::cout << " N1p1h Initialize for nuclei " << nuclei << std::endl;
	Nuclei[nuclei] = new Nucleus(nuclei);
	ComputeIntegrals(nuclei);*/
  }
  
}

double N1p1h::EmuMax(double E,double xcos,double ml, int id,double pfermi, double vc,double qval){
  //double mp = protonmass;
  //double mn = neutronmass;
  double mt,ms;
  if( id < 0.){
    mt=protonmass;
    ms=neutronmass;
  }else{
    mt=neutronmass;
    ms=protonmass;
  }

  double E2=E*E;
  double ml2=ml*ml; 
  double mt2=mt*mt;
  double ms2=ms*ms;
  double pfermi2=pfermi*pfermi; 

  //  double Enf = sqrt(pfermi*pfermi+mn*mn); 
  //  double Epf = sqrt(pfermi*pfermi+mp*mp); 

  double Emax = E - qval + sqrt(pfermi2+mt2);
  double Emumax = Emax-sqrt(pfermi2+ms2);
  
  // Solve Equations...Newton Method
 
  double emin = ml+vc;
  double emax = Emumax; 

  double d1,d2; 
  
  double delta = 0.; 
  int i=0;
  double Ep,pmu,p;
  double ecos2 = 2.*E*xcos;


    
  Ep = Emax - emin;
  double pmu2;
  pmu2 = (emin-vc)*(emin-vc)-ml2; 
  if ( pmu2 > 0 ){
	pmu = sqrt(pmu2);
  } else {
	pmu2=0.;
	pmu = 0;
  }

  if( std::isnan(pmu) || (emin-vc) < 0. ) { pmu2=pmu=0; }
  p = sqrt(E2+pmu2-pmu*ecos2)-pfermi;
  
  d1 = Ep-sqrt(p*p+ms2);
  
  Ep = Emax - emax;
  pmu2=(emax-vc)*(emax-vc)-ml2;
  if ( pmu2 > 0 ){
	pmu = sqrt(pmu2);
  }else{
	pmu2 = 0;
	pmu = 0;
  }
  
  p = sqrt(E2+pmu2-pmu*ecos2)-pfermi;
  
  if( p < pfermi ) { emin = emax = Emumax; return emax; }
  
  d2 = Ep-sqrt(p*p+ms2);
  
  do{

    if( d1*d2 > 0. ){
      emin = emax = Emumax;
      i=99.;
    }
    else if( d1*d2 < 0.0 ) { 
      double emid = (emax+emin)/2.;
      
      Ep = Emax - emid;
      pmu2 =(emid-vc)*(emid-vc)-ml2;
      pmu = sqrt(pmu2);
      p = sqrt(E2+pmu2-pmu*ecos2)-pfermi;

      double d3 = Ep-sqrt(p*p+ms2);
      
      if( d1*d3 < 0. ) {
	emax = emid;
	d2 = d3; 
      }
      else if( d2*d3 < 0. ){
	emin = emid;
	d1 = d3; 
      }
    }


    delta = (emax-emin)/(emax+emin)*2.; 
    i++;
  } while ( delta > 1.e-4 && i < 1e+2 ); 

  return emax; 
}

double N1p1h::EmuMin(double E,double xcos,double ml, int id,double pfermi, double vc, double qval){
  //double mp = protonmass;
  //double mn = neutronmass;
  double mt,ms;
  if( id < 0.){
    mt=protonmass;
    ms=neutronmass;
  } else {
    mt=neutronmass;
    ms=protonmass;
  }

  double pfermi2 = pfermi*pfermi;
  double mt2 = mt*mt; 
  double ms2 = ms*ms;
  double ml2 = ml*ml; 

  double E2 = E*E; 

  double Emax = E -qval + sqrt(pfermi2+mt2);
  double Emumax = Emax - sqrt(pfermi2+ms2);
 
  // Solve Equations...Newton Method
  
  double emin = ml+vc;
  double emax = Emumax; 

  double d1,d2,d3; 
  
  double delta = 0.; 
  int i=0;
  double Ep,pmu,p;

  double ecos2 = 2*E*xcos;

   Ep = Emax - emin;
   double pmu2;
   pmu2 = (emin-vc)*(emin-vc)-ml2; 
   if ( pmu2 > 0 ){
	 pmu = sqrt(pmu2);	 
   } else {
	 pmu2=0.;
	 pmu =0.;
   }

   if(std::isnan(pmu)){pmu2=pmu=0;}
   p = sqrt(E2+pmu2-pmu*ecos2)+pfermi;
   
   d1 = Ep-sqrt(p*p+ms2);
   
   Ep = Emax - emax;
   pmu2 = (emax-vc)*(emax-vc)-ml2; 
   if ( pmu2 > 0 ){
	 pmu = sqrt(pmu2);
   }else{
	 pmu2 = 0;
	 pmu  = 0;
   }

   p = sqrt(E2+pmu2-pmu*ecos2)+pfermi;
   
   d2 = Ep-sqrt(p*p+ms2);

  do{
    
    if(/*E-ml<ml && */d1*d2>0.0){//std::cout<<"pmu d1 "<<sqrt((emin-vc)*(emin-vc)-ml2) <<" pmu d2 "<<pmu<<" d1 "<<d1<<" d2 "<<d2<<" test "<<emin<<" "<<(emin+emax)/2<<" "<<emax<<" cos "<<xcos<<std::endl;
      i=99;
      emax=emin=ml;
    }    
    else if( d1*d2 < 0.0 ) {
      
      double emid = 0.5*(emax+emin);
      
      Ep = Emax - emid;

      double pmu2 = (emid-vc)*(emid-vc)-ml2;
      pmu = sqrt(pmu2);
      p = sqrt(E2+pmu2-pmu*ecos2)+pfermi;

      d3 = Ep-sqrt(p*p+ms2);

      if( d1*d3 < 0. ) {
	emax = emid;
	d2 = d3; 
      }
      else if( d2*d3 < 0. ) {
	emin = emid;
	d1=d3; 
      }

    }

    delta = (emax-emin)/(emax+emin)*2.; 
    i++;
  } while ( delta > 1.e-4 && i < 1e+2 ); 

  return emin; 
}


void N1p1h::SelectNucleusLepton(int leptonid,int nuclei) {
  static int Savedleptonid = 0; 
  static int Savednulei = 0; 
    
  // if( leptonid == Savedleptonid && Savednulei == nuclei ) return; 
  
  Savedleptonid = leptonid; 
  Savednulei = nuclei;
  
  CheckNuclei(nuclei);
  
  double ml = leptonmass[abs(leptonid)];
  double massMeV = ml*1000.;
    
  setleptonmass_(&leptonid,&massMeV);

  if( nuclei == 1 ) {
    if( leptonid < 0 ) 
      testieta_.ieta = -1;
    else 
      testieta_.ieta = 1;
    return;
  } 

  if( leptonid < 0 ) 
    Nuclei[nuclei]->CopyNucleiStructure(-1);
  else 
    Nuclei[nuclei]->CopyNucleiStructure(1);

  if( !ApplyBindEnergy  )  Nuclei[nuclei]->NoBindEnergy();

}

void N1p1h::ComputeHintegrals(int nuclei){
  
  if(nuclei !=1) std::cout<<"Entering wrong loop. This is only for Hydrogen.Requesed: "<< nuclei << std::endl;

  double hbarc = 197.3269602;
  
  int binmax = Emax/Enubin;
  //  double rmax = Nuclei[nuclei]->FindRMaximumDensity();
  int id,lid;
  double ml;

  double Enu,Enu2;
  double Xsect,Xsect2,MaxXsect;
  double maxcos,maxtl;
  int ntot;
  double q2max,emumaxemumin;
  double qval=0.;
  double lr,ang,xcos,tmommax;
  //double pmommax; 
  double Emumax,Emumin;	
  double Tlepmax,Tlepmin,tl;
  double El;
  double tmommin;
  double ldp,ldpsaved;
  double correctionangle;
  double xs,xs2,dd;
  double xsect2d;

  if( !ReadFromTable(nuclei) ) {
    
    double pi = acos(-1.);
    
    for( unsigned int il = 0; il < neutrinoIdlist.size(); il++ ) { // Looping for all neutrinos and antineutrinos.
      id = neutrinoIdlist[il];
      lid = daughter[id];
      Precompindx pindx(id,nuclei);
      std::vector<double> integral;
      std::vector<double> maximal;
      ml = leptonmass[abs(lid)];
      double ml2 = ml*ml;    
      //      binmax = 20; 

      SelectNucleusLepton(id,1); 
      
      std::cout << " Neutrino " << id << " Lepton  " << lid << " mass " << ml << std::endl; 
      
      for( int i = 0; i < binmax; i++ ) {  //looping over the neutrino energy
	Enu = (double)i*Enubin+Enubin*0.5;

	Enu2 = Enu*Enu;
	Xsect = Xsect2 = 0.;
	MaxXsect = 0.;
	maxcos   = 0.;
	maxtl    = 0.; 
	xsect2d  = 0.;
	
	ntot = 1000000; 

	q2max =0.;
	emumaxemumin = 0.; 

	if( id > 0 ) {
	  integral.push_back(0.);
	  maximal.push_back(0.);
	  continue; 
	}

	for( int  iloop = 0 ; iloop < ntot; iloop++ ){
	  //angle neutrino-lepton
	  //	  ang = Random()*pi;
	  xcos =1;// cos(ang);
	  
	  double intfactor4 = 1.; 
	  intfactor4 *= 2;;
	  
	  //lepton kinetic energy
	  Emumax=Enu;
	  	  	  
	  Tlepmax = Emumax-ml;

	  if( Tlepmax < 0.) continue;
	  
	  Emumin=ml;

	  Tlepmin = Emumin-ml;

	  if( Tlepmin > Tlepmax ) continue;

	  tl= Tlepmin+(Tlepmax-Tlepmin)*Random();
	  
	  intfactor4 *= (Tlepmax-Tlepmin);

	  // Four differential
	  
	  double xsl = DoubleDifferential(id,nuclei,Enu,tl,xcos);
	  xcos=cosH;
	  if( xsl < 0. || std::isnan(xsl)  ) xsl=0.;

	  xs =  xsl*intfactor4;
	  xs2 =  xsl*xsl*intfactor4;

	  if( xs > MaxXsect ) { MaxXsect = xs; maxcos = xcos; maxtl = tl; emumaxemumin = Emumax-Emumin;
	  }

	  Xsect   += xs;
	  Xsect2   += xs2;
	}

	Xsect   /= (double)ntot;
	Xsect2  /= (double)ntot;
	
	integral.push_back(Xsect);
	maximal.push_back(MaxXsect);
	
        if( debug ) 
	std::cout << " nuclei " << nuclei << " Lepton " << lid << " Enu " << Enu << " xsect " << Xsect << "+-" << sqrt(Xsect2-Xsect*Xsect)/sqrt(ntot) << " max " << MaxXsect << " ( " << maxcos << " , " << q2max  << " , " <<  emumaxemumin << " ) " << std::endl; 

      }
      
      
      IntCrossSection[pindx] = integral;
      MaxCrossSection[pindx] = maximal; 
      
      
    }

    RecomputeMaximum(nuclei); 

    WriteTable(nuclei);
  }

return; 
}

void N1p1h::GenerateHVectors(int id,int nuclei,double pnu[4],  double p[4][4], int idpart[4], int parent[4], double &R, double &xsec){

  double hbarc = 197.3269602;

  Precompindx pindx(id,nuclei);
  
  double Enu = pnu[0];
  
  double Enu2 = Enu*Enu; 

  double dP,vc; 

  int lid = daughter[id];

  double mass = leptonmass[abs(lid)];
  double mass2 = mass*mass; 

  double Tlepmax = Enu-mass;
  
  int ibin = (int) ((Enu+Enubin*0.001)/Enubin);

  double cmax = CrossSectionmax(id,nuclei,Enu);
    
  double cx,val; 

  int    it = 0; 
  double pi = acos(-1.);

  numberofiterationslastevent = 0; 
  double qval=0.;

  SelectNucleusLepton(lid,nuclei);  

  //  Nuclei[nuclei]->DumpInfo();
  
  int r0 = 0; 
  int r1 = 0; 
  int r2 = 0; 
  int r3 = 0;
  int r4 = 0; 
 
  do {

    double intfactor4=1;
    
    it++; 
    //    double ang = Random()*pi;
    double xcos = 1;//cos(ang); 

    intfactor4 *= 2;//sqrt(1.-xcos*xcos)*pi;

    R = 0.;
    dP=0.;
    
    double vc = 0.; 
        
    //lepton kinetic energy
    double Emumax = Enu;
    double Emumin = mass; 
    
    // Emumax = Enu; Emumin=mass; 
    
    double Tlepmax = Emumax-mass;
    double Tlepmin = Emumin-mass;

    // Keep the values for debugging 

    TlepmaxKept = Tlepmax; 
    TlepminKept = Tlepmin; 
        
    if( Tlepmax < 0. || Tlepmin > Tlepmax ) { 
      cx = 0; val = 1.;
      r0++;
      continue;
    } 
    
    double Tlep = Tlepmin+(Tlepmax-Tlepmin)*Random();

    intfactor4 *= (Tlepmax-Tlepmin);

    cx =  DoubleDifferential(id,nuclei,Enu,Tlep,xcos)*intfactor4;
    xcos=cosH;
    if( cx <= 0.0 ) { 
      cx = 0;
      val = 1.;
      r3++;
      continue; 

    }
    
    if( cx > cmax ) { 
      //      std::cout << " Enu " << Enu << " max " << cmax << " < " << cx << "( " <<  MaxCrossSection[pindx][ibin] << " , " <<  MaxCrossSection[pindx][ibin+1] << " ) " << std::endl; 
      cmax = cx; 
    }

    val = cmax*Random(); 

  } while ( cx < val && it < 50e+4 ); 


  //  std::cout << (double)r0/(double)it << "  " << 100.*(double)r1/(double)it << "  " << 100.*(double)r2/(double)it << "  " << 100.*(double)r3/(double)it << "  " << std::endl;

  xsec = cx;

  if( xsec == 0  )std::cout << it << std::endl;
  
  double z[3];

  z[0] = pnu[1]/Enu;
  z[1] = pnu[2]/Enu;
  z[2] = pnu[3]/Enu;

  comptmomentum_(z);
  
  for(int i=0;i<4;i++){
    p[0][i]=fourvectors_.vPnu[i];
    p[1][i]=fourvectors_.vPn[i];
    p[2][i]=fourvectors_.vPmu[i];
    p[3][i]=fourvectors_.vPp[i];
  }

  parent[0] = 0;
  parent[1] = 0;
  parent[2] = 1;
  parent[3] = 2;

  idpart[0] = id;
  if( id < 0 ) {
    idpart[1] = idproton;
    idpart[3] = idneutron; 
  }
  else {
    idpart[3] = idproton;
    idpart[1] = idneutron; 
  }
  idpart[2] = daughter[id]; 
  
  numberofiterationslastevent = it; 

  return; 
}
