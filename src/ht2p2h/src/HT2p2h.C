#include "HT2p2h.h"
#include <dirent.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fstream>

#define PNOFFSET 10000

HT2p2h::HT2p2h(std::string dirname, bool ApplyBind, bool d )
{
  Initialize(dirname, ApplyBind, d);
}

HT2p2h::HT2p2h()
{
}

void
HT2p2h::Initialize(std::string dirname, bool ApplyBind, bool d ){

  Enubin = 0.01;  // 10 MeV 
  Emax   = 30.;   // 30 GeV

  RFG = false; 

  HadKinMode = 0; 
  
  HFSparam1 = HFSparam2 = 0; 

  ApplyBindEnergy = ApplyBind; 
  
  debug = d; 

  std::cout << " RAND_MAX =   " << RAND_MAX << std::endl; 

  default_pn_nn_fraction = 0.9; 
  
  // Parameters 
  
  GFermi = 1.1664e-5;   // Gfermi is in GeV^(-2))
  facconv = 0.19733*0.19733*1.e+15;

  leptonmass[12] = 0.000510998910;  // lepton mass in GeV
  leptonmass[14] = 0.105658357; // muon mass in GeV
  leptonmass[16] = 1.77682; // tau mass in GeV 

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

  std::cout << " HT2p2h Reading files " << std::endl; 
  
  ReadTensors(dirname); 
  
  return; 
}


void HT2p2h::ReadTensors(std::string dirname ){
  // Read Tensor from a directory 
  DIR *dir;
  dirent *ent;
  
  if ( ( dir = opendir (dirname.c_str()) ) != 0 ) {

    std::cout << " HT2p2h reading directory " << dirname << std::endl; 
	
    while( ( ent = readdir(dir) ) != 0 ) {
	  
      if( ent->d_name[0] == '.' ) continue; 
	  
      HadronTensor *ht = new HadronTensor();
	  ht->Initialize(dirname+"/"+std::string(ent->d_name));
      int nuclei = ht->GetNuclei();
      if( nuclei == -1 ) continue;  // This is not a proper file. 
      if( ht->IsPN()  ){
		Tensor[nuclei+PNOFFSET] = ht;
		Tensor_fname[nuclei+PNOFFSET] = dirname+"/"+std::string(ent->d_name);
		Tensor_init[nuclei+PNOFFSET] = false;
		
		if( nuclei == 28 ) { // Create a copy for 27 (Al) 
		  Tensor[27+PNOFFSET] = ht;
		  Tensor_fname[nuclei+PNOFFSET] = dirname+"/"+std::string(ent->d_name);
		  Tensor_init[27+PNOFFSET] = false;
		}
      }
      else{
		Tensor[nuclei] = ht;
		Tensor_init[nuclei] = false;
		if( nuclei == 28 ) { // Create a copy for 27 (Al) 
		  Tensor[27] = ht;
		  Tensor_init[27] = false;
		}

      }
	  
      std::cout << " Tensor for nuclei " << nuclei << " has been read " << std::endl; 
	  
      /*
		InitializeNucleus(nuclei); 
	
		Precompindx pindx1(1,nuclei);
		
		Precompindx pindx2(-1,nuclei);
		
		if( !ht->IsPN() )   // Only for the total. 
		ComputeIntegrals(nuclei);
      */
      
    }
  }
  else {
    std::cout << " Directory " << dirname << " not found " << std::endl; 
    exit(0); 
  }
}


int HT2p2h::ReadIntegrals(int nuclei){
  
  std::ifstream infile;

  std::string fname;
  fname = Tensor_fname[nuclei+PNOFFSET];

#define LOCALBUFSIZ (1024*1024*10)
  char *line = NULL;
  char *item = NULL;
  int  ipar, ret, i;

  line = malloc(LOCALBUFSIZ);
  item = malloc(LOCALBUFSIZ);
  infile.open (fname.c_str());
  if (infile.fail()){
	std::cout << "Failed to open Tensor file" 
			  << fname 
			  << std::endl;
	exit(1);
  }

  std::cout << " Reading Cross-section from the Tensor File " 
			<< fname << std::endl; 

  int crstbl_exists = 0;
  /* search for the Cross-sections block  */
  while(infile.getline(line,LOCALBUFSIZ)){
	if (strstr (line,"#Cross-section Table")){
	  crstbl_exists = 1;
	  break;
	}
  }
  if (crstbl_exists == 0){
	free(line);
	free(item);
	return -1;
  }

  /* search for the parameters block  */
  int param_agreed = 0;
  while(infile.getline(line,LOCALBUFSIZ)){
	if ( strstr(line,"#PARAMSTART") ){
	  break;
	}
  }
  /* search for the table with same parameter value */
  while(infile.getline(line,LOCALBUFSIZ)){  
	if (line[0] == 0) continue;
	if (line[0] != '#') break;
	if (strncmp (line,"#NV2P2HQVAL",11)==0){
	  sscanf(line,"%s %d",item, &ipar);
	  if (ipar == nieves2p2hpar_.nv2p2hqval){
		param_agreed = 1;
		break;
	  }
	}
  }  
  if (param_agreed == 0){
	free(line);
	free(item);
	return -1;
  }

  int    nuclei_file;
  int    ebinmax;
  double enubin;

  /* read # of nuclei */
  if (infile.getline(line,LOCALBUFSIZ)){
	nuclei_file = atoi(line);
  }else{
	free(line);
	free(item);
	return -1;
  }

  if (nuclei_file != nuclei){
	std::cout << "Inconsistent nuclei was stored" << std::endl;
	free(line);
	free(item);
	return -2;
  }

  /* read # of energy bins */
  if (infile.getline(line,LOCALBUFSIZ)){
	ebinmax = atoi(line);
  }else{
	free(line);
	free(item);
	return -3;
  }

  /* read energy bin size */
  if (infile.getline(line,LOCALBUFSIZ)){
	enubin = atof(line);
  }else{
	free(line);
	free(item);
	return -4;
  }

  std::vector<double> integrals[neutrinoIdlist.size()];
  std::vector<double> maximals[neutrinoIdlist.size()];
  double crs[neutrinoIdlist.size()];
  double max[neutrinoIdlist.size()];

  double enu;

  for (i = 0 ; i < ebinmax ; i++){
	infile.getline(line,LOCALBUFSIZ);
	ret = sscanf(line,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
				 &ipar, &enu, 
				 &crs[0], &max[0], &crs[1], &max[1], &crs[2], &max[2],
				 &crs[3], &max[3], &crs[4], &max[4], &crs[5], &max[5])
	if (ret != 14){
	  free(line);
	  free(item);
	  return -2;
	}
	if ( enu - ((double)i*enubin+enubin*0.5) > enubin/100. ){
	  free(line);
	  free(item);
	  return -3;
	}

	for( unsigned int il = 0; il < 6;  il++ ) {	
	  integrals[il].push_back(crs[il]);
	  maximals[il].push_back(max[il]);
	}
  }

  for( unsigned int il = 0; il < 6;  il++ ) {	
	int id = neutrinoIdlist[il]; 
	Precompindx pindx(id,nuclei);
    IntCrossSection[pindx] = integrals[il];
    MaxCrossSection[pindx] = maximals[il]; 
  }

  free(line);
  free(item);

  return 0;

}


void HT2p2h::ComputeIntegrals(int nuclei){

  int ret;
  CheckNuclei_2(nuclei);  
  
  int binmax = Emax/Enubin;
  
  //  std::cout << " Looping for " << neutrinoIdlist.size() << " neutrinos with " << binmax << " bins " <<  std::endl; 

  ret = ReadIntegrals(nuclei);
  if (ret == 0){
	/* cross-section was stored in the tensor file */
	return;
  }


  for( unsigned int il = 0; il < neutrinoIdlist.size(); il++ ) {
    int id = neutrinoIdlist[il]; 
    Precompindx pindx(id,nuclei);
 
    Precompindx qindx(id<0? -1:1,nuclei); 

    std::vector<double> integral;
    std::vector<double> maximal;

    HadronTensor *tlc;
    tlc = Tensor[nuclei];
    if ( tlc == NULL )
      tlc = Tensor[PNOFFSET+nuclei];

    for( int i = 0; i < binmax; i++ ) {
      double ml = leptonmass[abs(id)];
      double Enu = (double)i*Enubin+Enubin*0.5;
      double Xsect = 0.;
      double MaxXsect = -1000.;

      double stepQ0stepQ3 = tlc->GetQ0Step()*tlc->GetQ3Step();
      
      for( double q0 = 0. ; q0 < tlc->GetMaxQ0(); q0 += tlc->GetQ0Step()  ) {
	double Q0 = q0+tlc->GetQ0Step()/2.;

	//	std::cout << " Enu " << Enu << " Q0  " << Q0 << "  " << tlc->GetQ0Step() << "  " << tlc->GetMaxQ0() << std::endl; 
	if( Enu-Q0 < 0. ) continue; 
	double El = Enu-Q0;
	if( El < ml ) continue; 
	double pl = sqrt(El*El-ml*ml); 
	double tl = El-ml;
	double Enu2 = Enu*Enu;
	double pl2 = pl*pl;
	double plEnu = pl*Enu;
	
	for( double q3 = q0; q3 < tlc->GetMaxQ3(); q3 += tlc->GetQ3Step()  ) {
	  double Q3 = q3+tlc->GetQ3Step()/2.;
	  double xcos = (Q3*Q3-Enu2-pl2)/(-2.*plEnu);
	  if ( xcos > 1. || xcos < -1. ) continue; 
	  
	  double xs = 0.;

	  for(  int i = 0; i < 20; i++ ) {
	    double R = GenerateR(nuclei,-1);
	    double lxs = DoubleDifferential(id,nuclei,Enu,tl,xcos,R);
	    if( lxs > 0. ) 
	      xs += lxs/20.;
	    if( lxs > MaxXsect ) MaxXsect = lxs*1.05; 
	  }

	  double Jacobian = Q3; 
	  
	  if( isnan(xs) ) { std::cout << " xs isnan " << xs << "  " << El << "  " << xcos << "  " << Q0 << "  " << Q3 << std::endl; continue; }

	  if( xs < 0. ) break; // IF this condition is fulfilled, next q3 will also, so we break. 
	  
	  
	  double alocal = xs*Jacobian/plEnu;
	  Xsect += alocal;
	}
      }
     
      integral.push_back(Xsect*stepQ0stepQ3); 
      maximal.push_back(MaxXsect);
      if( debug ) 
	std::cout << " nuclei " << nuclei << " Lepton " << id << " Enu " << Enu << " xsect " << Xsect << " max " << MaxXsect << std::endl; 
    }

    IntCrossSection[pindx] = integral;
    MaxCrossSection[pindx] = maximal; 
    
  }
  
  std::vector<double> integral_dump[neutrinoIdlist.size()];
  std::vector<double> maximal_dump[neutrinoIdlist.size()];

  for( unsigned int il = 0; il < neutrinoIdlist.size(); il++ ) {
    int id = neutrinoIdlist[il]; 
    Precompindx pindx(id,nuclei);
	
    integral_dump[il] = IntCrossSection[pindx];
    maximal_dump[il]  = MaxCrossSection[pindx];

  }  

  std::string fname;
  fname = Tensor_fname[nuclei+PNOFFSET];

  std::ofstream offile(fname.c_str(),std::ofstream::app);
  if (offile.fail()){
	return;
  }

  offile << "#Cross-section Table" << std::endl;
  offile << "#PARAMSTART" << std::endl;
  offile << "#NV2P2HQVAL   " << nieves2p2hpar_.nv2p2hqval << std::endl;
  offile << nuclei << std::endl;
  offile << binmax << std::endl;
  offile << Enubin << std::endl;

  for( int i = 0; i < binmax; i++ ) {
	offile << i << " ";
	offile << (double)i*Enubin+Enubin*0.5 << " ";
	for( unsigned int il = 0; il < neutrinoIdlist.size(); il++ ) {
	  offile << integral_dump[il][i] << " "
			 << maximal_dump[il][i] << " ";
	}
	offile << std::endl;
  }
  
  return; 
}


double HT2p2h::IntegralCrossSection(int id,int nuclei,double Enu){

  CheckNuclei(nuclei);  
  
  int ibin = (int) ((Enu+Enubin*0.001)/Enubin); 

  Precompindx pindx(id,nuclei);  

  if( Enu > Emax ) {
    int ibinmax = IntCrossSection[pindx].size();
    return IntCrossSection[pindx][ibinmax-1];
  }
  
  double a0 = IntCrossSection[pindx][ibin]; 
  double a1 = IntCrossSection[pindx][ibin+1];

  double cross_section = (a1-a0)/Enubin*(Enu-Enubin*(double)ibin)+a0;  

  return cross_section; 
}


int HT2p2h::GenerateLeptonKinematics(int id,int nuclei,double Enu,double &TLepton,double &xcos, double &R ) {
 
  CheckNuclei(nuclei);  

  Precompindx pindx(id,nuclei);
  double Tlepmax = Enu;

  double qval = GetQ0Fermivalue(nuclei,id,R); 

  if( ApplyBindEnergy )  
    Tlepmax -= (GetBind(nuclei,id)+qval);
  
  int ibin = (int) ((Enu+Enubin*0.001)/Enubin);

  double max;
  
  if( ibin >= MaxCrossSection[pindx].size()-1 ) {
    ibin =  MaxCrossSection[pindx].size()-2;
    max =  (MaxCrossSection[pindx][ibin+1]-MaxCrossSection[pindx][ibin])*(Enu-Enubin*(double)ibin)/Enubin+MaxCrossSection[pindx][ibin] ;
  }
  else 
    max =  (MaxCrossSection[pindx][ibin+1]-MaxCrossSection[pindx][ibin])*(Enu-Enubin*(double)ibin)/Enubin+MaxCrossSection[pindx][ibin] ;

  double tl,cos,xs;

  cos = 2.*Random()-1.;
  
  tl = Random()*Tlepmax;
  
  xs = DoubleDifferential(id,nuclei,Enu,tl,cos, R ); 
  
  if( xs == 0 ) return -1;
  
  double y = max * Random();
  
  if( debug ) 
    std::cout <<  "  y =  " << y << " max =  "<< max << " xs  " << xs  <<  std::endl; 
  
  if( xs > max ) {std::cout << " Cross-Section " << xs << " larger than maximum " << max << " cos " << cos << std::endl; MaxCrossSection[pindx][ibin] = max*1.1; }
  
  if( xs < y ) return -1;
  

  TLepton = tl; 
  xcos = cos; 

  return 1; 
}


int HT2p2h::GenerateHadronKinematics( int id, int nuclei, double frac, double q[4], double ni1[4], double ni2[4], double n1[4], double n2[4] , int idNucleon[4], double &R )  {

  CheckNuclei(nuclei);  

  double FermiLevel;
  
  double m1,m2;
  double mf1,mf2; 
  double q0;
  
  Precompindx pindx(id,nuclei);

  // We have to use locally the corrected q0 because the function might be called recursively and it substract values every time.
  if( ApplyBindEnergy ) 
    q0 = q[0] - GetBind(nuclei,id);       // Take the bind energy out of the equation.
  else
    q0 = q[0];
   
  if( id > 0 ) {   // Neutrinos n --> p 
    if ( Random() > frac ) {
      m1 = neutronmass;  // neutron 
      mf2 = m2 = neutronmass;  // neutron
      mf1 = protonmass; // proton
      idNucleon[0] = idneutron;
      idNucleon[1] = idNucleon[3] = idneutron;
      idNucleon[2] = idproton;
      if( RFG ) 
	FermiLevel = GetFermiRFG(nuclei);
      else 
	FermiLevel = GetFermiLFG(R,nuclei,-1);
    } else {
      m1 = neutronmass; // neutron
      mf2 = m2 = protonmass; // proton
      mf1 = protonmass; // proton
      idNucleon[0] = idneutron;
      idNucleon[1] = idNucleon[3] = idproton;
      idNucleon[2] = idproton;
      if( RFG ) 
	FermiLevel = GetFermiRFG(nuclei);
      else 
	FermiLevel = GetFermiLFG(R,nuclei,1);
    }
  }
  else {    // Anti-Neutrinos p --> n 
    if ( Random() > frac ) {
      m1 = protonmass; // proton
      mf2 = m2 = protonmass; // proton
      mf1 = neutronmass; // neutron
      idNucleon[0] = idproton;
      idNucleon[1] = idNucleon[3] = idproton;
      idNucleon[2] = idneutron;
      if( RFG ) 
	FermiLevel = GetFermiRFG(nuclei);
      else 
	FermiLevel = GetFermiLFG(R,nuclei,1);
    } else {
      m1 = protonmass; // proton
      mf2 = m2 = neutronmass; // neutron
      mf1 = neutronmass; // neutron
      idNucleon[0] = idproton;
      idNucleon[1] = idNucleon[3] = idneutron;
      idNucleon[2] = idneutron;
      if( RFG ) 
	FermiLevel = GetFermiRFG(nuclei);
      else 
	FermiLevel = GetFermiLFG(R,nuclei,-1);
    }
  }

  int ntries = 0; 

  double qm = sqrt( q[1]*q[1]+q[2]*q[2]+q[3]*q[3] ); 

  int code = 0; 

  double pmin = 0.;

  double pminsaved = 0.; 

  // There is not a possible solution.

  double Efmin1 = sqrt(FermiLevel*FermiLevel+mf1*mf1);
  double Efmin2 = sqrt(FermiLevel*FermiLevel+mf2*mf2);
   
  do{
    
    ntries++;

    //    if( ntries == 1e+4 ) return GenerateHadronKinematics(id,nuclei,frac, q, ni1, ni2, n1, n2 , idNucleon,R); // If error try a new radial position. 

    if( ntries == 1e+4 ) return -1;
    
    // q0+eN2+eN1>Efmin1+Efmin2  Minimum energy to produce two nucleons above the Fermi Level
    // eN2 > Efmin1+Efmin2-q0-eN1 and get eN1 = Efmin1

    if( (Efmin2-q0)*(Efmin2-q0)-m2*m2 < 0. )
      pmin = 0.; 
    else 
      pmin = sqrt((Efmin2-q0)*(Efmin2-q0)-m2*m2); 

    if( pmin > FermiLevel ) return -6; // we will never find a solution for this case. 
    
    double momN2 = pow(Random()*(pow(FermiLevel,3.)-pow(pmin,3.))+pow(pmin,3.),1./3.);
    double eN2   = sqrt( momN2*momN2 + m2*m2 );
    double cosN2 = 2.*Random()-1.;
    double sinN2 = sqrt(1.-cosN2*cosN2);
    double phiN2 = 2.*pi*Random();

    // q0+eN2+eN1>Efmin1+Efmin2  Minimum energy to produce two nucleons above the Fermi Level
    // eN1 > Efmin1+Efmin2-q0-eN2

    if( (Efmin1+Efmin2-q0-eN2)*(Efmin1+Efmin2-q0-eN2)-m1*m1 < 0. )
      pmin = 0.; 
    else 
      pmin = sqrt((Efmin1+Efmin2-q0-eN2)*(Efmin1+Efmin2-q0-eN2)-m1*m1); 

    if( pmin >= FermiLevel ) { code = -1; continue; }

    ni2[0] = eN2;
    ni2[1] = momN2*sinN2*cos(phiN2);
    ni2[2] = momN2*sinN2*sin(phiN2);
    ni2[3] = momN2*cosN2; 
    
    if( HadKinMode == 0  ) {
      double momN1 = pow(Random()*(pow(FermiLevel,3.)-pow(pmin,3.))+pow(pmin,3.),1./3.);
      double eN1   = sqrt(momN1*momN1 + m1*m1);   
      double cosN1 = 2.*Random()-1.;
      double sinN1 = sqrt(1.-cosN1*cosN1); 
      double phiN1 = 2.*pi*Random();
      ni1[0] = eN1;
      ni1[1] = momN1*sinN1*cos(phiN1); 
      ni1[2] = momN1*sinN1*sin(phiN1);
      ni1[3] = momN1*cosN1; 
    }
    else if( HadKinMode == 1  ){  // Back to Back 
      ni1[0] = ni2[0];
      for(int ii = 1; ii < 4;ii++ )
	ni1[ii] = -ni2[ii];
    }
    else if( HadKinMode == 2  ){  // In the same direction
      for(int ii = 0; ii < 4;ii++ )
	ni1[ii] = ni2[ii];
    }
    

    double qf[4];

    for( int i = 0; i < 4 ; i++ ) qf[i] = q[i] + ni1[i] + ni2[i];

    double qval = GetQ0Fermivalue(nuclei,id,R); 

    if( ApplyBindEnergy ) 
      qf[0] -= (GetBind(nuclei,id)+qval); 

    double qfm2 = qf[1]*qf[1] + qf[2]*qf[2] + qf[3]*qf[3];
    double qfm  = sqrt(qfm2);

    if( qfm == 0.0 ) {  code = -6; continue; }
    
    double c2 = qf[0]*qf[0] - qfm2;  // Mass^2 of the final hadron system == final state rest frame energy.
    
    if( c2 < 0.0 ) {  code = -2; continue; }

    double c = sqrt(c2); // Mass of the final hadronic system 

    if( mf1+mf2 > c ) { code = -3; continue; }  // Not enough energy in the rest frame. 
 
    double a2 = mf1*mf1;
    double b2 = mf2*mf2;

    // Working in the q rest frame. Find the momentum of the outgoing nucleons in the RF.
    // c = sqrt(pRF*pRF+m1*m1)+sqrt(pRF*pRF+m2*m2)
    // c - sqrt(pRF*pRF+m1*m1) = sqrt( pRF*pRF+m2*m2) 
    // c2 + pRF*pRF+m1*m1 - 2*c*sqrt(pRF*pRF+m1*m1) = pRF*pRF+m2*m2
    // c2 + m1*m1 -m2*m2 = 2*c*sqrt(pRF*pRF+m1*m1)
    // ( c2 + a2 - b2 )^2 / 4 c2 - m1*m1 = pRF*pRF 
    
    double pRF = sqrt( (c2+a2-b2)*(c2+a2-b2)/(4.*c2) - a2); 
        
    double e1RF = sqrt(pRF*pRF+a2);
    double e2RF = sqrt(pRF*pRF+b2);
    
    // Boost 

    double gamma = qf[0]/c;
    double beta  = sqrt(1.-1./(gamma*gamma));
     
    double cosangRF;

    if( HFSparam1 != 0 || HFSparam2 != 0 ){
      double aaaa;

      do {     
	do {
	  cosangRF = 2.*Random()-1.;
	  aaaa = 1.+4.*cosangRF*HFSparam1;
	} while( aaaa <= 0. ); 
	
	cosangRF = (-1. + sqrt(aaaa) )/2./HFSparam1;
      } while( cosangRF > 1. || cosangRF < -1. );
	
    } else
      cosangRF = 2.*Random()-1.;

    double sinangRF = sqrt(1.-cosangRF*cosangRF);
    
    double p1L = gamma*(pRF*cosangRF+beta*e1RF);
    double p2L = gamma*(-pRF*cosangRF+beta*e2RF);

    double p1 = sqrt(p1L*p1L+pRF*sinangRF*pRF*sinangRF); 
    double p2 = sqrt(p2L*p2L+pRF*sinangRF*pRF*sinangRF);
    
    if( isnan(p1) ) std::cout << " ERROR " << p1 << " " << gamma <<  "  " << pRF << "   " << beta << "   " << e1RF << std::endl;
    if( isnan(p2) ) std::cout << " ERROR " << p2 << " " << gamma <<  "  " << pRF << "   " << beta << "   " << e1RF << std::endl; 
    
    // Pauli blocking. Outgoing nucleons below fermi level. 

    if( p1 < FermiLevel || p2 < FermiLevel ) { code = -4;continue; }

    double tx[3],ty[3];

    if( HadKinMode == 0  ) {
      tx[0] = (qf[2]*q[3]-q[2]*qf[3]);
      tx[1] = (qf[3]*q[1]-q[3]*qf[1]);
      tx[2] = (qf[1]*q[2]-q[1]*qf[2]);
    }
    else { // In this case, qf == q. We take one of the initial nucleus to define the reference system. 
      tx[0] = (ni1[2]*qf[3]-qf[2]*ni1[3]);
      tx[1] = (ni1[3]*qf[1]-qf[3]*ni1[1]);
      tx[2] = (ni1[1]*qf[2]-qf[1]*ni1[2]);
    }
   
    double txm = sqrt(tx[0]*tx[0]+tx[1]*tx[1]+tx[2]*tx[2]);
    
    tx[0] /= txm; 
    tx[1] /= txm; 
    tx[2] /= txm; 

    ty[0] = (tx[1]*qf[3]-qf[2]*tx[2])/qfm;
    ty[1] = (tx[2]*qf[1]-qf[3]*tx[0])/qfm;
    ty[2] = (tx[0]*qf[2]-qf[1]*tx[1])/qfm;

    double phi = 2.*pi*Random(); 

    n1[0] = sqrt(p1*p1+mf1*mf1); 
    n1[1] = p1L*qf[1]/qfm+pRF*sinangRF*cos(phi)*tx[0]+pRF*sinangRF*sin(phi)*ty[0];
    n1[2] = p1L*qf[2]/qfm+pRF*sinangRF*cos(phi)*tx[1]+pRF*sinangRF*sin(phi)*ty[1];
    n1[3] = p1L*qf[3]/qfm+pRF*sinangRF*cos(phi)*tx[2]+pRF*sinangRF*sin(phi)*ty[2];
    
    n2[0] = sqrt(p2*p2+mf2*mf2); 
    n2[1] = p2L*qf[1]/qfm-pRF*sinangRF*cos(phi)*tx[0]-pRF*sinangRF*sin(phi)*ty[0];
    n2[2] = p2L*qf[2]/qfm-pRF*sinangRF*cos(phi)*tx[1]-pRF*sinangRF*sin(phi)*ty[1];
    n2[3] = p2L*qf[3]/qfm-pRF*sinangRF*cos(phi)*tx[2]-pRF*sinangRF*sin(phi)*ty[2];
    
    break;  // Good solution found. 

  } while(1); 

  
  return ntries; 
}





double   HT2p2h::DoubleDifferential(int id,int nuclei,double Enu,double TLep,double xcos, double R, bool pn ) {
 
  CheckNuclei_2(nuclei);  
  
  Precompindx pindx(id,nuclei);  

  Precompindx qindx(id<0?-1:1,nuclei);
  
  double xcos2  = xcos*xcos;
  double xsin2  = 1.-xcos2;
  double xsin   = sqrt(xsin2);
  double anghalf= acos(xcos)/2.;
  double xcosh  = cos(anghalf);
  double xcosh2 = xcosh*xcosh;
  double xsinh  = sin(anghalf);
  double xsinh2 = xsinh*xsinh;

  double xmlepton = leptonmass[abs(id)];
  double xmlepton2 = xmlepton*xmlepton;
  
  double energylepton   = TLep+xmlepton;         // final lepton energy       
  double energyneutrino = Enu;

  //     Four-momentum of the incoming neutrino

  double xk[4];
  
  xk[0]=energyneutrino;
  xk[1]=0.0;
  xk[2]=0.0;
  xk[3]=energyneutrino;

  double q0=energyneutrino-energylepton;
  double q0nucleus=q0;

  double qval = GetQ0Fermivalue(nuclei,id,R); 

  if( ApplyBindEnergy )
    q0nucleus-= (GetBind(nuclei,id)+qval);
  
  double q02 = q0*q0; 

  double xkprima[4];
  double q[4];
  double diffCS; 
  
  if( q0nucleus > 0.0 ){
    xkprima[0]=energylepton;
    double xmodkprima=sqrt(energylepton*energylepton-xmlepton2);
    xkprima[1]=xmodkprima*xsin;
    xkprima[2]=0.0;
    xkprima[3]=xmodkprima*xcos;

    q[0] = q0nucleus;
    q[1] = xk[1]-xkprima[1];
    q[2] = xk[2]-xkprima[2];
    q[3] = xk[3]-xkprima[3];
    
    double qm2 = q[1]*q[1]+q[2]*q[2]+q[3]*q[3];
    double qm = sqrt(qm2); 

    HadronTensor *tlc = NULL;

    if( !pn ) 
	tlc = Tensor[nuclei];
    if( tlc == NULL ){
      tlc = Tensor[PNOFFSET+nuclei];
      if( tlc == NULL ) return -1; 
    }
	  
    double T00 = tlc->Interpolate(0,0,qm,q[0]);    
    double T11 = tlc->Interpolate(1,1,qm,q[0]);
    double T33 = tlc->Interpolate(3,3,qm,q[0]);
    double T12 = tlc->Interpolate(1,2,qm,q[0]);
    double T03 = tlc->Interpolate(0,3,qm,q[0]);

    double w1= T11/2.;
    double w2= (T00+T11+q02/qm2*(T33-T11)-2.*T03*q0/qm)/2.;
    double w3= T12/qm;
    double w4= (T33-T11)/(2.*qm2);
    double w5= (T03-q0/qm*(T33-T11))/qm;
      
    if( id < 0)  w3 = -w3;

    if( debug ) {
      std::cout << " w1 " << w1 << std::endl; 
      std::cout << " w2 " << w2 << std::endl; 
      std::cout << " w3 " << w3 << std::endl; 
      std::cout << " w4 " << w4 << std::endl; 
      std::cout << " w5 " << w5 << std::endl; 
    }    

    double aux1= w1*xcos
               - w2/2.*xcos  
               + w3/2.*(energylepton+xmodkprima-(energylepton+energyneutrino)*xcos)
               + w4/2.*(xmlepton2*xcos+2.*energylepton*(energylepton+xmodkprima)*xsinh2)
               - w5*(energylepton+xmodkprima)/2.;


    if( debug ) {
      std::cout << "energylepton "<< energylepton <<std::endl; 

      std::cout << "xmodkprima " << xmodkprima << std::endl; 
      std::cout << "xsinh " << xsinh << std::endl; 

      std::cout << " aux1 " << aux1 << std::endl; 
      std::cout << " aux1 = " << w1*xcos << " - " << w2/2.*xcos << " +  "  << w3/2.*(energylepton+xmodkprima-(energylepton+energyneutrino)*xcos) << " + " << w4/2.*(xmlepton2*xcos+2.*energylepton*(energylepton+xmodkprima)*xsinh2) << " - " << w5*(energylepton+xmodkprima)/2. << std::endl; 
      
    }

    double aux0= 2.*w1*xsinh2
               +    w2*xcosh2
               -    w3*(energylepton+energyneutrino)*xsinh2+xmlepton2/(energylepton*(energylepton+xmodkprima))*aux1;


    if( debug ) std::cout << " aux0 " << aux0 << std::endl; 

    diffCS=xmodkprima*energylepton*GFermi*GFermi*2./pi*aux0*facconv;  // diffCS=dsigma/( dcos(theta) d TLep)  in 10^{-41} cm^2/GeV

    if( debug ) {
      std::cout << " diffCS " << diffCS << std::endl; 
      std::cout << " xmodkprima " << xmodkprima << std::endl; 
      std::cout << " energylepton " << energylepton << std::endl; 
      std::cout << " facconv      " << facconv  << std::endl; 
      std::cout << " GFermi      " << GFermi  << std::endl; 
    }

    diffCS /= 1000.; // in 10^{-38} cm^2/GeV

  }
  else
    diffCS=0.0;
  
  return diffCS; 
}


int HT2p2h::GenerateVectors(int id, int nuclei, double pnu[4],  double p[6][4], int idpart[6], int parent[6], double &R ){

  CheckNuclei(nuclei);  
  
  double Enu = pnu[0];
  
  double TLepton,coslepton; 
  
  double pl[4],q[4],qi1[4],qi2[4],qf1[4],qf2[4];

  double x[3],y[3],z[3];

  z[0] = pnu[1]/Enu;
  z[1] = pnu[2]/Enu;
  z[2] = pnu[3]/Enu;

  x[0] = 1.; x[1]=0.; x[2] = 0.;

  y[0] = z[1]*x[2]-x[1]*z[2];
  y[1] = z[2]*x[0]-x[2]*z[0];
  y[2] = z[0]*x[1]-x[0]*z[1];

  double yn = sqrt(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]);
  y[0] /= yn;y[1] /= yn;y[2] /= yn;
  
  x[0] = y[1]*z[2]-z[1]*y[2];
  x[1] = y[2]*z[0]-z[2]*y[0];
  x[2] = y[0]*z[1]-z[0]*y[1];
  double xn = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  x[0] /= xn;x[1] /= xn;x[2] /= xn;

  // pnu[0] = pnu[3] = Enu;
  // pnu[1] = pnu[2] = 0; 
  
  int idNucleon[4];
  double sq0;
  double mass;
  double elepton;
  double plepton;
  double phiang;
  double sinlepton;
  double cosphi;
  double sinphi;
  double fractionpn;

  int ntry;
  int hadrontries;

  while(1){

    double R = GenerateR(nuclei,-1);

    if( GenerateLeptonKinematics(id,nuclei,Enu,TLepton,coslepton,R) < 0 ) continue;
   
    mass = leptonmass[abs(id)];
    
    elepton = TLepton+mass; 
    
    plepton = sqrt( elepton*elepton - mass*mass  ); 
    
    phiang = 2.*3.141592*Random();
    
    sinlepton = sqrt(1.-coslepton*coslepton);
    
    cosphi = cos(phiang);
    sinphi = sin(phiang);
    
    pl[0] = elepton; 
    pl[1] = plepton*(coslepton*z[0]+sinlepton*(x[0]*cosphi+y[0]*sinphi));
    pl[2] = plepton*(coslepton*z[1]+sinlepton*(x[1]*cosphi+y[1]*sinphi));
    pl[3] = plepton*(coslepton*z[2]+sinlepton*(x[2]*cosphi+y[2]*sinphi));
    
    for(int i = 0; i < 4; i++ ) q[i] = pnu[i]-pl[i];
    sq0 = q[0];    
    
    //  std::cout << " Generating the hadron kinematics " << std::endl; 
    
    fractionpn = GetFraction(id,nuclei,Enu,TLepton,coslepton,R);

    while(1){
      hadrontries = GenerateHadronKinematics(id,nuclei,fractionpn,q,qi1,qi2,qf1,qf2,idNucleon,R);

      if ( hadrontries >= 0 ){
	break;
      }else{
	ntry++;
	if ( ntry > 10 ){
	  break;
	}
      }
      std::cout << "HT2p2h: Retry different R" << std::endl;
    }

    if ( hadrontries >= 0 ){
      break;
    }
    //    std::cout << "HT2p2h: Retry different lepton kinematics" << std::endl;
  }

  for(int i = 0; i < 4; i++ ){
    p[0][i] = pnu[i];
    p[1][i] = qi1[i];
    p[2][i] = qi2[i];
    p[3][i] = pl[i];
    p[4][i] = qf1[i];
    p[5][i] = qf2[i];
  }


  double pavx = 0;
  double pavy = 0;
  double pavz = 0;

			
  parent[0] = parent[1] = parent[2] = 0; // Original particles. 
  parent[3] = 1;
  parent[4] = 2;
  parent[5] = 3; 


  idpart[0] = id;
  idpart[1] = idNucleon[0];
  idpart[2] = idNucleon[1];
  idpart[3] = daughter[id];
  idpart[4] = idNucleon[2];
  idpart[5] = idNucleon[3];

  return hadrontries;
}


double HT2p2h::GetFraction(int id, int nuclei,double Enu,double TLepton,double coslepton, double R){

  CheckNuclei(nuclei);  
  
  double fraction = default_pn_nn_fraction;
  
  double totalCrossSection = DoubleDifferential(id,nuclei,Enu,TLepton,coslepton,R);
  
  double pnCrossSection = DoubleDifferential(id,nuclei,Enu,TLepton,coslepton,R,true);
  
  if( pnCrossSection >= 0. )     fraction =  pnCrossSection/totalCrossSection;
  
  //  std::cout << " FRACTION " << fraction << "  " << pnCrossSection << "  " <<  totalCrossSection <<  std::endl; 
  
  return fraction;
}

void HT2p2h::CheckNuclei(int nuclei){

  int offset = 0;
  if (Tensor[nuclei] == NULL){
	offset = PNOFFSET;
  }
  if (Tensor[nuclei+offset] == NULL){
    std::cout << " HT2p2h Error: nuclei " << nuclei << " not available " << std::endl;
    exit(1);
  }

  if( Tensor_init[nuclei+offset] == false ) {
	InitializeNucleus(nuclei); 
	Precompindx pindx1(1,nuclei);
	Precompindx pindx2(-1,nuclei);
	if(!Tensor[nuclei+offset]->IsPN() )   // Only for the total. 
	  ComputeIntegrals(nuclei);
	Tensor_init[nuclei+offset]=true;
  }
}

void HT2p2h::CheckNuclei_2(int nuclei){

  if( Tensor_init[nuclei] == false ) {
	/*
    std::cout << " HT2p2h Error: nuclei " << nuclei << " not available " << std::endl;
    exit(1);
	*/
	int offset = 0;
	if (Tensor[nuclei+PNOFFSET] != NULL){
	  if (Tensor[nuclei+PNOFFSET]->IsPN()){
		offset = PNOFFSET;
	  }
	}else{
	  std::cout << " HT2p2h Error: nuclei " << nuclei << " not available " << std::endl;
	  exit(1);
	}
  }
}


double HT2p2h::GetQ0Fermivalue(int nuclei,int id,double R){
  
  double pfermi = GetFermiLFG(R,nuclei,-1);
  
  double qval = pfermi*pfermi/2.*(1./protonmass+1./neutronmass);

  if ((nieves2p2hpar_.nv2p2hqval != 1) &&
	  (nieves2p2hpar_.nv2p2hqval != 2)){
	std::cout << "HT2p2h Error : Config param. nv2p2hqval is not properly set"
			  << "nv2p2hqval is " << nieves2p2hpar_.nv2p2hqval
			  << std::endl;
	exit(1);
  }

  if (nieves2p2hpar_.nv2p2hqval == 1){
	qval = 0.;
  }else

  return qval; 
}

