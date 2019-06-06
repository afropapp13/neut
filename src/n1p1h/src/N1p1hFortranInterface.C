#include "N1p1h.h"
#include <string> 
#include <math.h>
#include <stdio.h>


N1p1h *n1p1h;

void   N1p1hInitialize(void);

int    AddNuclei(int);

double N1p1hGetCrossSection(int id,int nuclei,double E);

double N1p1hGetKinematics(int id,int nuclei,double E, double dir[3], double p[4][4], int idpart[4], int parent[4] ); 

int AddNuclei(int nuclei);

extern "C" {

  extern struct{
    int numbndn;
    int numbndp;
    int numfrep;
    int numatom;
  } neuttarget_;

#define MAXNE 100 
  
  extern struct{
    int modene;
    int numne;
    int ipne[MAXNE];
    float pne[MAXNE][3];
    int iorgne[MAXNE];
    int iflgne[MAXNE];
    int icrnne[MAXNE];
  } nework_;

  #define MAXVC 100 
  
   extern struct{
     int ibound;
     float pos[MAXVC][3];
   } posinnuc_;

  extern struct {
    float pfsurf;
    float pfmax;
    float vnuini;
    float vnufin;
    float iformlen;
    float fzmu2; 
  } nenupr_;

  extern struct {
    char crstblpath[1024]; 
  }neutfilepath_;

  
  float fn_1p1h_n_(float &E,int &ip, int &mode) {

	double xsec;
	float E_tmp;

    N1p1hInitialize(); 
    
    if ( ( abs(ip) != 12 ) && ( abs(ip) != 14 ) && ( abs(ip) != 16 ) ) {
      printf(" Input IP ( %d ) is not neutrino ",ip);  
      exit(0); 
    }

    if( E <= 0.0 ) {
      printf(" Input Energy ( %f ) is negative ",E);  
      exit(0); 
    }

	if ( E > 9.5 ){
	  E_tmp = 9.5;
	}else{
	  E_tmp = E;
	}
	  
    int nuclei = neuttarget_.numatom;

    AddNuclei(nuclei); 
    
    xsec = N1p1hGetCrossSection(ip,nuclei,E_tmp);
	xsec = xsec * 1e38;

	if ( abs(ip) == 12 ){
	  if (E < 30e-3 ){
		xsec = 0.;
	  }
	}else if ( abs(ip) == 14 ){
	  if (E < 0.150 ){
		xsec = 0.;
	  }
	}else if ( abs(ip) == 16 ){
	  if (E < 3.5 ){
		xsec = 0.;
	  }
	}
	return (float) xsec;
  }


  float ne_1p1h_n_(int &ip,int &mode,float &E,float dir[3],int &ierr) {

   if ( ( abs(ip) != 12 ) && ( abs(ip) != 14 ) && ( abs(ip) != 16 ) ) {
      printf(" Input IP ( %d ) is not neutrino \n ",ip);  
      exit(0); 
    }

    if( E <= 0.0 ) {
      printf(" Input Energy ( %f ) is negative \n",E);  
      exit(0); 
    }
    
    int nuclei = neuttarget_.numatom;

    double p[6][4];
    int idpart[6];
    int parent[6];

    double dd[3];
    dd[0]=dir[0];    dd[1]=dir[1];    dd[2]=dir[2];
    
    ierr = -N1p1hGetKinematics(ip,nuclei,E,dd,p,idpart,parent);

    if ( ierr <= 0 ) ierr = 0; 
    
	if ( ip > 0 ){
	  nework_.modene = 1;
	}else{
	  nework_.modene = -1;
	}
    nework_.numne = 4;
    for( int i = 0; i < nework_.numne; i++ ) {
      nework_.ipne[i] = idpart[i];
      nework_.pne[i][0] = p[i][1]/1000.;
      nework_.pne[i][1] = p[i][2]/1000.;
      nework_.pne[i][2] = p[i][3]/1000.;
      nework_.iorgne[i] = parent[i]+1;
      if( i < 2 ) {
	nework_.iflgne[i] = -1;
	nework_.icrnne[i] = 0; 
      }
      else{
	nework_.iflgne[i] = 0;
	nework_.icrnne[i] = 1;
      }
    }
    
    return 0;
  }

}


void N1p1hInitialize(void) {
  static bool initialized = false; 

  if( !initialized ) {
#ifdef WITH_NEUT
	int i = 0;
	char dirname[1024],tmpstr[1024];
	std::string dirname_s;

	snprintf(tmpstr,sizeof(tmpstr),"%s",neutfilepath_.crstblpath);
	while(tmpstr[i] != ' '){
	  i++;
	}
	if ( i > 1022 ){
	  i = 1023;
	}
	tmpstr[i] = 0;
	snprintf(dirname,sizeof(dirname),"%s/N1p1h/Tables_std",tmpstr);
	dirname_s = dirname;
#endif
    bool applybind = true; 
    n1p1h = new N1p1h();
#ifdef WITH_NEUT
	n1p1h->Initialize(dirname_s);
#else
	n1p1h->Initialize(false);
#endif
	//	n1p1h->InitializeAllNuclei();
  }
    
  initialized = true;

  return; 
}

double N1p1hGetCrossSection(int id,int nuclei,double E){
  N1p1hInitialize(); 
  return n1p1h->IntegralCrossSection(id,nuclei,E);
}

double N1p1hGetKinematics(int id,int nuclei,double Enu, double dir[3], double p[4][4], int idpart[4], int parent[4] ){
  N1p1hInitialize();

  double dd = sqrt(dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2]);
  
  double pnu[4];
  pnu[0] = Enu;
  pnu[1] = Enu*dir[0]/dd;
  pnu[2] = Enu*dir[1]/dd;
  pnu[3] = Enu*dir[2]/dd;

  double R,xsect; 


  if(nuclei==1)
    n1p1h->GenerateHVectors(id,nuclei,pnu,p,idpart,parent,R,xsect);
  else
    n1p1h->GenerateVectors(id,nuclei,pnu,p,idpart,parent,R,xsect);
  
  double phi = (double)random()/(double)RAND_MAX*2.*M_PI;
  double costheta = 2.*(double)random()/(double)RAND_MAX-1;
  double sintheta = sqrt(1.-costheta*costheta); 
  posinnuc_.ibound = 1;
  for( int i = 0; i < 4;i++) { 
    posinnuc_.pos[i][0] = R*sintheta*cos(phi); 
    posinnuc_.pos[i][1] = R*sintheta*sin(phi); 
    posinnuc_.pos[i][2] = R*costheta; 
  }

  return 0.; 
}


int AddNuclei(int nuclei){
  if( !n1p1h->IsNucleiInitialized(nuclei) ) 
    n1p1h->InitializeNucleus(nuclei,true);
}
