#include "HT2p2h.h"
#include <string> 
#include <math.h>
#include <stdio.h>

HT2p2h *ht;

void ht2p2h_initialize(void);

double GetCrossSection(int id,int nuclei,double E);

double GetKinematics(int id,int nuclei,double E, double dir[3], double p[6][4], int idpart[6], int parent[6] ); 

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


  // Specific HT2p2h common block.
  extern struct {
    char model[64];
    int  LFG;
    int  hadrestframemodel; // Which model to use for the hadron rest frame definition. 
    double alphaiso; // These parameters changes the direction of the leading nucleon in the restframe (alpha=beta=0 is isotropic).
    double betaiso;
  }ht2p2h_;


  
  float fn_2p2h_ht_(float &E,int &ip, int &mode) {

	float cross_section;

    ht2p2h_initialize(); 
    
    if ( ( abs(ip) != 12 ) && ( abs(ip) != 14 ) && ( abs(ip) != 16 ) ) {
      printf(" Input IP ( %d ) is not neutrino ",ip);  
      exit(0); 
    }

    if( E <= 0.0 ) {
      printf(" Input Energy ( %f ) is negative ",E);  
      exit(0); 
    }


    int nuclei = neuttarget_.numatom;
    
    if( E < 0.2 ) {
	  cross_section = 0.;
	}else if( E > 30.0 ) {
      //	  cross_section = (float) GetCrossSection(ip,nuclei,30.0);     	  cross_section = 0;
	  
    }else{
	  cross_section = (float) GetCrossSection(ip,nuclei,E);
	}

	cross_section = cross_section / (float)(nuclei);

    return cross_section;
  }


  float ne_2p2h_ht_(int &ip,int &mode,float &E,float dir[3],int &ierr) {

   if ( ( abs(ip) != 12 ) && ( abs(ip) != 14 ) && ( abs(ip) != 16 ) ) {
      printf(" Input IP ( %d ) is not neutrino \n ",ip);  
      exit(0); 
    }

    if( E <= 0.0 ) {
      printf(" Input Energy ( %f ) is negative \n",E);  
      exit(0); 
    }

    //printf( " >>>>>>>>>> FERMI LEVEL %f \n ", nenupr_.pfmax); 
    
    int nuclei = neuttarget_.numatom;

    double p[6][4];
    int idpart[6];
    int parent[6];

    double dd[3];
    dd[0]=dir[0];    dd[1]=dir[1];    dd[2]=dir[2];
    
    ierr = -GetKinematics(ip,nuclei,E,dd,p,idpart,parent);

    if ( ierr <= 0 ) ierr = 0; 
    
    if ( ip > 0){
      nework_.modene = 2;
    }else{
      nework_.modene = -2;
    }
    nework_.numne = 6;
    for( int i = 0; i < nework_.numne; i++ ) {
      nework_.ipne[i] = idpart[i];
      nework_.pne[i][0] = p[i][1];
      nework_.pne[i][1] = p[i][2];
      nework_.pne[i][2] = p[i][3];
      nework_.iorgne[i] = parent[i];
      if( i < 3 ) {
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


void ht2p2h_initialize(void) {
  static bool initialized = false; 

  if( !initialized ) {
	char tmpstr[1024];
	int i=0;

	snprintf(tmpstr,sizeof(tmpstr),"%s",neutfilepath_.crstblpath);
	while(tmpstr[i] != ' '){
	  i++;
	}
	if ( i > 1022 ){
	  i = 1023;
	}
	tmpstr[i] = 0;

    std::string directorypath = std::string(tmpstr)+"/HT2p2h/Nieves";
    std::cout << " Reading from " << directorypath << std::endl;
    ht = new HT2p2h();
	ht->Initialize(directorypath,true, false);
  }
    
  initialized = true;

  return; 
}

double GetCrossSection(int id,int nuclei,double E){
  ht2p2h_initialize(); 
  return ht->IntegralCrossSection(id,nuclei,E);
}

double GetKinematics(int id,int nuclei,double Enu, double dir[3], double p[6][4], int idpart[6], int parent[6] ){
  ht2p2h_initialize();

  double dd = sqrt(dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2]);
  
  double pnu[4];
  pnu[0] = Enu;
  pnu[1] = Enu*dir[0]/dd;
  pnu[2] = Enu*dir[1]/dd;
  pnu[3] = Enu*dir[2]/dd;

  double R; 
  
  ht->GenerateVectors(id,nuclei,pnu,p,idpart,parent,R);

  double phi = (double)random()/(double)RAND_MAX*2.*pi;
  double costheta = 2.*(double)random()/(double)RAND_MAX-1;
  double sintheta = sqrt(1.-costheta*costheta); 
  posinnuc_.ibound = 1;
  for( int i = 0; i < 6;i++) { 
    posinnuc_.pos[i][0] = R*sintheta*cos(phi); 
    posinnuc_.pos[i][1] = R*sintheta*sin(phi); 
    posinnuc_.pos[i][2] = R*costheta; 
  }

  return 0.; 
}
