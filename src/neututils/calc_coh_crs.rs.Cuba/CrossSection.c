#include "cuba.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double crosin_(double *);
double crosto_(double *);
double getr_(double *, double *, double *);
double dbesi0_(double*);

const double PI     = 3.14159265;
const double HBARC  = 0.19732705;
const double GFERMI = 0.89719167E-7;
const double FPIFAC = 0.93;
const double PIMASS = 0.135;
const double PRMASS = 0.939;
const double AM     = 1.;
const double ATOMN  = 16.;
/*const double ATOMN  = 27.;*/
const double RAD0NU = 1.;
const double BRNTC2 = 1.E-24;
const double BRNTF2 = 1.E+2;
const double FMTCTM = 1.E-13;
const double P2OS   = 1./4.;

double mass_lep;
double glob_x, glob_y, glob_z;

double glob_energy;

int    ITYPE, CCNC;
double ccnc_ratio;
int    LeptonMass;


void   Initialize(int , int );
double Coefficient();
int    SetParameters(const int *, const double *);
static void DifferentialCrossSection(const int    *, 
									 const double *,
									 const int    *, 
									 double *);
double TotalCrossSection(const double);


int
main( int argc, char **argv )
{
  
  int nu_type = atoi(argv[1]);
  int current; //  = atoi(argv[2]); // CC(1), NC(0)

  current = 0;

  const int MaxiEnergy = 370;

  int    iEnergy;
  double Energy = 0;
  double xsec;

  FILE *output;

  char fname[128];

  for ( current = 1 ; current < 2 ; current++){
	
	if (current == 0){
	  snprintf(fname,sizeof(fname),"cross.%.2d.NC.dat",nu_type);
	}else if (current == 1){
	  snprintf(fname,sizeof(fname),"cross.%.2d.CC.dat",nu_type);
	}else{
	  printf("err?: current=%d\n",current);
	  exit(1);
	}

	output = fopen(fname,"w");
	
	Initialize(nu_type, current);
	
	for (iEnergy=1; iEnergy<=MaxiEnergy; iEnergy++) {
	  if (iEnergy < 100) {
		Energy = (double)(iEnergy)*0.01;           // 0 - 1GeV
	  } else if (iEnergy >= 100 && iEnergy < 190) {
		Energy = (double)(iEnergy-100)*0.1 + 1.0;  // 1 - 10GeV
	  } else if (iEnergy >= 190 && iEnergy < 280) {
		Energy = (double)(iEnergy-190)*1. + 10.0;  // 10 - 100GeV
	  } else if (iEnergy >= 280 && iEnergy <= 370) {
		Energy = (double)(iEnergy-280)*10. + 100.; // 100 - 1000GeV
	  }      
	  
	  if (Energy != 0.) {
		xsec = TotalCrossSection(Energy);
	  } else {
		xsec = 0.;
	  }
	  printf("%9.5e %9.5e\n", Energy, xsec);
	  fprintf(output, "%9.5e %9.5e\n", Energy, xsec);
	  fflush(output);
	}
	fclose(output);
  }
  
  exit(0);
}


void 
Initialize(int nu_type, int current){

  ITYPE = nu_type;
  CCNC  = current;

  if (ITYPE == 12) {
    mass_lep = 0.000511;
  } else if (ITYPE == 14) {
    mass_lep = 0.1057;
  } else if (ITYPE == 16) {
    mass_lep = 1.7770;
  }

  if (CCNC == 1) {
    ccnc_ratio = 2.;
  } else if (CCNC == 0) {
    ccnc_ratio = 1.;
  } else {
    ccnc_ratio = 0.;
  }

  LeptonMass = 1;  // Lepton mass effect on(1), off(0)

}

///////////////////////////////////////////////////
double 
Coefficient(){

  double bte, coeff;

  bte = pow(RAD0NU,2.)*pow(ATOMN,2./3.)*pow(glob_energy,2.)/(3.*pow(HBARC,2.));

  coeff = PRMASS*pow(ATOMN*PIMASS*GFERMI,2.)*pow(glob_energy,3.)*pow(FPIFAC,2.)
    /pow(PI,2.)/pow(HBARC,8.)*exp(bte*pow(PIMASS/glob_energy,2.))
	*pow(FMTCTM,2.);

  return coeff;
}


///////////////////////////////////////////////////
int 
SetParameters(const int *ndim, 
			  const double xx[]){

  double xhig, xlow;
  double yhig, ylow;
  double zhig, zlow;

  // z; minimum, maximum
  zhig =  1.;
  zlow = -1.;
  glob_z = (zhig - zlow)*xx[0] + zlow;

  // x; minimum, maximum
  xhig = 1.;
  xlow = 0.;
  glob_x = (xhig - xlow)*xx[1] + xlow;
  if (glob_x == 0.) glob_x = 0.0001;
  if (glob_x == 1.) glob_x = 0.9999;

  // y; minimum, maximum
  yhig = 1./(1.+PRMASS*glob_x/(2.*glob_energy));
  ylow = PIMASS/glob_energy;
  glob_y = (yhig - ylow)*xx[2] + ylow;

  return 0;

}

static void
DifferentialCrossSection(const int    *ndim, 
						 const double *xx,
						 const int    *ncomp, 
						 double *ff){

  double dsigdt;
  double qsmx,qsmn,xmax,xmin,y2mpe2,bte,soln,soln_tmp,epitot,epikin;
  double arg1,arg2,arg3,femto1,femto2;
  double ymin,ymax,q2min,q2max,q2,correction;
  double r, bess;

  int    ret;

  ff[0] = 0.;
	
  ret = SetParameters(ndim, xx);

  if(glob_y < PIMASS/glob_energy) return;
  if(glob_y > 1./(1.+PRMASS*glob_x/(2.*glob_energy))) return;

  qsmx = 4.*glob_energy*glob_energy*(1.-glob_y);
  qsmn = 0.;
  if(qsmx < 0.) return;

  xmax = qsmx/(2.*PRMASS*glob_y*glob_energy);
  xmin = qsmn/(2.*PRMASS*glob_y*glob_energy);
  if (glob_x > xmax) return;
  if (glob_x < xmin) return;

  //--- Lepton mass effects
  ymin = PIMASS/glob_energy;
  ymax = 1.-mass_lep/glob_energy;

  q2min = pow(mass_lep,2.)*glob_y/(1-glob_y);
  q2max = 2.*PRMASS*glob_energy*(1.-mass_lep/glob_energy);
  q2 = 2.*PRMASS*glob_energy*glob_x*glob_y;

  correction 
	= ( 1.-q2min/(2.*(q2+pow(PIMASS,2.))) *
	    1.-q2min/(2.*(q2+pow(PIMASS,2.))) )
    + glob_y/4.*q2min*(q2-q2min)/pow(q2+pow(PIMASS,2.),2.);

  if (glob_y < ymin || glob_y > ymax 
	  || q2 < q2min || q2 > q2max) correction = 0.;
  //-----

  y2mpe2 = pow(glob_y,2.) - pow(PIMASS/glob_energy,2.);
  bte = pow(RAD0NU,2.)*pow(ATOMN,2./3.)*pow(glob_energy,2.)/(3.*pow(HBARC,2.));

  soln_tmp 
	= sqrt(fabs(glob_y*(glob_y+2.*PRMASS*glob_x/glob_energy)
				*fabs(y2mpe2)))
    *(1.-glob_y)/pow(1.+2.*PRMASS*glob_energy*glob_x*glob_y/pow(AM,2.),2.);
  epitot = glob_energy*glob_y;
  epikin = glob_energy*glob_y-PIMASS;
  if(epikin < 0.) return;

  arg1 
	= -2.*bte*glob_y*(glob_y+PRMASS*glob_x/glob_energy-
					  glob_z*(1+PRMASS*glob_x/glob_energy)*sqrt(fabs(y2mpe2)));
  arg2 
	= -2.*bte*sqrt(fabs(((1.-pow(glob_z,2.))*y2mpe2
						 *(2.*PRMASS/glob_energy)*glob_x*glob_y
						 *(1.-glob_y*(1.+PRMASS*glob_x/(2.*glob_energy))))));
  
  femto1 = crosin_(&epikin);
  arg3 = -9.*pow(ATOMN,1./3.)*femto1/(16.*PI*pow(RAD0NU,2.));

  soln = soln_tmp*exp(arg1+arg3);

  femto2 = crosto_(&epikin);
  r = getr_(&epikin,&epitot,&femto2);
  dsigdt = P2OS*pow(femto2,2.)*(1.+pow(r,2.))/(4.*PI);

  if (arg2 < -50.) arg2 = -50.;
  bess = dbesi0_(&arg2);

  if (CCNC != 1 || LeptonMass != 1) correction = 1.;
  ff[0] = ccnc_ratio*soln*bess*dsigdt*correction;

  ff[0] = ff[0] * (1./(1.+PRMASS*glob_x/(2.*glob_energy))
				   -PIMASS/glob_energy);

  ff[0] = ff[0] * Coefficient()*1.E+38*2.;

  return;

}

///////////////////////////////////////////////////
double 
TotalCrossSection(const double arg_energy) {

#define NDIM 3
#define NCOMP 1
#define EPSREL 5e-4
#define EPSABS 1e-12
#define VERBOSE 0
#define LAST 4
#define MINEVAL 0
#define MAXEVAL 65536*32767

#define KEY 0

  int    nregions, neval, fail;
  double integral, error, prob;

  glob_energy = arg_energy;

  Cuhre(NDIM, NCOMP, DifferentialCrossSection,
    EPSREL, EPSABS, VERBOSE | LAST, MINEVAL, MAXEVAL,
    KEY,
    &nregions, &neval, &fail, &integral, &error, &prob);

  printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n",
    nregions, neval, fail);
  printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      integral, error, prob);

  return integral;

}
