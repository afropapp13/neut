#include <TH1D.h>
#include <iostream>
#include <t2kflux_sk.h>
#include <stdlib.h>
#include <TFile.h>


#ifdef FLUX_10A
#include "uhdef.h"
#define FLXCOM  nusk_
#else
#ifdef FLUX_10C
#include "uhdef_10c.h"
#define FLXCOM  nusk_
#else
#ifdef FLUX_11A
#include "uhdef_11a.h"
#define FLXCOM  nusk_
#else
#include "beamntplC.h"
#define FLXCOM  fxvcsk_
#endif
#endif
#endif

using namespace std;

extern "C"
{
  void neopskfxv_(int *, const char *, int *, int *, int );

  double fntotpau_(int *,float *);

  void hnoent_(int *, int *);
  void hgntf_(int *, int *, int *);
  void hrend_(const char *, int);
  
#include <rfa.h>

};


T2Kflux_SK::T2Kflux_SK(TRandom3 *tr3_src = NULL)
{
  memset(fluxhisto,0,sizeof(fluxhisto));
  memset(ratehisto,0,sizeof(ratehisto));
  flux_loaded = 0;
  
  if (tr3_src == NULL){
	cout << "Please initialize tr3_src before calling " 
		 << __PRETTY_FUNCTION__;
	exit(-1);
  }

  tr3 = tr3_src;

};

double
T2Kflux_SK::rand_enu(int ip_flux,int ip_crs, int *err)
{
  double enu;
  int flux_idx, crs_idx;

  if (flux_loaded == 0){
	*err = load_flux();
	if (*err < 0){
	  return 0.;
	}
  }
  flux_idx = pid_idx(ip_flux);
  crs_idx  = pid_idx(ip_crs);
  enu = (unr_rate[flux_idx][crs_idx])->Sample();
  return enu;
};

double
T2Kflux_SK::rand_enuflux(int ip_flux,int *err)
{
  double enu;
  int flux_idx;
  
  if (flux_loaded == 0){
	*err = load_flux();
	if (*err < 0){
	  return 0.;
	}
  }
  flux_idx = pid_idx(ip_flux);
  enu = (unr_flux[flux_idx])->Sample();
  return enu;
};

int
T2Kflux_SK::load_flux_histograms()
{

  Rfile rfile;
  int  ierr,i,j, status;
  char htitle[256],hname[256];
  int  cnpt = 1;
  TFile *histf;

  memset(&rfile,sizeof(rfile),0);

  rgetvalue(T2Kflux_SK::histo_luni,&rfile,cnpt,&status); 

  if (status != 0){
	return -1;
  }
				
  if(rfile.fname[0] == '\0'){
	return -1;
  }
  
  histf = new TFile(rfile.fname,"READ");
  if (histf->IsZombie()){
	cout << "Specified flux histogram file "
		 << rfile.fname
		 << "does not exist." << endl;
	cout << "Try flux hbook file." << endl;
	return -1;
  }

  ierr = 0;
  for ( i = 0 ; i < 4 ; i++ ){
	if (fluxhisto[i] != NULL){
	  delete fluxhisto[i];
	}
	snprintf(hname,sizeof(htitle),"t2k_skflux_%s",
			 flavor_string[i]);
	snprintf(htitle,sizeof(htitle),"t2k_skflux %s;energy",
			 flavor_string[i]);
	fluxhisto[i] = (TH1D *)(histf->Get(hname));
	if (fluxhisto[i] == NULL){
	  ierr = 1;
	}
 	
	for( j = 0 ; j < 4 ; j++){
	  if (ratehisto[i][j] != NULL){
		delete ratehisto[i][j];
	  }
	  snprintf(hname,sizeof(htitle),"skrate_%s_x_%s",
			   flavor_string[i],flavor_string[j]);
	  snprintf(htitle,sizeof(htitle),"skrate %s x %s ;energy",
			   flavor_string[i],flavor_string[j]);
	  ratehisto[i][j] = (TH1D *)(histf->Get(hname));
	  if (ratehisto[i][j] == NULL){
		ierr = 1;
	  }
	}
  }
  if (ierr != 0){
	cout << "Specified flux histogram file"
		 << rfile.fname
		 << "seems to be corrupted." << endl;
	cout << "Try flux hbook file." << endl;
	return -1;
  }
  return 0;
}

int
T2Kflux_SK::load_flux()
{

  int errcode = 0;
  
  int n_ntfiles = 0;
  int nt_entries;
  int crs_index,flux_index,crs_pid;
  int i;

  int histo_loaded = 0;

  Double_t enu_dbl,flux_weight_dbl;
  double crssect,crs_norm;

  static const char chdir[]="skbeam\0";

  reset_histos();

  errcode = load_flux_histograms();
  if (errcode == 0){
	histo_loaded = 1;
  }else{
	histo_loaded = 0;
  }

  while(histo_loaded != 1){
	n_ntfiles++;
	neopskfxv_((int *)(&T2Kflux_SK::luni), chdir, &n_ntfiles, &errcode, 
			   strlen(chdir));
	if (errcode != 0){
	  if (n_ntfiles == 1){
		cout << __PRETTY_FUNCTION__
			 << ": No file was specified in RFLIST";
		return -1;
	  } else {
		n_ntfiles--;
		break;
	  }
	}
	
	hnoent_((int *)(&t2k_skflux_ntid),&nt_entries);
	cout << "Number of entries for this file : " << nt_entries << "\n"; 

	for ( i = 1 ; i <= nt_entries ; i++){
	  hgntf_((int *)(&t2k_skflux_ntid), &i, &errcode);
	  if (errcode != 0){
		cout << __PRETTY_FUNCTION__
			 << ": Failed to read flux ntuple file.";
		if (( i == 1 )&&( n_ntfiles == 1)){
		  return -1;
		}

		cout << __PRETTY_FUNCTION__ 
		     << ": Expected # of events was " << nt_entries
			 << " but only " << i-1 << "events were read in."
			 << " from file #" << n_ntfiles;
		cout << " try to read the next file. ";
		break;
	  }else{
		if (i % 10000 == 0){
		  cout << "now loading event #" << i << "\n";
		}
	  }
	  
	  if (FLXCOM.modesk != 0){
		flux_index = (FLXCOM.modesk/10)-1;
		if (FLXCOM.normsk > 0){
		  fluxhisto[flux_index]->Fill(FLXCOM.enusk,FLXCOM.normsk);
		  for ( crs_index = 0; crs_index < 4 ; crs_index++){
			crs_pid  = pidtbl[crs_index];
			crssect  = fntotpau_(&crs_pid,&(FLXCOM.enusk));
			crs_norm = crssect * FLXCOM.normsk;
			ratehisto[flux_index][crs_index]->Fill(FLXCOM.enusk,crs_norm);
		  }
		}
	  }else{
		cout << "Found FLXCOM.modesk == 0 at Event #" << i << "\n";
	  }
	}
	hrend_(chdir,strlen(chdir));

  }

  for (flux_index = 0 ; flux_index < 4 ; flux_index ++){
	unr_fluxdist[flux_index] = new TUnuranEmpDist(fluxhisto[flux_index]);
	unr_flux[flux_index] = new TUnuran(tr3);
	
	(unr_flux[flux_index])->Init(*(unr_fluxdist[flux_index]),"empk");

	for (crs_index = 0 ; crs_index < 4 ; crs_index ++){
	  unr_ratedist[flux_index][crs_index] 
		= new TUnuranEmpDist(ratehisto[flux_index][crs_index]);
	  unr_rate[flux_index][crs_index] = new TUnuran(tr3);
	  (unr_rate[flux_index][crs_index])
		->Init(*(unr_ratedist[flux_index][crs_index]),"empk");
	}
  }


  flux_loaded = 1;
  return 0;
};

int
T2Kflux_SK::reset_histos()
{
  int i,j;

  char htitle[256],hname[256];

  for ( i = 0 ; i < 4 ; i++ ){
	if (fluxhisto[i] == NULL){
	  snprintf(hname,sizeof(htitle),"t2k_skflux_%s",
			   flavor_string[i]);
	  snprintf(htitle,sizeof(htitle),"t2k_skflux %s;energy",
			   flavor_string[i]);
	  fluxhisto[i] = new TH1D(hname,htitle,SKFLUX_NBINS,0,SKFLUX_EMAX);
	}else{
	  fluxhisto[i]->Reset();
	}
 	
	for( j = 0 ; j < 4 ; j++){
	  if (ratehisto[i][j] == NULL){
		snprintf(hname,sizeof(htitle),"skrate_%s_x_%s",
				 flavor_string[i],flavor_string[j]);
		snprintf(htitle,sizeof(htitle),"skrate %s x %s ;energy",
				 flavor_string[i],flavor_string[j]);
		ratehisto[i][j] = new TH1D(hname,htitle,SKFLUX_NBINS,0,SKFLUX_EMAX);
	  }else{
		ratehisto[i][j]->Reset();
	  }
	}
  }

};

int
T2Kflux_SK::pid_idx(int pid)
{
  int i;
  for ( i = 0 ; i < 4 ; i++ ){
	if (pid == pidtbl[i]) {
	  return i;
	}
  }
  return -1;
};


const char
T2Kflux_SK::flavor_string[4][10] =
{
  "numu","numu_bar","nue","nue_bar"
};

const int
T2Kflux_SK::pidtbl[4] =
{
  14,-14,12,-12
};

const int T2Kflux_SK::luni        = 10;
const int T2Kflux_SK::histo_luni   = 11;
const int T2Kflux_SK::t2k_skflux_ntid = 2000; 
