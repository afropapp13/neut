#include <TH1.h>
#include <TH1D.h>
#include <TRandom3.h>
/*
#include <TUnuran.h>
#include <TUnuranEmpDist.h>
*/

class T2Kflux_SK{
 public:
  T2Kflux_SK(TRandom3 *tr3_src);

  double rand_enu(int ip_flux, int ip_crs, int *err);
  double rand_enuflux(int ip_flux, int *err);

  /* taken from root code but use tr3. 
	 ROOT6 uses TR3 by default but ROOT5 does not.
	 Therefore, we need a custom function to use TR3.
  */
  double GetRandom(TH1D *h); 

  /*                     14 , -14,  12, -12 */
  TH1D *fluxhisto[4];

  /* 2019/09/14 : Stop using TUnuran 
	 because this does not bsupport var. bin widths.
  TUnuranEmpDist *unr_fluxdist[4];
  TUnuran *unr_flux[4];
  */

  /*    first  index : flux
        second index : cross-section
  */

#define SKFLUX_EMAX  30.   /* Maximum energy in GEV */
#define SKFLUX_NBINS 600   /* # of bins */

  TH1D *ratehisto[4][4];

  /* 2019/09/14 : Stop using TUnuran 
	 because this does not bsupport var. bin widths.
  TUnuranEmpDist *unr_ratedist[4][4];
  TUnuran *unr_rate[4][4];
  */

 private:
  static const char flavor_string[4][10];
  static const char beam_flavor_string[4][10];
  static const int pidtbl[4];
  static const int luni;
  static const int histo_luni;
  static const int t2k_skflux_ntid;

  int flux_loaded;
  int load_flux();
  int load_flux_histograms();
  int load_beam_flux_histograms();
  int reset_histos();
  int pid_idx(int pid);

  TRandom3 *tr3;

};
