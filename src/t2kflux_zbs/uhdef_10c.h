/*
 * Generated automatically by fh2h.pl
 * !!! DO NOT EDIT !!!
 * Edit the original fortran header file instead
 * or fix fh2h.pl if there is a translation bug.
 */


#ifndef FH2H___UHDEF_10C_NOEQ_H
#define FH2H___UHDEF_10C_NOEQ_H


#ifdef __cplusplus
extern "C" {
#endif


#ifndef IMPLICIT
#define IMPLICIT  /* Only to point out implicit types */
#endif


/*------ fortran header (without commons and data statements) ----------*/

#define FDLOOP 1000
#define MAX_FD 10

#define HID_NUSK 2000
#define HID_NUFD0 3000
#define HID_NUFD1 3001
#define HID_NUFD2 3002
#define HID_ENEDEP_TRGT 5000
#define HID_ENEDEP_DGAS 5001
#define HID_ENEDEP_DUMP 5004
#define HID_ENEDEP_DV 5100
#define HID_MUMON 8000

/*HGEN*ECHO(FORTRAN)*/

/* 9001-9010 for testplane ntuple*/

/* Ntuple variables for neutrino at SK*/
/* Modified by NA - 25/03/09*/
/* Modified by NA - 17/04/09*/
/* Modified by NA - 20/04/09*/
/* Modified by NA - 03/07/09*/
/* Modified by HK - 22/07/10*/

/*---*/

/*---  10a default variables*/
/*common nusk was here*/

/*---  10a optional variables*/
/*common nuskopt was here*/

/*---  10b new variables*/
/*common nusk10b was here*/
        
/* common variables for the Detector setting*/
/* beam center at ZFD*/
/*common deffd was here*/

/* Ntuple variables for neutrino at Near detectors*/
/* Modified by NA - 25/03/09*/
/* Modified by NA - 17/04/09*/
/* Modified by NA - 20/04/09*/
/* Modified by NA - 03/07/09*/
/* Modified by H.K 2010/07/22 */


/*---*/

/*---  10b new variables*/

/*common nufd was here*/

/*HGEN*END*/


/*------ common blocks -------------------------------------------------*/

extern struct nusk_common {
  float  enusk;
  int    ppidsk;
  int    modesk;
  float  ppisk;
  float  xpisk[3];
  float  npisk[3];
  float  cospibmsk;
  float  normsk;
  int    nvtx0sk;
  float  ppi0sk;
  float  xpi0sk[3];
  float  npi0sk[3];
  float  cospi0bmsk;
  int    gipartsk;
  float  gpos0sk[3];
  float  gvec0sk[3];
  float  gamom0sk;
} nusk_;
#ifndef NO_EXTERN_COMMON_POINTERS
extern struct nusk_common *nusk;
#endif
#ifdef STATIC_COMMON_POINTERS
static struct nusk_common *nusk = &nusk_;
#endif

extern struct nuskopt_common {
  int    spidsk;
  int    pgensk;
  float  psi0sk;
  float  xsi0sk[3];
  float  nsi0sk[3];
  float  cossi0bmsk;
  float  xsisk[3];
  int    smechsk;
  float  prvtxsk[3];
  int    intgtsk;
  int    smedsk;
  int    gppidsk;
  float  xgpi0sk[3];
  float  xgpisk[3];
  float  pgpi0sk;
  int    gpmechsk;
  int    gpmedsk;
  int    prmechsk;
  int    prmedsk;
  int    prdghtsk;
  int    sdghtsk;
  int    gpdghtsk;
  int    chainsk;
} nuskopt_;
#ifndef NO_EXTERN_COMMON_POINTERS
extern struct nuskopt_common *nuskopt;
#endif
#ifdef STATIC_COMMON_POINTERS
static struct nuskopt_common *nuskopt = &nuskopt_;
#endif

extern struct nusk10b_common {
  int    ngsk;
  float  gpxsk[20];
  float  gpysk[20];
  float  gpzsk[20];
  float  gcosbmsk[20];
  float  gvxsk[20];
  float  gvysk[20];
  float  gvzsk[20];
  int    gpidsk[20];
  int    gmecsk[20];
} nusk10b_;
#ifndef NO_EXTERN_COMMON_POINTERS
extern struct nusk10b_common *nusk10b;
#endif
#ifdef STATIC_COMMON_POINTERS
static struct nusk10b_common *nusk10b = &nusk10b_;
#endif

extern struct deffd_common {
  int    nfd;
  float  bxfd[MAX_FD];
  float  byfd[MAX_FD];
  float  xfd[MAX_FD];
  float  yfd[MAX_FD];
  float  zfd[MAX_FD];
  float  hfd[MAX_FD];
  float  vfd[MAX_FD];
} deffd_;
#ifndef NO_EXTERN_COMMON_POINTERS
extern struct deffd_common *deffd;
#endif
#ifdef STATIC_COMMON_POINTERS
static struct deffd_common *deffd = &deffd_;
#endif

extern struct nufd_common {
  float  enufd;
  int    ppidfd;
  int    modefd;
  float  ppifd;
  float  xpifd[3];
  float  npifd[3];
  float  cospibmfd;
  float  normfd;
  int    nvtx0fd;
  float  ppi0fd;
  float  xpi0fd[3];
  float  npi0fd[3];
  float  cospi0bmfd;
  float  rfd;
  float  xnufd;
  float  ynufd;
  float  nnufd[3];
  int    idfd;
  int    gipartfd;
  float  gpos0fd[3];
  float  gvec0fd[3];
  float  gamom0fd;
  int    ngfd;
  float  gpxfd[20];
  float  gpyfd[20];
  float  gpzfd[20];
  float  gcosbmfd[20];
  float  gvxfd[20];
  float  gvyfd[20];
  float  gvzfd[20];
  int    gpidfd[20];
  int    gmecfd[20];
  float  enuskfd;
  float  normskfd;
  float  anormfd;
} nufd_;
#ifndef NO_EXTERN_COMMON_POINTERS
extern struct nufd_common *nufd;
#endif
#ifdef STATIC_COMMON_POINTERS
static struct nufd_common *nufd = &nufd_;
#endif


/*------ data statements -----------------------------------------------*/


#ifndef NO_STATIC_DATA


#endif  /* #ifndef NO_STATIC_DATA */


/*------ end of fortran header -----------------------------------------*/


#ifdef __cplusplus
}
#endif


#endif  /* #ifndef FH2H___UHDEF_10C_NOEQ_H */
