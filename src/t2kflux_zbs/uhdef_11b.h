/*
 * Generated automatically by fh2h.pl
 * !!! DO NOT EDIT !!!
 * Edit the original fortran header file instead
 * or fix fh2h.pl if there is a translation bug.
 */


#ifndef FH2H_UHDEF_11B_FH_FOR_FH2H
#define FH2H_UHDEF_11B_FH_FOR_FH2H


#ifdef __cplusplus
extern "C" {
#endif


#ifndef IMPLICIT
#define IMPLICIT  /* Only to point out implicit types */
#endif


/*------ fortran header (without commons and data statements) ----------*/

#define FDLOOP 1000
#define MAX_FD 10

#define HID_VERSET 1000
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


/* common variables for version-tag & setting*/
      

/*common verset was here*/

/*data statement for versettags was here*/

/* Ntuple variables for neutrino at SK*/
/* Modified by NA - 25/03/09*/
/* Modified by NA - 17/04/09*/
/* Modified by NA - 20/04/09*/
/* Modified by NA - 03/07/09*/
/* Modified by HK - 22/07/10*/
/* Modified by M.H. 2011/01/21*/

/*---*/

/*---  10a default variables*/
/*common nusk was here*/

/*data statement for nusktags was here*/

/*---  10a optional variables*/
/*common nuskopt was here*/

/*data statement for nusktagsopt was here*/


/*---  10b new variables*/
/*common nusk10b was here*/

/*data statement for nusktags10b was here*/

/*--- 11a new variables*/
/*common nusk11a was here*/

/*data statement for nusktags11a was here*/


/* common variables for the Detector setting*/
/* beam center at ZFD*/
/*common deffd was here*/
/*data statement for deffdtags was here*/


/* Ntuple variables for neutrino at Near detectors*/
/* Modified by NA - 25/03/09*/
/* Modified by NA - 17/04/09*/
/* Modified by NA - 20/04/09*/
/* Modified by NA - 03/07/09*/
/* Modified by H.K 2010/07/22 */
/* Modified by M.H. 2011/01/21*/
/* Modified by Y.H. 2016/11/05*/


/*---*/
#define NGFDMAX (12)
                

/*---  11b default variables*/

/*common nufd was here*/

/*data statement for nufdtags was here*/

/*---  10a optional variables*/

/*common nufdopt was here*/

/*data statement for nufdtagsopt was here*/

/*HGEN*END*/


/*------ common blocks -------------------------------------------------*/

extern struct verset_common {
  float  set_version;
  int    set_tune;
  int    set_ntrig;
  int    set_pint;
  float  set_bpos[2];
  float  set_btilt[2];
  float  set_brms[2];
  float  set_emit[2];
  float  set_alpha[2];
  float  set_hcur[3];
  int    set_rand;
  int    set_rseed[2];
} verset_;
#ifndef NO_EXTERN_COMMON_POINTERS
extern struct verset_common *verset;
#endif
#ifdef STATIC_COMMON_POINTERS
static struct verset_common *verset = &verset_;
#endif

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
  float  gpxsk[25];
  float  gpysk[25];
  float  gpzsk[25];
  float  gcosbmsk[25];
  float  gvxsk[25];
  float  gvysk[25];
  float  gvzsk[25];
  int    gpidsk[25];
  int    gmecsk[25];
} nusk10b_;
#ifndef NO_EXTERN_COMMON_POINTERS
extern struct nusk10b_common *nusk10b;
#endif
#ifdef STATIC_COMMON_POINTERS
static struct nusk10b_common *nusk10b = &nusk10b_;
#endif

extern struct nusk11a_common {
  int    gmatsk[25];
  float  gdistcsk[25];
  float  gdistalsk[25];
  float  gdisttisk[25];
  float  gdistfesk[25];
} nusk11a_;
#ifndef NO_EXTERN_COMMON_POINTERS
extern struct nusk11a_common *nusk11a;
#endif
#ifdef STATIC_COMMON_POINTERS
static struct nusk11a_common *nusk11a = &nusk11a_;
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
  float  gpxfd[12];
  float  gpyfd[12];
  float  gpzfd[12];
  float  gcosbmfd[12];
  float  gvxfd[12];
  float  gvyfd[12];
  float  gvzfd[12];
  int    gpidfd[12];
  int    gmecfd[12];
  float  enuskfd;
  float  normskfd;
  float  anormfd;
  int    gmatfd[12];
  float  gdistcfd[12];
  float  gdistalfd[12];
  float  gdisttifd[12];
  float  gdistfefd[12];
} nufd_;
#ifndef NO_EXTERN_COMMON_POINTERS
extern struct nufd_common *nufd;
#endif
#ifdef STATIC_COMMON_POINTERS
static struct nufd_common *nufd = &nufd_;
#endif

extern struct nufdopt_common {
  int    spidfd;
  int    pgenfd;
  float  psi0fd;
  float  xsi0fd[3];
  float  nsi0fd[3];
  float  cossi0bmfd;
  float  xsifd[3];
  int    smechfd;
  int    intgtfd;
  float  prvtxfd[3];
  int    smedfd;
  int    gppidfd;
  float  xgpi0fd[3];
  float  xgpifd[3];
  float  pgpi0fd;
  int    gpmechfd;
  int    gpmedfd;
  int    prmechfd;
  int    prmedfd;
  int    prdghtfd;
  int    sdghtfd;
  int    gpdghtfd;
  int    chainfd;
} nufdopt_;
#ifndef NO_EXTERN_COMMON_POINTERS
extern struct nufdopt_common *nufdopt;
#endif
#ifdef STATIC_COMMON_POINTERS
static struct nufdopt_common *nufdopt = &nufdopt_;
#endif


/*------ data statements -----------------------------------------------*/


#ifndef NO_STATIC_DATA

#define VERSETTAGS_MAX_LENGTH (60)
static char   *versettags[4] = {      "version:R         , tuneID:I          , ntrig:I           , ",       "pint:I            , bpos(2):R         , btilt(2):R        , ",       "brms(2):R         , emit(2):R         , alpha(2):R        , ",       "hcur(3):R         , rand:I            , rseed(2):I         "};

#define NUSKTAGS_MAX_LENGTH (60)
static char   *nusktags[6] = {     "Enu:R             , ppid:I            , mode:I            , ",       "ppi:R             , xpi(3):R          , npi(3):R          , ",       "cospibm:R         , norm:R            , nvtx0:I           , ",       "ppi0:R            , xpi0(3):R         , npi0(3):R         , ",       "cospi0bm:R        , gipart[0, 100]:I   , gpos0(3):R        , ",       "gvec0(3):R        , gamom0:R "};

#define NUSKTAGSOPT_MAX_LENGTH (60)
static char   *nusktagsopt[8] = {      "spid:I            , pgen:I            , psi0:R            , ",       "xsi0(3):R         , nsi0(3):R         , cossi0bm:R        , ",       "xsi(3):R          , smech:I           , prvtx(3):R        , ",       "intgt:I           , smed:I            , gppid:I           , ",       "xgpi0(3):R        , xgpi(3):R         , pgpi0:R           , ",       "gpmech:I          , gpmed:I           , prmech:I          , ",       "prmed:I           , prdgh:I           , sdght:I           , ",       "gpdght:I          , chain:I "};

#define NUSKTAGS10B_MAX_LENGTH (60)
static char   *nusktags10b[4] = {      "ng:I::[0, 25]      , gpx(ng):R         , gpy(ng):R         , ",       "gpz(ng):R         , gcosbm(ng):R      , gvx(ng):R         , ",       "gvy(ng):R         , gvz(ng):R         , gpid(ng):I        , ",       "gmec(ng):I   "};

#define NUSKTAGS11A_MAX_LENGTH (60)
static char   *nusktags11a[2] = {      "gmat(ng):I        , gdistc(ng):R      , gdistal(ng):R     , ",       "gdistti(ng):R     , gdistfe(ng):R "};

#define DEFFDTAGS_MAX_LENGTH (60)
static char   *deffdtags[3] = {      "NFD[0, 10]:I       , BXFD(NFD):R       , BYFD(NFD):R       , ",       "XFD(NFD):R        , YFD(NFD):R        , ZFD(NFD):R        , ",       "HFD(NFD):R        , VFD(NFD):R                            "};

#define NUFDTAGS_MAX_LENGTH (60)
static char   *nufdtags[15] = {      "Enu:R             , ppid:I            , modeFD:I          , ",       "ppi:R             , xpi(3):R          , npi(3):R          , ",       "cospibm:R         , norm:R            , nvtx0:I           , ",       "ppi0:R            , xpi0(3):R         , npi0(3):R         , ",       "cospi0bm:R        , rnu:R             , xnu:R             , ",       "ynu:R             , nnu(3):R          , idFD:I            , ",       "gipart[0, 100]:I   , gpos0(3):R        , gvec0(3):R        , ",       "gamom0:R          ,                                       ",       "ng[0, 12]:I        , gpx(ng):R         , gpy(ng):R         , ",       "gpz(ng):R         , gcosbm(ng):R      , gvx(ng):R         , ",       "gvy(ng):R         , gvz(ng):R         , gpid(ng):I        , ",       "gmec(ng):I        , EnuSK:R           , normSK:R          , ",       "anorm:R           ,                                       ",       "gmat(ng):I        , gdistc(ng):R      , gdistal(ng):R     , ",       "gdistti(ng):R     , gdistfe(ng):R                         "};

#define NUFDTAGSOPT_MAX_LENGTH (60)
static char   *nufdtagsopt[8] = {      "                  , spid:I            , pgen:I            , ",       "psi0:R            , xsi0(3):R         , nsi0(3):R         , ",       "cossi0bm:R        , xsi(3):R          , smech:I           , ",       "intgt:I           , prvtx(3):R        , smed:I            , ",       "gppid:I           , xgpi0(3):R        , xgpi(3):R         , ",       "pgpi0:R           , gpmech:I          , gpmed:I           , ",       "prmech:I          , prmed:I           , prdght:I          , ",       "sdght:I           , gpdght:I          , chain:I            "};


#endif  /* #ifndef NO_STATIC_DATA */


/*------ end of fortran header -----------------------------------------*/


#ifdef __cplusplus
}
#endif


#endif  /* #ifndef FH2H_UHDEF_11B_FH_FOR_FH2H */
