/*
 * Generated automatically by fh2h.pl
 * !!! DO NOT EDIT !!!
 * Edit the original fortran header file instead
 * or fix fh2h.pl if there is a translation bug.
 */


#ifndef FH2H_UHDEF_INC
#define FH2H_UHDEF_INC


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

/*HGEN*ECHO(FORTRAN)*/

/* 9001-9010 for testplane ntuple*/

/* Ntuple variables for neutrino at SK*/
/*      real EnuSK, ppiSK, xpiSK(3), npiSK(3), cospibmSK, normSK*/
/*      real ppi0SK, xpi0SK(3), npi0SK(3), cospi0bmSK*/
/*      integer modeSK, ppidSK, nvtx0SK*/
/*      common /nusk/ EnuSK, ppidSK, modeSK, ppiSK, xpiSK, npiSK, cospibmSK, */
/*     &	normSK, nvtx0SK, ppi0SK, xpi0SK, npi0SK, cospi0bmSK*/

/* Ntuple variables for neutrino at SK*/
/* Modified by H.Kubo 2009/01/04 make it optional*/
/* Modified by NA - 25/03/09*/
/* Modified by NA - 17/04/09*/
/* Modified by NA - 20/04/09*/
/* Modified by NA - 03/07/09*/
/*common nusk was here*/

/*      character*60 nusktags(6)*/
/*      character*360 nusktag*/
/*      equivalence (nusktag, nusktags)*/
/*      data nusktags/ */
/*     &	   'Enu:R             ,ppid:I            ,mode:I            ,',*/
/*     &     'ppi:R             ,xpi(3):R          ,npi(3):R          ,',*/
/*     &     'cospibm:R         ,norm:R            ,nvtx0:I           ,',*/
/*     &     'ppi0:R            ,xpi0(3):R         ,npi0(3):R         ,',*/
/*     &     'cospi0bm:R        ,gipart[0,100]:I   ,gpos0(3):R        ,',*/
/*     &     'gvec0(3):R        ,gamom0:R '/ */
/**/
/*      character*60 nusktags2(9)*/
/*      character*900 nusktag2*/
/*      equivalence (nusktag2, nusktag)*/
/*      equivalence (nusktag2(361:900), nusktags2)*/
/*      data nusktags2/ */
/*     &     '                                     ,spid:I            ,',*/
/*     &     'pgen:I            ,psi0:R            ,xsi0(3):R         ,',*/
/*     &     'nsi0(3):R         ,cossi0bm:R        ,xsi(3):R          ,',*/
/*     &     'smech:I           ,prvtx(3):R        ,intgt:I           ,',*/
/*     &     'smed:I            ,gppid:I           ,xgpi0(3):R        ,',*/
/*     &     'xgpi(3):R         ,pgpi0:R           ,gpmech:I          ,',*/
/*     &     'gpmed:I           ,prmech:I          ,prmed:I           ,',*/
/*     &     'prdght:I          ,sdght:I           ,gpdght:I          ,',*/
/*     &     'chain:I            '/   */

/* common variables for the Detector setting*/
/* beam center at ZFD*/
/*common deffd was here*/
/*      character*180 defFDtag*/
/*      character*60  defFDtags(3)*/
/*      equivalence (defFDtag, defFDtags)*/
/*      data defFDtags/ */
/*     &     'NFD[0,10]:I       ,BXFD(NFD):R       ,BYFD(NFD):R       ,',*/
/*     &     'XFD(NFD):R        ,YFD(NFD):R        ,ZFD(NFD):R        ,',*/
/*     &     'HFD(NFD):R        ,VFD(NFD):R                            '/ */

/* Ntuple variables for neutrino at Near detectors*/
/*      real EnuFD, ppiFD, xpiFD(3), npiFD(3), cospibmFD, normFD*/
/*      real ppi0FD, xpi0FD(3), npi0FD(3), cospi0bmFD*/
/*      real rFD, xnuFD, ynuFD, nnuFD(3)*/
/*      integer modeFD, ppidFD, nvtx0FD*/
/*      integer idFD*/
/*      common /nuFD/ EnuFD, ppidFD, modeFD, ppiFD, xpiFD, npiFD, */
/*     &	cospibmFD, normFD, nvtx0FD, ppi0FD, xpi0FD, npi0FD, cospi0bmFD, */
/*     &	rFD, xnuFD, ynuFD, nnuFD, idFD*/

/* Ntuple variables for neutrino at Near detectors*/
/* Modified by H.Kubo 2009/01/04 make it optional*/
/* Modified by NA - 25/03/09*/
/* Modified by NA - 17/04/09*/
/* Modified by NA - 20/04/09*/
/* Modified by NA - 03/07/09*/

/*common nufd was here*/

/*      character*60 nuFDtags(8)*/
/*      character*480 nuFDtag*/
/*      equivalence (nuFDtag, nuFDtags)*/
/*      data nuFDtags/ */
/*     &	   'Enu:R             ,ppid:I            ,mode:I            ,',*/
/*     &     'ppi:R             ,xpi(3):R          ,npi(3):R          ,',*/
/*     &     'cospibm:R         ,norm:R            ,nvtx0:I           ,',*/
/*     &     'ppi0:R            ,xpi0(3):R         ,npi0(3):R         ,',*/
/*     &     'cospi0bm:R        ,rnu:R             ,xnu:R             ,',*/
/*     &     'ynu:R             ,nnu(3):R          ,idFD:I            ,',*/
/*     &     'gipart[0,100]:I   ,gpos0(3):R        ,gvec0(3):R        ,',*/
/*     &     'gamom0:R '/ */

/*      character*60 nuFDtags2(8)*/
/*      character*960 nuFDtag2*/
/*      equivalence (nuFDtag2, nuFDtag)*/
/*      equivalence (nuFDtag2(481:960), nuFDtags2)*/
/*      data nuFDtags2/ */
/*     &     '                  ,spid:I            ,pgen:I            ,',*/
/*     &     'psi0:R            ,xsi0(3):R         ,nsi0(3):R         ,',*/
/*     &     'cossi0bm:R        ,xsi(3):R          ,smech:I           ,',*/
/*     &     'intgt:I           ,prvtx(3):R        ,smed:I            ,',*/
/*     &     'gppid:I           ,xgpi0(3):R        ,xgpi(3):R         ,',*/
/*     &     'pgpi0:R           ,gpmech:I          ,gpmed:I           ,',*/
/*     &     'prmech:I          ,prmed:I           ,prdght:I          ,',*/
/*     &     'sdght:I           ,gpdght:I          ,chain:I            '/ */

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
} nusk_;
#ifndef NO_EXTERN_COMMON_POINTERS
extern struct nusk_common *nusk;
#endif
#ifdef STATIC_COMMON_POINTERS
static struct nusk_common *nusk = &nusk_;
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


#endif  /* #ifndef FH2H_UHDEF_INC */
