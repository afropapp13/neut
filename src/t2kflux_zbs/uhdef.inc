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

CHGEN*ECHO(FORTRAN)

C 9001-9010 for testplane ntuple

C Ntuple variables for neutrino at SK
c      real EnuSK, ppiSK, xpiSK(3), npiSK(3), cospibmSK, normSK
c      real ppi0SK, xpi0SK(3), npi0SK(3), cospi0bmSK
c      integer modeSK, ppidSK, nvtx0SK
c      common /nusk/ EnuSK, ppidSK, modeSK, ppiSK, xpiSK, npiSK, cospibmSK, 
c     &	normSK, nvtx0SK, ppi0SK, xpi0SK, npi0SK, cospi0bmSK

C Ntuple variables for neutrino at SK
c Modified by H.Kubo 2009/01/04 make it optional
c Modified by NA - 25/03/09
c Modified by NA - 17/04/09
c Modified by NA - 20/04/09
c Modified by NA - 03/07/09
      real EnuSK, ppiSK, xpiSK(3), npiSK(3), cospibmSK, normSK
      real ppi0SK, xpi0SK(3), npi0SK(3), cospi0bmSK
      real psi0SK, xsi0SK(3), nsi0SK(3), cossi0bmSK, xsiSK(3)
      integer modeSK, ppidSK, nvtx0SK, spidSK, pgenSK
      integer smechSK, intgtSK, smedSK
      real prvtxSK(3)
      integer gppidSK
      real xgpi0SK(3), xgpiSK(3), pgpi0SK
      integer gpmechSK, gpmedSK, prmechSK, prmedSK
      integer prdghtSK, sdghtSK, gpdghtSK, chainSK
      integer gipartSK
      real gpos0SK(3), gvec0SK(3), gamom0SK
      common /nusk/ EnuSK, ppidSK, modeSK, ppiSK, xpiSK, npiSK, cospibmSK, 
     &	normSK, nvtx0SK, ppi0SK, xpi0SK, npi0SK, cospi0bmSK, 
     &  gipartSK, gpos0SK, gvec0SK, gamom0SK, spidSK,
     &  pgenSK, psi0SK, xsi0SK, nsi0SK, cossi0bmSK, xsiSK, smechSK,
     &  prvtxSK, intgtSK, smedSK, gppidSK, xgpi0SK, xgpiSK, pgpi0SK,
     &  gpmechSK, gpmedSK, prmechSK, prmedSK,
     &  prdghtSK, sdghtSK, gpdghtSK, chainSK

C      character*60 nusktags(6)
C      character*360 nusktag
C      equivalence (nusktag, nusktags)
C      data nusktags/
C     &	   'Enu:R             ,ppid:I            ,mode:I            ,',
C     &     'ppi:R             ,xpi(3):R          ,npi(3):R          ,',
C     &     'cospibm:R         ,norm:R            ,nvtx0:I           ,',
C     &     'ppi0:R            ,xpi0(3):R         ,npi0(3):R         ,',
C     &     'cospi0bm:R        ,gipart[0,100]:I   ,gpos0(3):R        ,',
C     &     'gvec0(3):R        ,gamom0:R '/
C
C      character*60 nusktags2(9)
C      character*900 nusktag2
C      equivalence (nusktag2, nusktag)
C      equivalence (nusktag2(361:900), nusktags2)
C      data nusktags2/
C     &     '                                     ,spid:I            ,',
C     &     'pgen:I            ,psi0:R            ,xsi0(3):R         ,',
C     &     'nsi0(3):R         ,cossi0bm:R        ,xsi(3):R          ,',
C     &     'smech:I           ,prvtx(3):R        ,intgt:I           ,',
C     &     'smed:I            ,gppid:I           ,xgpi0(3):R        ,',
C     &     'xgpi(3):R         ,pgpi0:R           ,gpmech:I          ,',
C     &     'gpmed:I           ,prmech:I          ,prmed:I           ,',
C     &     'prdght:I          ,sdght:I           ,gpdght:I          ,',
C     &     'chain:I            '/   

C common variables for the Detector setting
      integer NFD	
      real BXFD(MAX_FD), BYFD(MAX_FD) ! beam center at ZFD
      real XFD(MAX_FD), YFD(MAX_FD), ZFD(MAX_FD), HFD(MAX_FD), VFD(MAX_FD)
      common /DEFFD/NFD, BXFD, BYFD, XFD, YFD, ZFD, HFD, VFD
C      character*180 defFDtag
C      character*60  defFDtags(3)
C      equivalence (defFDtag, defFDtags)
C      data defFDtags/
C     &     'NFD[0,10]:I       ,BXFD(NFD):R       ,BYFD(NFD):R       ,',
C     &     'XFD(NFD):R        ,YFD(NFD):R        ,ZFD(NFD):R        ,',
C     &     'HFD(NFD):R        ,VFD(NFD):R                            '/

C Ntuple variables for neutrino at Near detectors
c      real EnuFD, ppiFD, xpiFD(3), npiFD(3), cospibmFD, normFD
c      real ppi0FD, xpi0FD(3), npi0FD(3), cospi0bmFD
c      real rFD, xnuFD, ynuFD, nnuFD(3)
c      integer modeFD, ppidFD, nvtx0FD
c      integer idFD
c      common /nuFD/ EnuFD, ppidFD, modeFD, ppiFD, xpiFD, npiFD, 
c     &	cospibmFD, normFD, nvtx0FD, ppi0FD, xpi0FD, npi0FD, cospi0bmFD, 
c     &	rFD, xnuFD, ynuFD, nnuFD, idFD

C Ntuple variables for neutrino at Near detectors
c Modified by H.Kubo 2009/01/04 make it optional
c Modified by NA - 25/03/09
c Modified by NA - 17/04/09
c Modified by NA - 20/04/09
c Modified by NA - 03/07/09
      real EnuFD, ppiFD, xpiFD(3), npiFD(3), cospibmFD, normFD
      real ppi0FD, xpi0FD(3), npi0FD(3), cospi0bmFD
      real psi0FD, xsi0FD(3), nsi0FD(3), cossi0bmFD, xsiFD(3)
      real rFD, xnuFD, ynuFD, nnuFD(3)
      integer modeFD, ppidFD, nvtx0FD, spidFD
      integer idFD
      integer pgenFD
      integer intgtFD, smechFD, smedFD
      real prvtxFD(3)
      integer gppidFD
      real xgpi0FD(3), xgpiFD(3), pgpi0FD
      integer gpmechFD, gpmedFD, prmechFD, prmedFD
      integer prdghtFD, sdghtFD, gpdghtFD, chainFD
      integer gipartFD
      real gpos0FD(3), gvec0FD(3), gamom0FD

      common /nuFD/ EnuFD, ppidFD, modeFD, ppiFD, xpiFD, npiFD, 
     &	cospibmFD, normFD, nvtx0FD, ppi0FD, xpi0FD, npi0FD, cospi0bmFD, 
     &	rFD, xnuFD, ynuFD, nnuFD, idFD,
     &  gipartFD, gpos0FD, gvec0FD, gamom0FD,
     &	spidFD, pgenFD, psi0FD, xsi0FD, nsi0FD, cossi0bmFD, xsiFD, smechFD,
     &	intgtFD, prvtxFD, smedFD, gppidFD, xgpi0FD, xgpiFD, pgpi0FD,
     &  gpmechFD, gpmedFD, prmechFD, prmedFD,
     &  prdghtFD, sdghtFD, gpdghtFD, chainFD

C      character*60 nuFDtags(8)
C      character*480 nuFDtag
C      equivalence (nuFDtag, nuFDtags)
C      data nuFDtags/
C     &	   'Enu:R             ,ppid:I            ,mode:I            ,',
C     &     'ppi:R             ,xpi(3):R          ,npi(3):R          ,',
C     &     'cospibm:R         ,norm:R            ,nvtx0:I           ,',
C     &     'ppi0:R            ,xpi0(3):R         ,npi0(3):R         ,',
C     &     'cospi0bm:R        ,rnu:R             ,xnu:R             ,',
C     &     'ynu:R             ,nnu(3):R          ,idFD:I            ,',
C     &     'gipart[0,100]:I   ,gpos0(3):R        ,gvec0(3):R        ,',
C     &     'gamom0:R '/

C      character*60 nuFDtags2(8)
C      character*960 nuFDtag2
C      equivalence (nuFDtag2, nuFDtag)
C      equivalence (nuFDtag2(481:960), nuFDtags2)
C      data nuFDtags2/
C     &     '                  ,spid:I            ,pgen:I            ,',
C     &     'psi0:R            ,xsi0(3):R         ,nsi0(3):R         ,',
C     &     'cossi0bm:R        ,xsi(3):R          ,smech:I           ,',
C     &     'intgt:I           ,prvtx(3):R        ,smed:I            ,',
C     &     'gppid:I           ,xgpi0(3):R        ,xgpi(3):R         ,',
C     &     'pgpi0:R           ,gpmech:I          ,gpmed:I           ,',
C     &     'prmech:I          ,prmed:I           ,prdght:I          ,',
C     &     'sdght:I           ,gpdght:I          ,chain:I            '/

CHGEN*END
