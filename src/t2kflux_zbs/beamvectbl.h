      integer*4 ITBLMAX
      parameter (ITBLMAX=1000000)
      real Enutbl(ITBLMAX),     ppitbl(ITBLMAX),    xpitbl(3,ITBLMAX)
      real npitbl(3,ITBLMAX),   cospibmtbl(ITBLMAX),normtbl(ITBLMAX)
      real ppi0tbl(ITBLMAX),    xpi0tbl(3,ITBLMAX), npi0tbl(3,ITBLMAX)
      real cospi0bmtbl(ITBLMAX)

      integer modetbl(ITBLMAX), ppidtbl(ITBLMAX),   nvtx0tbl(ITBLMAX)
      integer idFDtbl(ITBLMAX)

#ifdef FLUX_10A
      integer spidtbl(ITBLMAX)
      real psi0tbl(ITBLMAX)
	  real xsi0tbl(3,ITBLMAX), nsi0tbl(3,ITBLMAX)
      real cossi0bmtbl(ITBLMAX),xsitbl(3,ITBLMAX)
      real rnutbl(ITBLMAX),     xnutbl(ITBLMAX),    ynutbl(ITBLMAX)
      real nnutbl(3,ITBLMAX)

      integer pgentbl(ITBLMAX)
      integer intgttbl(ITBLMAX),smechtbl(ITBLMAX), smedtbl(ITBLMAX)
      real prvtxtbl(3,ITBLMAX)
      integer gppidtbl(ITBLMAX)
      real xgpi0tbl(3,ITBLMAX),xgpitbl(3,ITBLMAX), pgpi0tbl(ITBLMAX)
      integer gpmechtbl(ITBLMAX), gpmedtbl(ITBLMAX), prmechtbl(ITBLMAX)
      integer prmedtbl(ITBLMAX)
      integer prdghttbl(ITBLMAX), sdghttbl(ITBLMAX), gpdghttbl(ITBLMAX)
      integer chaintbl(ITBLMAX)
      integer giparttbl(ITBLMAX)
      real gpos0tbl(3,ITBLMAX), gvec0tbl(3,ITBLMAX)
      real gamom0tbl(ITBLMAX)
#endif
#ifdef FLUX_10C
	  real xsi0tbl(3,ITBLMAX), nsi0tbl(3,ITBLMAX)
      real cossi0bmtbl(ITBLMAX),xsitbl(3,ITBLMAX)
      integer pgentbl(ITBLMAX)

      real rnutbl(ITBLMAX),     xnutbl(ITBLMAX),    ynutbl(ITBLMAX)
      real nnutbl(3,ITBLMAX)

      integer giparttbl(ITBLMAX)
      real gpos0tbl(3,ITBLMAX), gvec0tbl(3,ITBLMAX), gamom0tbl(ITBLMAX)
      integer ngtbl(ITBLMAX)

      real    gpxtbl(20,ITBLMAX),gpytbl(20,ITBLMAX)
      real    gpztbl(20,ITBLMAX),gcosbmtbl(20,ITBLMAX)
      real    gvxtbl(20,ITBLMAX),gvytbl(20,ITBLMAX),gvztbl(20,ITBLMAX)
      integer gpidtbl(20,ITBLMAX),gmectbl(20,ITBLMAX)

      real    EnuSKtbl(ITBLMAX), normSKtbl(ITBLMAX)
      real    anormtbl(ITBLMAX)
#endif      
#ifdef FLUX_11A
	  real xsi0tbl(3,ITBLMAX), nsi0tbl(3,ITBLMAX)
      real cossi0bmtbl(ITBLMAX),xsitbl(3,ITBLMAX)
      integer pgentbl(ITBLMAX)

      real rnutbl(ITBLMAX),     xnutbl(ITBLMAX),    ynutbl(ITBLMAX)
      real nnutbl(3,ITBLMAX)

      integer giparttbl(ITBLMAX)
      real gpos0tbl(3,ITBLMAX), gvec0tbl(3,ITBLMAX), gamom0tbl(ITBLMAX)
      integer ngtbl(ITBLMAX)

      real    gpxtbl(25,ITBLMAX),gpytbl(25,ITBLMAX)
      real    gpztbl(25,ITBLMAX),gcosbmtbl(25,ITBLMAX)
      real    gvxtbl(25,ITBLMAX),gvytbl(25,ITBLMAX),gvztbl(25,ITBLMAX)
      integer gpidtbl(25,ITBLMAX),gmectbl(25,ITBLMAX)

      real    EnuSKtbl(ITBLMAX), normSKtbl(ITBLMAX)
      real    anormtbl(ITBLMAX)

      integer gmatFDtbl(25,ITBLMAX)
      real    gdistcFDtbl(25,ITBLMAX),gdistalFDtbl(25,ITBLMAX),
     $     gdisttiFDtbl(25,ITBLMAX),gdistfeFDtbl(25,ITBLMAX)
#endif      
#ifdef FLUX_11B
	  real xsi0tbl(3,ITBLMAX), nsi0tbl(3,ITBLMAX)
      real cossi0bmtbl(ITBLMAX),xsitbl(3,ITBLMAX)
      integer pgentbl(ITBLMAX)

      real rnutbl(ITBLMAX),     xnutbl(ITBLMAX),    ynutbl(ITBLMAX)
      real nnutbl(3,ITBLMAX)

      integer giparttbl(ITBLMAX)
      real gpos0tbl(3,ITBLMAX), gvec0tbl(3,ITBLMAX), gamom0tbl(ITBLMAX)
      integer ngtbl(ITBLMAX)

      real    gpxtbl(12,ITBLMAX),gpytbl(12,ITBLMAX)
      real    gpztbl(12,ITBLMAX),gcosbmtbl(12,ITBLMAX)
      real    gvxtbl(12,ITBLMAX),gvytbl(12,ITBLMAX),gvztbl(12,ITBLMAX)
      integer gpidtbl(12,ITBLMAX),gmectbl(12,ITBLMAX)

      real    EnuSKtbl(ITBLMAX), normSKtbl(ITBLMAX)
      real    anormtbl(ITBLMAX)

      integer gmatFDtbl(12,ITBLMAX)
      real    gdistcFDtbl(12,ITBLMAX),gdistalFDtbl(12,ITBLMAX),
     $     gdisttiFDtbl(12,ITBLMAX),gdistfeFDtbl(12,ITBLMAX)
#endif      
#ifdef FLUX_13
      real rnutbl(ITBLMAX),     xnutbl(ITBLMAX),    ynutbl(ITBLMAX)
      real nnutbl(3,ITBLMAX)

      integer giparttbl(ITBLMAX)
      real gpos0tbl(3,ITBLMAX), gvec0tbl(3,ITBLMAX), gamom0tbl(ITBLMAX)
      integer ngtbl(ITBLMAX)

      real    gpxtbl(12,ITBLMAX),gpytbl(12,ITBLMAX)
      real    gpztbl(12,ITBLMAX),gcosbmtbl(12,ITBLMAX)
      real    gvxtbl(12,ITBLMAX),gvytbl(12,ITBLMAX),gvztbl(12,ITBLMAX)
      integer gpidtbl(12,ITBLMAX),gmectbl(12,ITBLMAX)

      real    EnuSKtbl(ITBLMAX), normSKtbl(ITBLMAX)
      real    anormtbl(ITBLMAX)

      integer gmatFDtbl(12,ITBLMAX)
      real    gdistcFDtbl(12,ITBLMAX),gdistalFDtbl(12,ITBLMAX),
     $     gdisttiFDtbl(12,ITBLMAX),gdistfeFDtbl(12,ITBLMAX)

      real posExitFDtbl(3,ITBLMAX), momExitFDtbl(3,ITBLMAX)
      integer idExitFDtbl(ITBLMAX), ngExitFDtbl(ITBLMAX)

#endif      

#ifdef FLUX_10A
      common /beamvecs/Enutbl,ppidtbl,modetbl,ppitbl,xpitbl,npitbl,
     $     cospibmtbl,normtbl,ppi0tbl,xpi0tbl,npi0tbl,
     $     cospi0bmtbl,rnutbl,xnutbl,ynutbl,nnutbl,
     $     idFDtbl,
     $     xsi0tbl,  nsi0tbl, psi0tbl,  cossi0bmtbl,xsitbl,   nvtx0tbl,
     $     spidtbl,  pgentbl, intgttbl, smechtbl,   smedtbl,  prvtxtbl,
     $     gppidtbl, xgpi0tbl,xgpitbl,  pgpi0tbl,   gpmechtbl,gpmedtbl,
     $     prmechtbl,prmedtbl,prdghttbl,sdghttbl,   gpdghttbl,chaintbl,
     $     giparttbl,gpos0tbl,gvec0tbl, gamom0tbl
#else
#ifdef FLUX_10C
      common /beamvecs/Enutbl,ppidtbl,modetbl,ppitbl,xpitbl,npitbl,
     $     cospibmtbl,normtbl,ppi0tbl,xpi0tbl,npi0tbl,
     $     cospi0bmtbl,rnutbl,xnutbl,ynutbl,nnutbl,
     $     idFDtbl,nvtx0tbl,
     &     giparttbl, gpos0tbl, gvec0tbl, gamom0tbl,
     $     ngtbl, gpxtbl, gpytbl, gpztbl, gcosbmtbl, 
     $     gvxtbl, gvytbl, gvztbl, 
     $     gpidtbl, gmectbl, EnuSKtbl, normSKtbl, anormtbl
#else
#ifdef FLUX_11A
      common /beamvecs/Enutbl,ppidtbl,modetbl,ppitbl,xpitbl,npitbl,
     $     cospibmtbl,normtbl,ppi0tbl,xpi0tbl,npi0tbl,
     $     cospi0bmtbl,rnutbl,xnutbl,ynutbl,nnutbl,
     $     idFDtbl,nvtx0tbl,
     &     giparttbl, gpos0tbl, gvec0tbl, gamom0tbl,
     $     ngtbl, gpxtbl, gpytbl, gpztbl, gcosbmtbl, 
     $     gvxtbl, gvytbl, gvztbl, 
     $     gpidtbl, gmectbl, EnuSKtbl, normSKtbl, anormtbl,
     $     gmatFDtbl, gdistcFDtbl, gdistalFDtbl, 
     $     gdisttiFDtbl, gdistfeFDtbl
#else
#ifdef FLUX_11B
      common /beamvecs/Enutbl,ppidtbl,modetbl,ppitbl,xpitbl,npitbl,
     $     cospibmtbl,normtbl,ppi0tbl,xpi0tbl,npi0tbl,
     $     cospi0bmtbl,rnutbl,xnutbl,ynutbl,nnutbl,
     $     idFDtbl,nvtx0tbl,
     &     giparttbl, gpos0tbl, gvec0tbl, gamom0tbl,
     $     ngtbl, gpxtbl, gpytbl, gpztbl, gcosbmtbl, 
     $     gvxtbl, gvytbl, gvztbl, 
     $     gpidtbl, gmectbl, EnuSKtbl, normSKtbl, anormtbl,
     $     gmatFDtbl, gdistcFDtbl, gdistalFDtbl, 
     $     gdisttiFDtbl, gdistfeFDtbl
#else
#ifdef FLUX_13
      common /beamvecs/Enutbl,ppidtbl,modetbl,ppitbl,xpitbl,npitbl,
     $     cospibmtbl,normtbl,
     $     nvtx0tbl,ppi0tbl,xpi0tbl,npi0tbl,
     $     cospi0bmtbl,rnutbl,xnutbl,ynutbl,nnutbl,
     $     idFDtbl,
     &     giparttbl, gpos0tbl, gvec0tbl, gamom0tbl,
     $     ngtbl, gpxtbl, gpytbl, gpztbl, gcosbmtbl, 
     $     gvxtbl, gvytbl, gvztbl, 
     $     gpidtbl, gmectbl, EnuSKtbl, normSKtbl, anormtbl,
     $     gmatFDtbl, gdistcFDtbl, gdistalFDtbl, 
     $     gdisttiFDtbl, gdistfeFDtbl,
     $     idExitFDtbl, posExitFDtbl, momExitFDtbl, ngExitFDtbl 
#else
      common /beamvecs/Enutbl,ppidtbl,modetbl,ppitbl,xpitbl,npitbl,
     $     cospibmtbl,normtbl,ppi0tbl,xpi0tbl,npi0tbl,
     $     cospi0bmtbl,rnutbl,xnutbl,ynutbl,nnutbl,
     $     idFDtbl
#endif
#endif
#endif
#endif
#endif

	  integer*4  maxfvnpt(4),idxtbl(4,ITBLMAX)
      real*8     totnormtbl(4),normfvtbl(4,ITBLMAX)
      common /beamveci/maxfvnpt,idxtbl,
     $     totnormtbl,normfvtbl
