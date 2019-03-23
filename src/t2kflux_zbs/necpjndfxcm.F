      SUBROUTINE necpjndfxcm(iptr)

      implicit none

      integer*4 iptr

#ifdef FLUX_10A
#include "uhdef.fh"
#else
#ifdef FLUX_10C
#include "uhdef_10c.fh"
#else
#ifdef FLUX_11A
#include "uhdef_11a.fh"
#else
#include "beamntpl.h"
#endif      
#endif      
#endif

#include "beamvectbl.h"

      integer*4 i

      if ((iptr.gt.ITBLMAX).or.(iptr.le.0)) then
         write(*,*) 'necpjndfxcm: subscript out of range...'
         write(*,*) 'necpjndfxcm: iptr=',iptr,
     $        ' must be smaller than',ITBLMAX
         stop
      endif

      EnuFD     = Enutbl(iptr)     
      ppiFD     = ppitbl(iptr)     
      CospibmFD = Cospibmtbl(iptr) 
      normFD    = normtbl(iptr)    
      ppi0FD    = ppi0tbl(iptr)    
      cospi0bmFD= cospi0bmtbl(iptr)
      rFD       = rnutbl(iptr)     
      xnuFD     = xnutbl(iptr)     
      ynuFD     = ynutbl(iptr)     
      ppidFD    = ppidtbl(iptr)    
      modeFD    = modetbl(iptr)    
      idFD      = idFDtbl(iptr)    

#ifdef FLUX_10A
      psi0FD=     psi0tbl(iptr)
      cossi0bmFD= cossi0bmtbl(iptr)

      nvtx0FD=    nvtx0tbl(iptr)
      spidFD=     spidtbl(iptr)
      pgenFD=     pgentbl(iptr)
      intgtFD=    intgttbl(iptr)
      smechFD=    smechtbl(iptr)
      smedFD=     smedtbl(iptr)
      gppidFD=    gppidtbl(iptr)
      pgpi0FD=    pgpi0tbl(iptr)
      gpmechFD=   gpmechtbl(iptr)
      gpmedFD=    gpmedtbl(iptr)
      prmechFD=   prmechtbl(iptr)
      prmedFD=    prmedtbl(iptr)
      prdghtFD=   prdghttbl(iptr)
      sdghtFD=    sdghttbl(iptr)
      gpdghtFD=   gpdghttbl(iptr)
      chainFD=    chaintbl(iptr)
      gipartFD=   giparttbl(iptr)
      gamom0FD=   gamom0tbl(iptr)
#endif
#if defined(FLUX_10C)||defined(FLUX_11A)
      nvtx0FD=    nvtx0tbl(iptr)
      gipartFD=   giparttbl(iptr)
      gamom0FD=   gamom0tbl(iptr)

      ngFD=       ngtbl(iptr)

      do 80 I=1,ngFD
         gpxFD(I)   =gpxtbl(I,iptr)
         gpyFD(I)   =gpytbl(I,iptr)
         gpzFD(I)   =gpztbl(I,iptr)
         gcosbmFD(I)=gcosbmtbl(I,iptr)
         gvxFD(I)   =gvxtbl(I,iptr)
         gvyFD(I)   =gvytbl(I,iptr)
         gvzFD(I)   =gvztbl(I,iptr)
         gpidFD(I)  =gpidtbl(I,iptr)
         gmecFD(I)  =gmectbl(I,iptr)
#ifdef FLUX_11A
         gmatFD(I)   =gmatFDtbl(I,iptr)
         gdistcFD(I) =gdistcFDtbl(I,iptr) 
         gdistalFD(I)=gdistalFDtbl(I,iptr)
         gdisttiFD(I)=gdisttiFDtbl(I,iptr)
         gdistfeFD(I)=gdistfeFDtbl(I,iptr)
#endif
 80   continue

      EnuSKFD=    EnuSKtbl(iptr)
      normSKFD=   normSKtbl(iptr)
      anormFD=    anormtbl(iptr)
          
#endif

      do 90 I=1,3

         npiFD(i)  = npitbl(i,iptr)     
         xpiFD(i)  = xpitbl(i,iptr)     
         npi0FD(i) = npi0tbl(i,iptr)    
         xpi0FD(i) = xpi0tbl(i,iptr)    
         nnuFD(i)  = nnutbl(i,iptr)     

#ifdef FLUX_10A
         xsi0FD(i)=  xsi0tbl(i,iptr)
         nsi0FD(i)=  nsi0tbl(i,iptr)
         xsiFD(i)=   xsitbl(i,iptr)
         prvtxFD(i)= prvtxtbl(i,iptr)
         xgpi0FD(i)= xgpi0tbl(i,iptr)
         xgpiFD(i)=  xgpitbl(i,iptr)
         gpos0FD(i)= gpos0tbl(i,iptr)
         gvec0FD(i)= gvec0tbl(i,iptr)
#endif
#ifdef FLUX_10C
         gpos0FD(I)=    gpos0tbl(i,iptr)
         gvec0FD(I)=    gvec0tbl(i,iptr)
#endif         

 90   continue

      return
      end
      
