      SUBROUTINE nerdnufxv(ipbm,idfdet,ierr)

      implicit none

      integer*4 ipbm,idfdet,ierr
      
#ifdef FLUX_10A
#include "uhdef.fh"
#define FLXCOM  nusk_
#define NDBMNTPLID HID_NUFD2
#else
#ifdef FLUX_10C
#include "uhdef_10c.fh"
#define FLXCOM  nusk_
#define NDBMNTPLID HID_NUFD2
#else
#ifdef FLUX_11A
#include "uhdef_11a.fh"
#define FLXCOM  nusk_
#define NDBMNTPLID HID_NUFD2
#else
#include "beamntpl.h"
#define FLXCOM  fxvcsk_
#define NDBMNTPLID HID_NUFD1
#endif      
#endif      
#endif

#include "detectordefs.h"
#include "fxdetgeo.h"
#include "beamvectbl.h"

      integer*4 ipartbl(4)
      data      ipartbl/14,-14,12,-12/

      integer*4 nmaxnd,ifd

      INTEGER*4 LUNI
      PARAMETER (LUNI=10)

      integer*4 lrec
      save      lrec
      data      lrec/1024/

      integer*4 npt,maxnpt,nent,nptall
      save      npt,maxnpt,nent,nptall
      data      npt/1/
      data      maxnpt/0/

      integer*4 MAXFILES
      parameter (MAXFILES=100)
      integer*4 nmaxent(MAXFILES)

      integer*4 i,j
      real*8    prev

      real*4    fntotpau
      external  fntotpau
C--------------

      ierr=0
      if ((abs(ipbm).ne.12).and.(abs(ipbm).ne.14)) then
         write(*,*) 'Invalid PID ',ipbm,' was given.\n'
         stop
      endif
      if (maxnpt.eq.0) then
         nptall=0
         npt=1

         do 100 i=1,4
            maxfvnpt(i)  =0
            totnormtbl(i)=0.
            do 200 j=1,ITBLMAX
               idxtbl(i,j)   =0
               normfvtbl(i,j)=0.
 200        continue
 100     continue

 10      continue
            maxnpt=maxnpt+1
            call neopnufxv(luni,maxnpt,ierr)
            if (ierr.ne.0) then
               if (maxnpt.eq.1) then
                  write(*,*) 
     $                 'negtnufx:No file was specified in RFLIST'
                  stop
               else
                  maxnpt=maxnpt-1
                  ierr=0
                  goto 20
               endif
            endif
            
            call hnoent(3000,nmaxnd)
            if (nmaxnd.le.0) then
               write(*,*) 'negtnufx:Can not get detector defs.'
               write(*,*) 'negtnufx:Number of definitions=',nmaxnd
               write(*,*) 'fileno #',maxnpt
               stop
            else
               write(*,*) 'Number of detector defs:',nmaxnd
            endif
            
            do 60 nent=1,nmaxnd
               call hgntf(3000,nent,ierr)
               if (ierr.ne.0) then
                  write(*,*) 'negtnufx:Can not read detector defs. #'
     $                 ,nent
                  write(*,*) 'file #',maxnpt
                  stop
               endif

               do 70 ifd=1,nfd
                  if (npt.ne.1) then
                     if ((cbxfd(ifd).ne.bxfd(ifd)).or.
     $                    (cbyfd(ifd).ne.byfd(ifd)).or.
     $                    (cxfd(ifd).ne.xfd(ifd)).or.
     $                    (cyfd(ifd).ne.yfd(ifd)).or.
     $                    (czfd(ifd).ne.zfd(ifd)).or.
     $                    (chfd(ifd).ne.hfd(ifd)).or.
     $                    (cvfd(ifd).ne.vfd(ifd)).or.
     $                    (ncfd.ne.nfd)) then
                        write(*,*) 
     $                       'Inconsistent detector def. was found.'
                        stop
                     endif
                  else
                     ncfd      =nfd
                     cbxfd(ncfd)=bxfd(ncfd)
                     cbyfd(ncfd)=byfd(ncfd)
                     cxfd(ncfd)=xfd(ncfd)
                     cyfd(ncfd)=yfd(ncfd)
                     czfd(ncfd)=zfd(ncfd)
                     chfd(ncfd)=hfd(ncfd)
                     cvfd(ncfd)=vfd(ncfd)
                  endif
 70            continue
 60         continue

            call hnoent(NDBMNTPLID,nmaxent(maxnpt))
            if (ierr.ne.0) then
               write(*,*) 'negtnufx:Can not read 1st vector.'
               write(*,*) 'fileno #',maxnpt
               stop
            endif
            write(*,*) 'Number of entries in file #',
     $                 maxnpt,':',nmaxent(maxnpt)
            do 30 nent=1,nmaxent(maxnpt)
               call hgntf(NDBMNTPLID,nent,ierr)
               if (ierr.ne.0) then
                  write(*,*) 'negtnufx:Can not read vector #',nent
                  write(*,*) 'file #',maxnpt
                  stop
               else

                  if (idfdet.ne.idFD) then
                     goto 30
                  endif

                  nptall=nptall+1

                  if (nptall.gt.ITBLMAX) then
                     write(*,*) 'negtnufxv: Too much vectors.'
                     write(*,*) 'Check!'
                     stop
                  endif

                  Enutbl(nptall)     =EnuFD
                  ppitbl(nptall)     =ppiFD
                  Cospibmtbl(nptall) =CospibmFD
                  normtbl(nptall)    =normFD
                  ppi0tbl(nptall)    =ppi0FD
                  cospi0bmtbl(nptall)=cospi0bmFD
                  rnutbl(nptall)     =rFD
                  xnutbl(nptall)     =xnuFD
                  ynutbl(nptall)     =ynuFD
                  ppidtbl(nptall)    =ppidFD
                  modetbl(nptall)    =modeFD
                  idFDtbl(nptall)    =idFD
#ifdef FLUX_10A
                  psi0tbl(nptall)    =psi0FD
                  cossi0bmtbl(nptall)=cossi0bmFD
                                    
                  nvtx0tbl(nptall)   =nvtx0FD
                  spidtbl(nptall)    =spidFD
                  pgentbl(nptall)    =pgenFD
                  intgttbl(nptall)   =intgtFD
                  smechtbl(nptall)   =smechFD
                  smedtbl(nptall)    =smedFD
                  gppidtbl(nptall)   =gppidFD
                  pgpi0tbl(nptall)   =pgpi0FD
                  gpmechtbl(nptall)  =gpmechFD
                  gpmedtbl(nptall)   =gpmedFD
                  prmechtbl(nptall)  =prmechFD
                  prmedtbl(nptall)   =prmedFD
                  prdghttbl(nptall)  =prdghtFD
                  sdghttbl(nptall)   =sdghtFD
                  gpdghttbl(nptall)  =gpdghtFD
                  chaintbl(nptall)   =chainFD
                  giparttbl(nptall)  =gipartFD
                  gamom0tbl(nptall)  =gamom0FD
#endif
#if defined(FLUX_10C)||defined(FLUX_11A)
                  nvtx0tbl(nptall)   =nvtx0FD   
                  giparttbl(nptall)  =gipartFD  
                  gamom0tbl(nptall)  =gamom0FD  
                  EnuSKtbl(nptall)   =EnuSKFD   
                  normSKtbl(nptall)  =normSKFD  
                  anormtbl(nptall)   =anormFD   
                  ngtbl(nptall)      =ngFD      
      
                  DO 80 I=1,ngFD
                     gpxtbl(I,nptall)   = gpxFD(I)
                     gpytbl(I,nptall)   = gpyFD(I)
                     gpztbl(I,nptall)   = gpzFD(I)
                     gcosbmtbl(I,nptall)= gcosbmFD(I)
                     gvxtbl(I,nptall)   = gvxFD(I)
                     gvytbl(I,nptall)   = gvyFD(I)
                     gvztbl(I,nptall)   = gvzFD(I)
                     gpidtbl(I,nptall)  = gpidFD(I)
                     gmectbl(I,nptall)  = gmecFD(I)
 80               continue
#endif                  
C                  if (nptall.ge.2) then
C                     if (idFDtbl(nptall-1).ne.idFD) then
C                        write(*,*) "Different detector ID!"
C                        write(*,*) "Prev:",idFDtbl(nptall-1)
C                        write(*,*) "Now :",idFD
C                       stop
C                    endif
C                 endif
                  do 90 I=1,3
                     npitbl(i,nptall)     =npiFD(i)
                     xpitbl(i,nptall)     =xpiFD(i)
                     npi0tbl(i,nptall)    =npi0FD(i)
                     xpi0tbl(i,nptall)    =xpi0FD(i)
                     nnutbl(i,nptall)     =nnuFD(i)
#ifdef FLUX_10A
                     xsi0tbl(i,nptall)    =xsi0FD(i)
                     nsi0tbl(i,nptall)    =nsi0FD(i)
                     xsitbl(i,nptall)     =xsiFD(i)
                     prvtxtbl(i,nptall)   =prvtxFD(i)
                     xgpi0tbl(i,nptall)   =xgpi0FD(i)
                     xgpitbl(i,nptall)    =xgpiFD(i)
                     gpos0tbl(i,nptall)   =gpos0FD(i)
                     gvec0tbl(i,nptall)   =gvec0FD(i)
#endif
#ifdef FLUX_10C
                     gpos0tbl(i,nptall)   =gpos0FD(i)
                     gvec0tbl(i,nptall)   =gvec0FD(i)
#endif
#ifdef FLUX_11A
                     gmatFDtbl(I,nptall)   =gmatFD(I)   
                     gdistcFDtbl(I,nptall) =gdistcFD(I) 
                     gdistalFDtbl(I,nptall)=gdistalFD(I)
                     gdisttiFDtbl(I,nptall)=gdisttiFD(I)
                     gdistfeFDtbl(I,nptall)=gdistfeFD(I)
#endif                     

 90               continue
                  
                  j=int(modetbl(nptall)/10)
                  maxfvnpt(j)=maxfvnpt(j)+1
                  idxtbl(j,maxfvnpt(j))=nptall
                  if (maxfvnpt(j).eq.1) then
                     prev=0.D0
                  else
                     prev=normfvtbl(j,maxfvnpt(j)-1)
                  endif
                  
                  normfvtbl(j,maxfvnpt(j))
     $                 = prev
     $                  + normtbl(nptall)
     $                   *fntotpau(ipbm,Enutbl(nptall))
                  
                  totnormtbl(j)=normfvtbl(j,maxfvnpt(j))

               endif
 30         continue

            call hrend('nubeam')
            if (maxnpt.ge.100) goto 20
         goto 10
            
 20      continue

      endif
      
      return

      end
