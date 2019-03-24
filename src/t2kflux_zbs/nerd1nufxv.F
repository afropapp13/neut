      SUBROUTINE nerd1nufxv(ipbm,idfdet,ierr)

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

      integer*4 nmaxnd,ifd

      INTEGER*4 LUNI
      PARAMETER (LUNI=10)

      integer*4 lrec
      save      lrec
      data      lrec/1024/

      integer*4 ndet,nev,maxnev,maxdet,nfile
      save      ndet,nev,maxnev,maxdet,nfile

      data      nev/1/
      data      ndet/1/
      data      nfile/0/

C--------------

      ierr=0
      if ((abs(ipbm).ne.12).and.(abs(ipbm).ne.14)) then
         write(*,*) 'Invalid PID ',ipbm,' was given.\n'
         stop
      endif

 10   if ((nfile.eq.0).or.(nev.gt.maxnev)) then

         if (nfile.ne.0) then
            call hrend('nubeam')
         endif

         nfile = nfile + 1

         call neopnufxv(luni,nfile,ierr)
         
         if (ierr.ne.0) then
            if (nfile.eq.1) then
               write(*,*) 
     $              'negtnufx:No file was specified in RFLIST'
               stop
            else
               ierr=1
               goto 20
            endif
         endif

         call hnoent(3000,maxdet)
         if (maxdet.le.0) then
            write(*,*) 'negtnufx:Can not get detector defs.'
            write(*,*) 'negtnufx:Number of definitions=',maxdet
            write(*,*) 'fileno #',nfile
            stop
         else
            write(*,*) 'Number of detector defs:',maxdet
         endif

         do 60 ndet=1,maxdet
            call hgntf(3000,ndet,ierr)
            if (ierr.ne.0) then
               write(*,*) 'negtnufx:Can not read detector defs. #'
     $              ,ndet
               write(*,*) 'file #',nfile
               stop
            endif
            
            do 70 ifd=1,nfd
               if (nfile.ne.1) then
                  if ((cbxfd(ifd).ne.bxfd(ifd)).or.
     $                 (cbyfd(ifd).ne.byfd(ifd)).or.
     $                 (cxfd(ifd).ne.xfd(ifd)).or.
     $                 (cyfd(ifd).ne.yfd(ifd)).or.
     $                 (czfd(ifd).ne.zfd(ifd)).or.
     $                 (chfd(ifd).ne.hfd(ifd)).or.
     $                 (cvfd(ifd).ne.vfd(ifd)).or.
     $                 (ncfd.ne.nfd)) then
                     write(*,*) 
     $                    'Inconsistent detector def. was found.'
                     stop
                  endif
               else
                  ncfd      =nfd
                  cbxfd(ifd)=bxfd(ifd)
                  cbyfd(ifd)=byfd(ifd)
                  cxfd(ifd)=xfd(ifd)
                  cyfd(ifd)=yfd(ifd)
                  czfd(ifd)=zfd(ifd)
                  chfd(ifd)=hfd(ifd)
                  cvfd(ifd)=vfd(ifd)
               endif
 70         continue
 60      continue
         
         call hnoent(NDBMNTPLID,maxnev)
         if (ierr.ne.0) then
            write(*,*) 'negtnufx:Can not read 1st vector.'
            write(*,*) 'fileno #',nfile
            stop
         endif
         write(*,*) 'Number of entries in file #',
     $        nfile,':',maxnev
         nev = 1
      endif

 30   call hgntf(NDBMNTPLID,nev,ierr)
      if (ierr.ne.0) then
         write(*,*) 'negtnufx:Can not read vector #',nev
         write(*,*) 'file #',nfile
         stop
      else
         nev = nev + 1
         if ((idfdet.ne.-1).and.(idfdet.ne.idFD)) then
            goto 10
         endif
      endif
         
      if (nev.gt.maxnev) then
         call hrend('nubeam')
      endif
            
 20   continue

      return

      end
