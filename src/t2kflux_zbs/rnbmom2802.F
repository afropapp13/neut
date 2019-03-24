      subroutine rnbmom2802(idfdet,iparfx,iparbm,iptr,pos,pmom)

      implicit none
      
C #include "beamntpl.h"
#include "beamvectbl.h"
#define IDFDET

      integer*4 idfdet,iparfx,iparbm,iptr
      real*4    pmom(3),pos(3)

      real*4    rlu,fntotpau
      external  rlu,fntotpau

      integer*4 norgvecmax
      integer*4 norgvec
      data norgvecmax/0/
      data norgvec/0/
      save norgvecmax,norgvec

      integer*4    numntvars
      data         numntvars/30/

C      integer*4    maxntvars
C      parameter    (maxntvars=30)

      real*8    rtry
      real*4    dum

      integer*4 ipartbl(4)
      data      ipartbl/14,-14,12,-12/
      
      integer*4 ipcond,i,ierr
      integer*4 moved,ipos

      call nerdnufxv(iparbm,idfdet,ierr)
      if (ierr.ne.0) then
         write(*,*) 'rnbmom280: Error in retriving flux info.'
         stop
      endif

      do 10 ipcond=1,4
         if (ipartbl(ipcond).eq.iparfx) goto 20
 10   continue

      write(*,*) 'rnbmom280: Invalid IPAR was given..'
      stop
      
 20   continue
      
      rtry=rlu(dum)*totnormtbl(ipcond)

      if (rtry.le.normfvtbl(ipcond,1)) then
         iptr = idxtbl(ipcond,1)
         goto 200
      endif
      if (rtry.ge.normfvtbl(ipcond,maxfvnpt(ipcond)-1)) then
         iptr = idxtbl(ipcond,maxfvnpt(ipcond))
         goto 200
      endif
         
      moved=maxfvnpt(ipcond)/2
      ipos =maxfvnpt(ipcond)-moved

 100  continue

C      write(*,'(A4,I8,A9,I8,A7,F10.8,A8,F10.8,A3,F10.8)') 
C     $     'ipos:',ipos,' / moved:',moved,' / try:',
C     $     rtry/totnormtbl(ipcond),' / norm:',
C     $     normfvtbl(ipcond,ipos-1)/totnormtbl(ipcond),' - ',
C     $     normfvtbl(ipcond,ipos)/totnormtbl(ipcond)

         if (ipos.eq.1) then
            write(*,*) 'rnbmom2802: Illegal condition ipos=1...'
            stop
         endif

         moved=moved/2
         if (moved.le.1) then
            moved=1
         endif
         if (rtry.le.normfvtbl(ipcond,ipos)) then
            if (rtry.ge.normfvtbl(ipcond,ipos-1)) then
               iptr = idxtbl(ipcond,ipos)
               goto 200
            endif
            ipos=ipos-moved
            if (ipos.le.1) then
               ipos=1
            endif
            goto 100
         endif

         if (rtry.ge.normfvtbl(ipcond,ipos)) then
            ipos=ipos+moved
            if (ipos.ge.maxfvnpt(ipcond)-1) then
               ipos=maxfvnpt(ipcond)-1
            endif
         endif
      goto 100
      
 200  continue

 50   continue
      
C      write(*,'(A4,I8,A9,I8,A7,F10.8,A8,F10.8,A3,F10.8)') 
C     $     'ipos:',ipos,' / moved:',moved,' / try:',
C     $     rtry/totnormtbl(ipcond),' / norm:',
C     $     normfvtbl(ipcond,ipos-1)/totnormtbl(ipcond),' - ',
C     $     normfvtbl(ipcond,ipos)/totnormtbl(ipcond)

      pos(1)=xnutbl(iptr)
      pos(2)=ynutbl(iptr)
      pos(3)=0.
      
      do 60 I=1,3
         pmom(I)=enutbl(iptr)*nnutbl(i,iptr)*1.e3
 60   continue
      
      if (idfdet.eq.0)  idfdet=idfdtbl(iptr)
      if (idfdet.ne.idfdtbl(iptr)) then
#ifndef IDFDET
         write(*,*) 
     $        "rnbmom280: Det.ID is not the one in the flux table"
         stop
#else
         goto 20
#endif         
      endif

      write(*,*) "Detector: ",idfdtbl(iptr)
      return

      end
