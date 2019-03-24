      subroutine neopskfxv(luni,chdir,npt,ierr)

      implicit none

#ifdef FLUX_10A
#include "uhdef.fh"
#else
#ifdef FLUX_10C
#include "uhdef_10c.fh"
#else
#ifdef FLUX_11A
#include "uhdef_11a.fh"
#else
#ifdef FLUX_11B
#include "uhdef_11b.fh"
#else
#ifdef FLUX_13
#include "uhdef_13_uwfunc.fh"
#else
#include "beamntpl.h"
#endif      
#endif      
#endif      
#endif
#endif      

      character*(*) chdir
      integer*4 luni,npt,ierr
C--------------
      
      
      character*160 fname,unit,lgname,access,stg,multi,
     &              comment,hname,fmode
      character*160 comnt,f_recl

      integer*4     lrec

      ierr=0

      call rgetv(luni,npt,fname,unit,lgname,
     &           access,stg,multi,comment,hname,fmode)
      if(fname.eq.' ')then
         write(*,*) 
     $        'neopskfxv:No file was specified in RFLIST'
         ierr=-1
         return
      endif

      comnt=comment
      call lcase(comnt)
      call getvar(comnt,'recl',f_recl)
      if(f_recl.ne.' ') then
         read(f_recl,*) lrec
      else
         lrec=5670
      endif

      call hropen(LUNI,'skbeam',fname,' PX',lrec,ierr)
      if (ierr.ne.0) then
         write(*,*) 'neopskfxv:file in RFLIST does not exist.'
         write(*,*) 'filename:',fname
         ierr=-1
         stop
      else
         write(*,*) 'Opened file : ', fname
      endif

      call hrin(2000,999,0)
      CALL HBNAME(2000,' ',0,'$CLEAR')
      CALL HBNAME(2000,'NU',EnuSK,'$SET')
      
      call hgnt(2000,1,ierr)

      return

      end
