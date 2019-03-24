      subroutine neopnufxv(luni,npt,ierr)

      implicit none

#ifdef FLUX_10A
#include "uhdef.fh"
#define NDBMNTPLID 3002
#else
#ifdef FLUX_10C
#include "uhdef_10c.fh"
#define NDBMNTPLID 3002
#else
#ifdef FLUX_11A
#include "uhdef_11a.fh"
#define NDBMNTPLID 3002
#else
#ifdef FLUX_11B
#include "uhdef_11b.fh"
#define NDBMNTPLID 3002
#else
#include "beamntpl.h"
#define NDBMNTPLID 3001
#endif      
#endif      
#endif      
#endif

#include "beamvectbl.h"

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
     $        'neopnufxv:No file was specified in RFLIST'
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

      call hropen(LUNI,'nubeam',fname,' PX',lrec,ierr)
      if (ierr.ne.0) then
         write(*,*) 'neopnufxv:file in RFLIST does not exist.'
         write(*,*) 'filename:',fname
         ierr=-1
         stop
      endif

      call hrin(NDBMNTPLID,999,0)
      CALL HBNAME(NDBMNTPLID,' ',0,'$CLEAR')
      CALL HBNAME(NDBMNTPLID,'NU',EnuFD,'$SET')
      
      call hrin(3000,999,0)
      CALL HBNAME(3000,' ',0,'$CLEAR')
      CALL HBNAME(3000,'DEFINE',NFD,'$SET')

      call hgnt(3000,1,ierr)
      call hgnt(NDBMNTPLID,1,ierr)


      return

      end
