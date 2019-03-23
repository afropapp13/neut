      PROGRAM chkfort

      implicit none

#include "nework.h"
#include "vcwork.h"
#include "vcvrtx.h"
#include "neutparams.h"      
#include "posinnuc.h"

      character*80 FILENAME,OPTION
      data FILENAME/'chkfort.root'/
      data OPTION/'recreate'/

      integer*4 rootopen,rootclose
      external rootopen,rootclose

      integer*4   iret

      iret = rootopen(FILENAME,OPTION)

      if (iret .ne. 0) then
         write(*,*) 'Failed to open file ',FILENAME,
     $              ' with option', OPTION ,'.'
         stop
      endif
         
      iret = rootclose(FILENAME)

      if (iret .ne. 0) then
         write(*,*) 'Failed to close file ',FILENAME, '.'
         stop
      endif
      
      write(*,*) 'succeed.'
      stop
      end
      
