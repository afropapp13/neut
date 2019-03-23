***********************************************************************
*     ------------------------
      SUBROUTINE nefillver
*     ------------------------
*     
*     ( purpose )
*       Store library version in "NEUTVERSION"
*
*     ( input ) 
*       None
*
*     ( output )
*       Common block "NEUTVERSION"
*
*     ( creation date and author )
*       2007.12.01 G.Mitsuka
*
***********************************************************************

#include"nefillver.h"

      integer corever, nucever, nuccver

C     Get version number of each library
      call nelibver(corever)
      call eflibver(nucever)
      call nrlibver(nuccver)

C     Store in common block "NEUTVERSION"
      corev = corever
      nucev = nucever
      nuccv = nuccver

      print*, ''
      print*, '======= Library version ======='
      print*, 'NEUTCORE  :',corev
      print*, 'NUCEFF    :',nucev
      print*, 'NUCCORSPL :',nuccv

      return
      end
