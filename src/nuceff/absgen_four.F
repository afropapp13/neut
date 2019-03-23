
c ---------------------------------------------------------------

      subroutine absgen_four(s)

      implicit none

      real*4 s

C      real*4 EandP(4)
C      real*4 p1cm(4),p2cm(4),p3cm(4),p4cm(4)
C      real*4 p1lab(4),p2lab(4),p3lab(4),p4lab(4)
      real*4 mn/0.93893/

C      real*4 tecm,amass(18),pcm(5,18),wt

C      integer*4 lpdg,npdg(4)
C      integer*4 np,kgenev

      integer*4 iseq(4,24)/1,2,3,4,1,2,4,3,1,3,2,4,
     &                     1,3,4,2,1,4,2,3,1,4,3,2,
     &                     2,1,3,4,2,1,4,3,2,3,1,4,
     &                     2,3,4,1,2,4,1,3,2,4,3,1,
     &                     3,1,2,4,3,1,4,2,3,2,1,4,
     &                     3,2,4,1,3,4,1,2,3,4,2,1,
     &                     4,1,2,3,4,1,3,2,4,2,1,3,
     &                     4,2,3,1,4,3,1,2,4,3,2,1/

      logical ierr

#include "absneject.h"      
C      common /labEandP/ lpdg,EandP
C      common /cmnucleons/ p1cm,p2cm,p3cm,p4cm
C      common /labnucleons/ npdg,p1lab,p2lab,p3lab,p4lab

#include "absgeninout.h"
C      COMMON /NEGENIN / NP, TECM, AMASS, KGENEV
C      COMMON /NEGENOUT/ PCM, WT

      real*4    wtmax,xr
      integer*4 ixr
      integer*4 j1,j2,j3,j4

      np = 4
      kgenev = 0

      amass(1) = mn
      amass(2) = mn
      amass(3) = mn
      amass(4) = mn

      tecm = sqrt(s)

      wtmax = 0.19

 160  call absgenbod(ierr)
      if (ierr) then
         write (6,*) "error from absgenbod"
         return
      end if

      call ranlux(xr,1)
      xr = wtmax*xr
      if (wt.lt.xr) go to 160

c  Randomize to make sure there is no bias

      call ranlux(xr,1)
      ixr = 24.*xr + 1.
      j1  = iseq(1,ixr)
      j2  = iseq(2,ixr)
      j3  = iseq(3,ixr)
      j4  = iseq(4,ixr)

      p1cm(1) = pcm(4,j1)
      p1cm(2) = pcm(1,j1)
      p1cm(3) = pcm(2,j1)
      p1cm(4) = pcm(3,j1)

      p2cm(1) = pcm(4,j2)
      p2cm(2) = pcm(1,j2)
      p2cm(3) = pcm(2,j2)
      p2cm(4) = pcm(3,j2)     

      p3cm(1) = pcm(4,j3)
      p3cm(2) = pcm(1,j3)
      p3cm(3) = pcm(2,j3)
      p3cm(4) = pcm(3,j3)

      p4cm(1) = pcm(4,j4)
      p4cm(2) = pcm(1,j4)
      p4cm(3) = pcm(2,j4)
      p4cm(4) = pcm(3,j4)

      end
