      subroutine NRSETXPOS

      implicit none

#include "nrnuclparam.h"

      integer*4 ino16
      real*4    xin
      COMMON /SAVXIN/INO16,XIN(3)

      integer*4 ndata,lbank

      CHARACTER*32 CNAME
      INTEGER IDATA(20)
      REAL*4  RDATA(20)
      EQUIVALENCE (IDATA(1),RDATA(1))

      real*4  X(3)

C
C     -- SET INTERACTION POSITION
C

      CNAME='NUCPOS'
      CALL KZBLOC(CNAME,LBANK)
      IF (LBANK.eq.0) then
         write(*,*) '*** Can not locate NUCPOS ***'
         stop
      endif
      CALL KZGET1(CNAME,1,NDATA,IDATA)
      IF (NDATA.ne.3) then
         write(*,*) '*** NUCPOT: Length is incorrect ***'
         stop
      endif
      xin(1)=rdata(1)
      xin(2)=rdata(2)
      xin(3)=rdata(3)
      if ((XIN(1)**2+XIN(2)**2+XIN(3)**2).gt.9999.) then
         write(*,*) '*** NUCPOT: Location is incorrect ***'
         CALL EFRANBLL(APN,X)
         xin(1)=x(1)
         xin(2)=x(2)
         xin(3)=x(3)
      endif

      write(*,*)
     $     'Position:(',xin(1),',',xin(2),',',xin(3),')'

      return

      end
