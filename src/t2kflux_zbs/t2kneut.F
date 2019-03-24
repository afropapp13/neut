      PROGRAM T2KNEUT

      IMPLICIT NONE

#include "nework.h"
#include "mcgenpar.h"      
#include "necard.h"
#include "neutmodel.h"
#include "necardbm.h"
#include "nefillver.h"
#include "beamntpl.h"

      REAL*4 POSI(3)
      REAL*4 DIR(3),GDIR(3)
      REAL*4 COST,PHI

      INTEGER*4  ISIZE
      PARAMETER (ISIZE=5000000)

      REAL*4     H(ISIZE)
      COMMON/PAWC/H

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER*4 LUNO

      INTEGER*4 I,IERR

      INTEGER*4 iptr
      REAL*4    POS(3),PMOM(3)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      LUNO = 20

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! load configuration for interactions
      CALL NECARD
      CALL NECARDBM

C--Fill interaction models

      call nefillmodel
      
      corev = NECORE_VERSION
      nucev = NENUCE_VERSION
      nuccv = NENUCC_VERSION
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC            

C--Initialize FILE
      CALL KZINIT
      call hlimit(-1*ISIZE)
      CALL SKOPENF(LUNO,1,'Z',IERR)
      IF (IERR.ne.0) THEN
         WRITE(*,*) "CAN not open FILE(LUN=",LUNO,")"
         STOP
      ENDIF

      DO 100 I=1,NEVBM
         write(*,*) 'Event #',I
         CALL KZECLR
         
 10      continue
         call rnbmom2802(idbmdet,idbmflx,idbmpid,iptr,pos,pmom)
         if (sqrt(pos(1)**2+pos(2)**2+pos(3)**2).gt.BMRADMX) then
            goto 10
         endif

         CALL nevect(idbmpid,pos,pmom,ierr)
         
         if (ierr.ne.0) then
            call kzeclr
            write(*,*) 'FAILED IN GENERATING EVENT #',I
            write(*,*) '                     IP     ',IDBMPID
            write(*,*) '                     PMOM(',
     $           PMOM(1),',',PMOM(2),',',PMOM(3),')'
            GOTO 100
         endif

C
C --  MAKE NEUT BANK
C
         CALL NEMKNEBK(POS)
         CALL NEMKMODELBK
         CALL NEMKFSIBK
         CALL NEMKCRSBK
C     
C --  MAKE VECT BANK
C
         CALL VCMKVC
         CALL VCMKVX
C     --- MAKE NUCPOS BANK
         CALL NEMKNPOS
C      
         if (iptr.ne.0) then
            call necpjndfxcm(iptr)
            call nemkjndfx
         endif
         if (idbmdet.ne.0) then
            call vcmkdetgm
         endif

         call nemknetarg

         CALL KZWRIT(LUNO)
         
 100  continue
         
      CALL KZEND(LUNO)
      close(LUNO)

      STOP
      END

      SUBROUTINE RMARIN(int i1,int i2,int i3)

      write(*,*) "RMARIN: Dummy routine"
      return
      end
      
