      PROGRAM T2KNEUT

      IMPLICIT NONE

#include "nework.h"
#include "mcgenpar.h"      
#include "necard.h"
#include "neutmodel.h"
#include "necardbm.h"
#include "nefillver.h"
#include "vcwork.h"
#include "posinnuc.h"

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

      REAL*4 POSI(3)
      REAL*4 DIR(3),GDIR(3)
      REAL*4 COST,PHI

      INTEGER*4  ISIZE
      PARAMETER (ISIZE=5000000)

      REAL*4     H(ISIZE)
      COMMON/PAWC/H

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER*4 LUNO

      INTEGER*4 I,J,IERR

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
         call nerd1nufxv(idbmflx,idbmdet,ierr)

         if (ierr.ne.0) goto 110

C --  copy ntuple variables to common 

         pos(1)=xnufd
         pos(2)=ynufd
         pos(3)=0

         do 60 J=1,3
            pmom(J)=enufd*nnufd(J)*1.e3
 60      continue

         idfxdet = idFD
         bmcnt(1)=bxfd(idfd)
         bmcnt(2)=byfd(idfd)
         fxdetcnt(1)=xfd(idfd)
         fxdetcnt(2)=yfd(idfd)
         fxdetcnt(3)=zfd(idfd)
         fxdetsz(1)=hfd(idfd)
         fxdetsz(2)=vfd(idfd)

         if (sqrt(pos(1)**2+pos(2)**2+pos(3)**2).gt.BMRADMX) then
            goto 10
         endif

C --  call vector generation func.

         if (enufd.le.0.100) then
            call NECLPOS
            call NECLRWRK
            call necpnewk(pos)
            IBOUND = 0
            NUMNE=1
            IPNE(1)=idbmpid
            PNE(1,1)=PMOM(1)/1000.
            PNE(2,1)=PMOM(2)/1000.
            PNE(3,1)=PMOM(3)/1000.
            IORGNE(1)=0
            IFLGNE(1)=-1
            ICRNNE(1)=0
            goto 50
         endif

 150     CALL nevect(idbmpid,pos,pmom,ierr)
         
         if (ierr.ne.0) then
            call kzeclr
            write(*,*) 'FAILED IN GENERATING EVENT #',I
            write(*,*) '                     IP     ',IDBMPID
            write(*,*) '                     PMOM(',
     $           PMOM(1),',',PMOM(2),',',PMOM(3),')'
            GOTO 150
         endif

C
C --  MAKE NEUT BANK
C
 50      CALL NEMKNEBK(POS)
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
         call nemkjndfx
         call vcmkdetgm
         call nemknetarg

         CALL KZWRIT(LUNO)
         
 100  continue
         
 110  CALL KZEND(LUNO)
      close(LUNO)

      STOP
      END

      SUBROUTINE RMARIN(int i1,int i2,int i3)

      write(*,*) "RMARIN: Dummy routine"
      return
      end
      
