C*******************************************************
C    <Updates>
C
C    2010.06 P. de Perio - Get filename from command line
C    2011.08 S. Tobayama - Flux input option implemented
C
C******************************************************    

      PROGRAM NEUT

      IMPLICIT NONE

C#include <f77/f77_floatingpoint.h>

#include "necardev.h"
#include "neutmodel.h"

      REAL*4    POS(3),DIR(3),PNEU(3)
      REAL*4    ENEU

      REAL*4    rndenuevtrt
      external  rndenuevtrt

      REAL*4    fntotpau
      external  fntotpau

      INTEGER*4 loadflx
      external  loadflx

      REAL*4    RLU
      external  RLU

      REAL*4    RNDTMP

      INTEGER*4 I,J,IDUM,IERR
      INTEGER*4 IRETRY

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL*4 COST,PHI

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#include "nework.h"
#include "vcwork.h"
#include "neutparams.h"      
#include "posinnuc.h"

      character*80 FILENAME,OPTION,TREE,VECBRANCH,VTXBRANCH
      data FILENAME/'neutvect.root'/
      data OPTION/'recreate'/
      data TREE/'neuttree'/
      data VECBRANCH/'vectorbranch'/
      data VTXBRANCH/'vertexbranch'/

      integer*4 rootopen,rootwrite,rootclose
      external rootopen,rootwrite,rootclose

      integer*4 neutfill,neutfilldsdq2
      external  neutfill,neutfilldsdq2

      integer*4   iret

      integer*4 lux,iseed1,iseed2,k1,k2,k3,k4

      LUX = 3

C--Read card files

      CALL NECARD
      CALL NECARDEV

C--Fill interactio model
      call nefillmodel

C--open file
      IF (IARGC() .GE. 2) THEN
         CALL GETARG(2, FILENAME)
      ENDIF

      iret = rootopen(FILENAME,OPTION)

      if (iret .ne. 0) then
         write(*,*) 'Failed to open file ',FILENAME,
     $              ' with option', OPTION ,'.'
         stop
      endif

C--Read flux histogram from ROOT file
      if (MPVEVCT.eq.3) then
         iret = loadflx(IDPTEVCT,FILENMEVCT,HISTNMEVCT,INMEVEVCT)
         if (iret .ne. 0) then
            write(*,*) 'Failed in reading flux histogram.'
            stop
         endif
      endif

C--Main loop
      
      DO 100 I=1,NEVTEVCT

C--Set vertex
         IF (MPOSEVCT.eq.1) then
            pos(1)=posevct(1)
            pos(2)=posevct(2)
            pos(3)=posevct(3)
         else if (MPOSEVCT.eq.2) then
            call nesetpos(RADEVCT,0,pos)
         else
            write(*,*) 'NEUT : MPOS is not set to 1 or 2.'
            stop
         endif

C--Set direction
 20      if (MDIREVCT.eq.1) then
            dir(1)=direvct(1)
            dir(2)=direvct(2)
            dir(3)=direvct(3)
         else if (MDIREVCT.eq.2) then
            call rndir(dir)
         else
            write(*,*) 'NEUT : MDIR is not set to 1 or 2.'
            stop
         endif

C--Set energy
         IRETRY = 0
 30      continue
         if (MPVEVCT.eq.1) then
            ENEU=PVEVCT(1)
         else if (MPVEVCT.eq.2) then
            ENEU=(PVEVCT(2)-PVEVCT(1))*RLU(IDUM)+PVEVCT(1)
         else if (MPVEVCT.eq.3) then
C--Draw ENEU from event rate (flux*cross-section) distribution
            RNDTMP=RLU(IDUM)
            ENEU = rndenuevtrt(RNDTMP)
         else
            write(*,*) 'NEUT : MPV is not set to 1, 2 or 3.'
            stop
         endif
         if (fntotpau(IDPTEVCT,ENEU/1000.).le.1.e-7) then
            IRETRY = IRETRY + 1
            if (IRETRY .gt. 100000) then
               write(*,*) 
     $              'Too many tries to set neutrino energy'
               write(*,*) 
     $              'Check the energy and interaction modes in the card'
               stop
            endif
            goto 30
         endif
         
         DO 50 J=1,3
            PNEU(J)=ENEU*DIR(J)
 50      continue
         
C--Generate event
         CALL RLUXAT(LUX,ISEED1,K1,K2)
         CALL RM48UT(ISEED2,K3,K4)
         WRITE(6,'(A4,I7,A7,I10,I10,I10,I10,I10,I10)')
     $        'Ev.#',I, ' SEEDS ',ISEED1,K1,K2,ISEED2,K3,K4

         call NEVECGEN(IDPTEVCT,POS,PNEU,IERR)

         if (ierr.ne.0) then
            goto 20
         endif

         iret = neutfill(filename,tree,vecbranch,vtxbranch)
         if (iret .ne. 0) then
            write(*,*) 'Failed in filling vector.'
            stop
         endif
C         iret = neutfilldsdq2(filename,tree,vecbranch,vtxbranch)

 100  continue

      iret = rootwrite(FILENAME)

      if (iret .ne. 0) then
         write(*,*) 'Failed to write file ',FILENAME
         stop
      endif

      iret = rootclose(FILENAME)

      if (iret .ne. 0) then
         write(*,*) 'Failed to close file ',FILENAME
         stop
      endif

      CALL NECHADIR(COST,PHI,DIR,POS)

      STOP
      END

