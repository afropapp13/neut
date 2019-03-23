      PROGRAM NEUT

      IMPLICIT NONE

C#include <f77/f77_floatingpoint.h>

#include "necard.h"
#include "necardev.h"

      REAL*4    POS(3),DIR(3),PNEU(3)
      REAL*4    ENEU

      REAL*4    RLU
      external  RLU

      INTEGER*4 I,J,IDUM,IERR

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL*4 COST,PHI

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      INTEGER*4 LUNO
      LUNO = 20

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C      CALL ieee_handler('set','common',SIGFPE_ABORT)
c      CALL ieee_handler('set','division',SIGFPE_ABORT)
C      CALL ieee_handler('set','invalid',SIGFPE_ABORT)

C--Reac card files

      CALL NECARD
      CALL NECARDEV

C--Initialize FILE
C      CALL KZINIT
C      CALL SKOPENF(LUNO,1,'Z',IERR)
C      IF (IERR.ne.0) THEN
C         WRITE(*,*) "CAN not open FILE(LUN=",LUNO,")"
C         STOP
C      ENDIF

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
         if (MPVEVCT.eq.1) then
            ENEU=PVEVCT(1)
         else if (MPVEVCT.eq.2) then
            ENEU=(PVEVCT(2)-PVEVCT(1))*RLU(IDUM)+PVEVCT(1)
         else
            write(*,*) 'NEUT : MPV is not set to 1 or 2.'
            stop
         endif
         
         DO 50 J=1,3
            PNEU(J)=ENEU*DIR(J)
 50      continue
         
C--Generate event
         call NEVECGEN(IDPTEVCT,POS,PNEU,IERR)

         if (ierr.ne.0) then
            goto 20
         endif

 100  continue

      CALL NECHADIR(COST,PHI,DIR,POS)

      STOP
      END
