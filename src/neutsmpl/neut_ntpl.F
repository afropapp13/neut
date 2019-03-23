      PROGRAM NEUT

      IMPLICIT NONE

C#include <f77/f77_floatingpoint.h>

#include "necardev.h"
#include "nework.h"      
#include "vcwork.h"
#include "vcvrtx.h"

      REAL*4    POS(3),DIR(3),PNEU(3)
      REAL*4    ENEU

      REAL*4    RLU
      external  RLU

      REAL*4    fntotpau
      external  fntotpau

      INTEGER*4 I,J,IDUM,IERR,N1,N2
      INTEGER*4 IRETRY

      character*127 out_fname

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL*4 COST,PHI

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      INTEGER*4  ICYCLE
      INTEGER*4  LRECL
      PARAMETER (LRECL=1024)

      INTEGER*4  ISIZE
      PARAMETER (ISIZE=5000000)

      REAL*4     H(ISIZE)
      COMMON/PAWC/H

C     Neut PARAMETER COMMON
      INTEGER*4   IEVENT,MODE,NPARNEUT,IPNEUT(100)
      REAL*4      ABSPNEUT(100),PNEUT(3,100)
      COMMON/COM1/IEVENT,MODE,NPARNEUT,
     $            IPNEUT,ABSPNEUT,PNEUT

C
C     Vector Parameter COMMON
C
      INTEGER*4   NPARVEC,IORGVEC(100),IPVEC(100),ICHVEC(100),
     $            IFLVEC(100)
      REAL*4      ABSPVEC(100),PVEC(3,100),POSV(3)
      COMMON/COM2/NPARVEC,IPVEC,IORGVEC,ICHVEC,IFLVEC,
     $            ABSPVEC,PVEC,POSV

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C-- Initialize ntuple
      
      CALL HLIMIT(ISIZE)

      out_fname='sample.nt'
      call getarg(2,out_fname)

      CALL HROPEN(85,'CWNT',out_fname,'N',LRECL,IERR)

      IF (IERR.ne.0) THEN
         WRITE(6,*) 'HROPEN:FAILED TO OPEN FILE ', out_fname
         STOP
      ENDIF

      CALL HBNT(10,'NEUTVEC',' ')
      CALL HBNAME(10,'NEUTINFO',IEVENT,
     $'NEV,MODE,NUMNU[0,100],IPNU(NUMNU),ABSPNU(NUMNU),PNU(3,NUMNU)')

      CALL HBNAME(10,'VECINFO',NPARVEC,
     $'NPAR[0,100],IPV(NPAR),IORGV(NPAR),ICRNV(NPAR),IFLGV(NPAR),
     $ABSPV(NPAR),PMOMV(3,NPAR),POS(3)')

C--Reac card files

      CALL NECARD
      CALL NECARDEV

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
130      continue
         if (MPVEVCT.eq.1) then
            ENEU=PVEVCT(1)
         else if (MPVEVCT.eq.2) then
            ENEU=(PVEVCT(2)-PVEVCT(1))*RLU(IDUM)+PVEVCT(1)
         else
            write(*,*) 'NEUT : MPV is not set to 1 or 2.'
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
            goto 130
         endif
         
         DO 50 J=1,3
            PNEU(J)=ENEU*DIR(J)
 50      continue
         
C--Generate event
         call NEVECGEN(IDPTEVCT,POS,PNEU,IERR)

         if (ierr.ne.0) then
            goto 20
         endif

         DO 60 N1=1,50
            IPNEUT(N1)=0
            ABSPNEUT(N1)=0.
            IPVEC(N1)=0
            IORGVEC(N1)=0
            ICHVEC(N1)=0
            ABSPVEC(N1)=0.
            IFLVEC(N1)=0

            DO 70 N2=1,3
               PNEUT(N2,N1)=0.
               PVEC(N2,N1) =0.
 70         CONTINUE
 60      CONTINUE

         IEVENT=I
         MODE=MODENE         
         NPARNEUT = NUMNE
C         write(*,*) "MODE=",MODE," / NPAR=",NPAR

         DO 30 N1=1,NUMNE
            IPNEUT(N1)=IPNE(N1)
            ABSPNEUT(N1)=sqrt( PNE(1,N1)**2
     $                        +PNE(2,N1)**2
     $                        +PNE(3,N1)**2)
            
            PNEUT(1,N1)=PNE(1,N1)
            PNEUT(2,N1)=PNE(2,N1)
            PNEUT(3,N1)=PNE(3,N1)
            
C    write(*,*) "IP(",N,")=",IPNE(N),
C    $                 " / ABSP(",N,")=",ABSP(N),
C    $                 " / P(",N,")=(",PNE(1,N),",",PNE(2,N),",",
C    $                 PNE(3,N),")"
 30      CONTINUE

         NPARVEC=NVC

         DO 40 N1=1,NVC
            IPVEC(N1)    = IPVC(N1)
            ICHVEC(N1)   = ICRNVC(N1)
            IORGVEC(N1)  = IORGVC(N1)
            IFLVEC(N1)   = IFLGVC(N1)
            
            ABSPVEC(N1)  = sqrt( PVC(1,N1)**2
     $                          +PVC(2,N1)**2
     $                          +PVC(3,N1)**2)
            
            PVEC(1,N1)   = PVC(1,N1)
            PVEC(2,N1)   = PVC(2,N1)
            PVEC(3,N1)   = PVC(3,N1)

 40      CONTINUE

         posv(1)    = pvtxvc(1,1)
         posv(2)    = pvtxvc(2,1)
         posv(3)    = pvtxvc(3,1)

         CALL HFNT(10)

 100  continue

      CALL HROUT(10,ICYCLE,' ')
      CALL HREND('CWNT')

      CALL NECHADIR(COST,PHI,DIR,POS)

      STOP
      END

      
