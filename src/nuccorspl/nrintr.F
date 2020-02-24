**********************************************************************
*     ------------------
      SUBROUTINE nrintr
*     ------------------
*
*     ( purpose )
*       CONSIDER NUCLEAR EFFECT FOR EACH nucleon
*
*     ( input )
*       COMMON VCWORK
*
*     ( output )
*       COMMON VCWORK
*
*     ( creation date and author )
*        W.Gajewski
*        Y.Hayato
*
*        2010.05 P.de Perio
*            - Changed EFTRACE arguments NZIG -> NEL,NINEL
*            - elastic scatter: IFLGVC = 9
*
*        2010.09 P. de Perio
*            - Add flag for pion hadron production
*            - Make pion re-interaction similar to neffect code
*            - Assumes FSI vertex history vars (NVERT,NVCVERT) initialized
*              elsewhere (in NEFFECT or other calling program)
**********************************************************************
      IMPLICIT NONE

#include "nework.h"
#include "vcwork.h"
#include "vcvrtx.h"
#include "nrcard.h"
#include "nrnuclparam.h"
#include "neutparams.h"

      REAL*4    PIN(3),POUT(3)
C      integer*4 npartl,ipint,ipiori,looppi
C      real*4    pint,xint
C
C      COMMON /PIINTH/ NPARTL,IPINT(50),PINT(3,50),XINT(3,50),IPIORI(50),
C     &     LOOPPI

#include "efpiinth.h"

#include "posinnuc.h"

#include "nucleonfsihist.h"

C      integer*4 ino16
C      real*4    xin
C      COMMON /SAVXIN/INO16,XIN(3)

      integer*4 i,j

      integer*4 nvc2,invc

      real*4    pabs
      integer*4 kin,kfi,nabs,nzig,nel,ninel

      real*4    xpos(3)

C
C Check whether use or not
C
      if (NUCRESCAT.eq.0) then
         return
      endif

      if (ibound.ne.1) then
         return
      endif

C
C     ++ ADDED BY M.N ON 23-OCT-88
C
      NPARTL=0
C
C     Check Free or not
C
      IF  ((PVC(1,2)**2+PVC(2,2)**2+PVC(3,2)**2).lt.1.e-5) THEN
         return
      endif

!
!     reset nucleon fsi common block
!
      NFnvert = 0
      NFnstep = 0
      RW_step_i= 0

  100 CONTINUE

      NPARTL=0
      NVC2=NVC
C     from 3rd particle to end of vect
      DO 1000 I=3,NVC


CCCC ------ Skip unnecessary particles
         IF (ICRNVC(I).eq.0)  goto 1000

CCCC ------  for protons and neutrons   W. Gajewski 3-25-2001
C         IF (ABS(IPVC(I)).EQ.2212 .OR. IPVC(I).EQ.2112) THEN
         IF (IPVC(I).EQ.2212 .OR. IPVC(I).EQ.2112) THEN
C
C- copy interaction position
C
            DO 110 J=1,3
               XPOS(J)=POSNUC(J,I)
               if (XPOS(J).eq.-9999.) then
                  write(*,*) 'nrintr: Position for the particle #',
     $                 I,'(PID=',ipvc(i),') seems to be strange.'
                  stop
               endif
 110        continue

C            write(*,'(A,I3,A,I5,A,F12.7,A,F12.7,A,F12.7,A)')
C     $           'Initial:nvc=',nvc,':IP=',ipvc(I),
C     $           ':(',pvc(1,I),',',pvc(2,I),',',pvc(3,I),')'
             write(*,*) 'WGPROTON HAPPENING'
             write(*,*) "nvc = ",nvc
             write(*,*) "nvc2 = ",nvc2
             write(*,*) "I = ",I
            CALL wgproton(i,xpos,apn,nvc2)

C            write(*,'(A,I3,A,I5,A,F12.7,A,F12.7,A,F12.7,A)')
C     $           'Final  :nvc=',nvc2,':IP=',ipvc(I),
C     $           ':(',pvc(1,I),',',pvc(2,I),',',pvc(3,I),')'

            if (nvc.ne.nvc2) then
               do 950 invc=nvc+1,nvc2
C                  write(*,'(A8,I7,A,I5,A,F12.7,A,F12.7,A,F12.7,A)')
C     $                 'Final  :',invc,':IP=',ipvc(invc),
C     $                 ':(',pvc(1,invc),',',pvc(2,invc),
C     $                 ',',pvc(3,invc),')'
 950           continue
            endif
         END IF

 1000 CONTINUE
C
C     ++ ADDED BY W. Gajewski (3-29-2001)
C     -- RE-CONSIDER PION INTERACTION FOR PIONS PRODUCED IN N-N INT.
C     IN 16O
C
      IF (nvc2.ge.MAXVC) then
         nvc2=MAXVC-1
         write(*,*) 'nrintr: Too many particles...'
      endif
      do i=nvc+1,nvc2
CCCC ------- for PION ----------
         IF (ABS(IPVC(I)).EQ.211 .OR. IPVC(I).EQ.111) THEN
            CALL UCOPY(PVC(1,I),PIN,3)
            PABS=SQRT(PIN(1)**2+PIN(2)**2+PIN(3)**2)

            KIN=IPVC(I)
            LOOPPI=I

            DO 120 J=1,3
               XPOS(J)=POSNUC(J,I)
               if (XPOS(J).eq.-9999.) then
                  write(*,*) 'nrintr: 120:Position for the particle #',
     $                 I,'(PID=',ipvc(i),') seems to be strange.'
                  stop
               endif
 120        continue

            CALL EFTRACE(APN,AP,KIN,XPOS,PIN,NABS,NEL,NINEL,KFI,POUT)
            NZIG = NEL+NINEL
C
C --        ABSORPTION
C
C            IF (KFI.EQ.100000) THEN
            IF (NABS.EQ.0) THEN
               IF (NPARTL.GT.0) THEN
                  IFLGVC(I)=7
               else
                  IFLGVC(I)=3
               endif
               ICRNVC(I)=0
               CALL UCOPY(POSVC,POSFVC(1,I),3)
               IVTFVC(I)=1
C
C --        CHARGE EXCHANGE and scattering
C
C            ELSE IF (KIN.NE.KFI) THEN
            ELSE IF ((KIN.NE.KFI).or.(NZIG.gt.0)) THEN
               IF (KIN.NE.KFI) THEN
                  IFLGVC(I)=4
               ELSE
                  IF (NPARTL.GT.0) THEN
                     IFLGVC(I)=7
                  ELSE IF (NINEL.gt.0) THEN
                     IFLGVC(I)=8
                  ELSE
                     IFLGVC(I)=9
                  ENDIF
               ENDIF
               ICRNVC(I)=0
               CALL UCOPY(POSVC,POSFVC(1,I),3)
               IVTFVC(I)=1
               IF (NVC2.lt.MAXVC) THEN
                  NVC2=NVC2+1
                  IPVC(NVC2)=KFI
                  CALL UCOPY(POUT,PVC(1,NVC2),3)
                  IORGVC(NVC2)=I
                  IFLGVC(NVC2)=0
                  ICRNVC(NVC2)=1
                  CALL UCOPY(POSVC,POSIVC(1,NVC2),3)
                  IVTIVC(NVC2)=1
                  CALL MCMASS(IPVC(NVC2),AMASVC(NVC2))
                  DO J=1,3
                     POSNUC(J,NVC2)=9999.
                  enddo
               ENDIF
C
C --        OTHER CASE (Momentum can change if Coulomb correction is on)
C
            ELSE IF (FEFCOUL.eq.0) then
C               CALL UCOPY(POUT,PVC(1,I),3)
               DO J=1,3
                  if (POUT(J).ne.PVC(J,I)) then
                     write(*,*)
     $                    'NRINTR: Momentum changed without interaction'
                     write(*,*)
     $                    'NRINTR: Mode=',modene
                     stop
                  endif
               enddo
            ENDIF
         else
            DO J=1,3
               POSNUC(J,I)=9999.
            enddo
         ENDIF
      end do

      NVC=NVC2

#ifdef SUBTRACT
      do 200 i=1,nvc
         IF (ICRNVC(I).eq.0)  goto 200
         if (((IPVC(I).eq.2112)).or.((IPVC(I).eq.2212))) then
            if (ipvc(I).eq.2112) then
               AMASS=939.566
            endif
            if (ipvc(I).eq.2212) then
               AMASS=938.272
            endif
            P2OLD=(PVC(1,i)**2+PVC(2,i)**2+pvc(3,i)**2)
            EOLD=sqrt(P2OLD+AMASS**2)
            write(*,'(A,I3,A,I5,A,F12.7,A,F12.7,A,F12.7,A,F12.7,A)')
     $           'Initial:nvc=',i,':IP=',ipvc(I),
     $           ':E=',EOLD,
     $           ':(',pvc(1,I),',',pvc(2,I),',',pvc(3,I),')'
            P2NEW=(EOLD-25.)**2-AMASS**2
            if (p2new.le.0) then
               p2new=0
            endif
            do 300 JJ=1,3
               pvc(jj,i)=sqrt(p2new/p2old)*pvc(jj,i)
 300        continue
            write(*,'(A,I3,A,I5,A,F12.7,A,F12.7,A,F12.7,A,F12.7,A)')
     $           'Corr.  :nvc=',i,':IP=',ipvc(I),
     $           'E=',EOLD-25.,
     $           ':(',pvc(1,I),',',pvc(2,I),',',pvc(3,I),')'
         endif
 200  continue
#endif
C
C
      RETURN
      END


