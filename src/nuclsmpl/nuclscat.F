      PROGRAM NEUT

      IMPLICIT NONE

C#include <f77/f77_floatingpoint.h>

#include "necard.h"
#include "necardev.h"
#include "nework.h"
#include "vcwork.h"
#include "vcvrtx.h"
#include "posinnuc.h"
#include "efpion.h"
#include "nucleonfsihist.h"

      REAL*4    POS(3),DIR(3),PNEU(3)
      REAL*4    ENEU

      REAL*4    RLU
      external  RLU

      DOUBLE PRECISION EVPIPROB
      external EVPIPROB

      INTEGER*4 I,J,IDUM,IERR,N1,N2,IVERT,IVCVERT

      character*127 out_fname,strpid

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL*4 COST,PHI

      REAL*4 RMAX,RAD,PTOT

      INTEGER*4 IDNUCLEON

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      INTEGER*4  ICYCLE
      INTEGER*4  LRECL
      PARAMETER (LRECL=1024)

      INTEGER*4  ISIZE
      PARAMETER (ISIZE=5000000)

      REAL*4     H(ISIZE)
      COMMON/PAWC/H

C     Neut PARAMETER COMMON
      INTEGER*4   IEVENT,MODE,NPARNEUT,IPNEUT(MAXNE)
      REAL*4      ABSPNEUT(MAXNE),PNEUT(3,MAXNE)
      COMMON/COM1/IEVENT,MODE,NPARNEUT,
     $            IPNEUT,ABSPNEUT,PNEUT

C
C     Vector Parameter COMMON
C
      INTEGER*4   NPARVEC,IORGVEC(MAXVC),IPVEC(MAXVC),ICHVEC(MAXVC),
     $            IFLVEC(MAXVC)
      REAL*4      ABSPVEC(MAXVC),PVEC(3,MAXVC),POSV(3)
      COMMON/COM2/NPARVEC,IPVEC,IORGVEC,ICHVEC,IFLVEC,
     $            ABSPVEC,PVEC,POSV

C
C     FSI history COMMON
C
        integer*4   oNFnvert, oNFiflag,oNFfirststep
        real*4      oNFx,oNFy,oNFz,oNFpx,oNFpy,oNFpz,oNFe

        common /COM3/ oNFnvert, oNFiflag(maxnucleonvert),
     &       oNFx(maxnucleonvert),oNFy(maxnucleonvert),
     $       oNFz(maxnucleonvert),
     &       oNFpx(maxnucleonvert),oNFpy(maxnucleonvert),
     $       oNFpz(maxnucleonvert),
     &       oNFe(maxnucleonvert),oNFfirststep(maxnucleonvert)

        integer*4   oNFnstep
        real*4      oNFecms2,oNFptot

        COMMON /COM4/oNFnstep, oNFecms2(maxnucleonstep),
     $       oNFptot(maxnucleonstep)

C
C     Nucleus Target COMMON
C
      INTEGER*4   oNUMBNDN,oNUMBNDP,oNUMFREP,oNUMATOM
      COMMON /COM5/ oNUMBNDN,oNUMBNDP,oNUMFREP,oNUMATOM

      INTEGER*4 IQUEST
      COMMON/QUEST/IQUEST(100)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C-- Initialize ntuple

      CALL HLIMIT(ISIZE)

      out_fname='sample.hb'
      call getarg(2,out_fname)

      call getarg(3,STRPID)
      READ(STRPID, FMT='(I5)') IDNUCLEON
C      IDPION = ATOI(STRPID)
      write(*,*) 'Trying to track PID: ',IDNUCLEON

      Iquest(10) = 256000
      CALL HROPEN(85,'CWNT',out_fname,'NQE',LRECL,IERR)

      IF (IERR.ne.0) THEN
         WRITE(6,*) 'HROPEN:FAILED TO OPEN FILE ', out_fname
         STOP
      ENDIF

      CALL HBNT(1,'NEUTVEC',' ')
      CALL HBNAME(1,'VECINFO',NPARVEC,
     $'NPVC[0,100],IPVC(NPVC),IORGVC(NPVC),ICHVC(NPVC),IFLVC(NPVC),
     $ABSPVC(NPVC),PVC(3,NPVC),POSVC(3)')

      CALL HBNAME(1,'NFSIVERT',oNFnvert,
     $     'NFNVERT[0,200],NFIFLAG(NFNVERT),
     &      X(NFNVERT),Y(NFNVERT),Z(NFNVERT),
     &      PX(NFNVERT),PY(NFNVERT),PZ(NFNVERT),
     &      E(NFNVERT),NFfirststep(NFNVERT)')


      CALL HBNAME(1,'NFSIPART',oNFnstep,
     &     'NFNSTEP[0,2000], ECMS2(NFNSTEP),PROB(NFNSTEP)')

      CALL HBNAME(1,'NETARG',oNUMBNDN,
     $' NUMBNDN, NUMBNDP, NUMFREP, NUMATOM')

C--Read card files

      CALL NECARD
      CALL NECARDEV

C--Main loop

C      RMAX = C * 2.5
      RMAX = C * 2.

      DO 100 I=1,NEVTEVCT

C      if ( MOD(I,1000).EQ.0 ) WRITE(*,*) 'Processing event: ',I

C--Set vertex
 110     continue
            posv(2)=-RMAX+RLU(IDUM)*RMAX*2.
            posv(3)=-RMAX+RLU(IDUM)*RMAX*2.
 2       RAD   =SQRT(POSV(2)**2+POSV(3)**2)
         IF(RAD.GE.RMAX)GO TO 110
         POSV(1)=-SQRT(RMAX**2-RAD**2)+0.0001

C--Set direction
         dir(1)=1.
         dir(2)=0.
         dir(3)=0.

C--Set energy
         if (MPVEVCT.eq.1) then
            PTOT=PVEVCT(1)
         else if (MPVEVCT.eq.2) then
            PTOT=(PVEVCT(2)-PVEVCT(1))*RLU(IDUM)+PVEVCT(1)
         else
            write(*,*) 'NEUT : MPV is not set to 1 or 2.'
            stop
         endif


         CALL VCCLCM
         NVC      = 4
         IPVC(4)  = IDNUCLEON
         call mcmass(IPVC(4),AMASVC(4))
         PVC(1,4) = PTOT
         PVC(2,4) = 0.
         PVC(3,4) = 0.
         POSNUC(1,4)=posv(1)
         POSNUC(2,4)=posv(2)
         POSNUC(3,4)=posv(3)

         IORGVC(4)= -1
         IFLGVC(4)= 0
         ICRNVC(4)= 1
         TIMVC(4) = 0

C        -- dummy parameter to pass the check
         IPVC(2)  = IDNUCLEON
         call mcmass(IPVC(2),AMASVC(2))
         PVC(1,2) = PTOT
         PVC(2,2) = 0.
         PVC(3,2) = 0.
         IORGVC(2)= -1
         IFLGVC(2)= 0
         ICRNVC(2)= 0
         TIMVC(2) = 0

         POSNUC(1,2)=posv(1)
         POSNUC(2,2)=posv(2)
         POSNUC(3,2)=posv(3)

         IBOUND   = 1

C-- Generate event
         CALL EFCLFSI
C         CALL NEFFECT
         CALL NRINTR

C --  Fill VCWORK
         NPARVEC=NVC

         DO 40 N1=2,NVC
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

C --   FILL FSI history
         oNFNVERT=NFnvert

         DO 50 IVERT=1,NFnvert
            oNFiflag(IVERT)     =NFiflag(IVERT)
            oNFx(IVERT)         =NFx(IVERT)
            oNFy(IVERT)         =NFy(IVERT)
            oNFz(IVERT)         =NFz(IVERT)
            oNFpx(IVERT)        =NFpx(IVERT)
            oNFpy(IVERT)        =NFpy(IVERT)
            oNFpz(IVERT)        =NFpz(IVERT)
            oNFe(IVERT)         =NFe(IVERT)
            oNFfirststep(IVERT) =NFfirststep(IVERT)
 50      CONTINUE

         oNFnstep = NFnstep

         DO 60 IVCVERT=1,NFnstep
            oNFecms2(IVCVERT) = NFecms2(IVCVERT)
            oNFptot(IVCVERT)  = NFptot(IVCVERT)
 60      CONTINUE

         oNUMBNDN = NUMBNDN
         oNUMBNDP = NUMBNDP
         oNUMFREP = NUMFREP
         oNUMATOM = NUMATOM

         CALL HFNT(1)

 100  continue

      CALL HROUT(1,ICYCLE,' ')
      CALL HREND('CWNT')

C      CALL NECHADIR(COST,PHI,DIR,POSV)

      STOP
      END
