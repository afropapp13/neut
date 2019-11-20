      PROGRAM NEUT

      IMPLICIT NONE

C#include <f77/f77_floatingpoint.h>

#include "necard.h"
#include "necardev.h"
#include "nework.h"
#include "vcwork.h"
#include "vcvrtx.h"
#include "efpion.h"
#include "posinnuc.h"
#include "fsihist.h"

      REAL*4    POS(3),DIR(3),PNEU(3)
      REAL*4    ENEU

      REAL*4    RLU
      external  RLU

      REAL*4 EVPIPROB
      external EVPIPROB

      INTEGER*4 I,J,IDUM,IERR,N1,N2,IVERT,IVCVERT

      character*127 out_fname,strpid

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL*4 COST,PHI

      REAL*4 RMAX,RAD,PTOT

      INTEGER*4 IDPION

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
      INTEGER*4   oNVERT, oNVCVERT,oIFLGVERT(MAXVERT),
     $     oIPVERT(MAXVCVERT), oIVERTI(MAXVCVERT), oIVERTF(MAXVCVERT)
      REAL*4  oPOSVERT(3,MAXVERT), oDIRVERT(3,MAXVCVERT),
     $        oABSPVERT(MAXVCVERT), oABSTPVERT(MAXVCVERT),oFSIPROB

      COMMON /COM3/ oNVERT,oPOSVERT,oIFLGVERT
      COMMON /COM4/ oNVCVERT, oDIRVERT,
     $    oABSPVERT, oABSTPVERT,
     $    oIPVERT,oIVERTI,oIVERTF,oFSIPROB


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
      READ(STRPID, FMT='(I5)') IDPION
C      IDPION = ATOI(STRPID)
      write(*,*) 'Trying to track PID: ',IDPION

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

      CALL HBNAME(1,'FSIVERT',oNVERT,
     $'NVERT[0,100],POSVERT(3,NVERT),IFLGVERT(NVERT)')

      CALL HBNAME(1,'FSIPART',oNVCVERT,
     $' NVCVERT[0,300], DIRVERT(3,NVCVERT),
     $    ABSPVERT(NVCVERT), ABSTPVERT(NVCVERT),
     $    IPVERT(NVCVERT),IVERTI(NVCVERT),IVERTF(NVCVERT),
     $    FSIPROB')

      CALL HBNAME(1,'NETARG',oNUMBNDN,
     $' NUMBNDN, NUMBNDP, NUMFREP, NUMATOM')

C--Read card files

      CALL NECARD
      CALL NECARDEV

C--Main loop

      RMAX = C * 2.5

      DO 100 I=1,NEVTEVCT

C      if ( MOD(I,1000).EQ.0 ) WRITE(*,*) 'Processing event: ',I

C--Set vertex
 110     continue
            posv(2)=real(-RMAX+RLU(IDUM)*RMAX*2.)
            posv(3)=real(-RMAX+RLU(IDUM)*RMAX*2.)
 2       RAD   =SQRT(POSV(2)**2+POSV(3)**2)
         IF(RAD.GE.RMAX)GO TO 110
         POSV(1)=real(-SQRT(RMAX**2-RAD**2)+0.0001)

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
         NVC      = 5
         IPVC(5)  = IDPION
         call mcmass(IPVC(5),AMASVC(5))
         PVC(1,5) = real(PTOT)
         PVC(2,5) = 0.
         PVC(3,5) = 0.

         IBOUND   = 1

         POSNUC(1,5)=posv(1)
         POSNUC(2,5)=posv(2)
         POSNUC(3,5)=posv(3)

         IORGVC(5)= -1
         IFLGVC(5)= 0
         ICRNVC(5)= 1
         TIMVC(5) = 0

C-- Generate event
         CALL EFCLFSI
         CALL NEFFECT

C --  nucleon ejection following pion/delta absorption
         call absneject


C --  Fill VCWORK
         NPARVEC=NVC

         DO 40 N1=5,NVC
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
         oNVERT=NVERT

         DO 50 IVERT=1,NVERT
            oIFLGVERT(IVERT) = IFLGVERT(IVERT)
            oPOSVERT(1,IVERT) = POSVERT(1,IVERT)
            oPOSVERT(2,IVERT) = POSVERT(2,IVERT)
            oPOSVERT(3,IVERT) = POSVERT(3,IVERT)

 50      CONTINUE

         oNVCVERT=NVCVERT

         DO 60 IVCVERT=1,NVCVERT
            oIPVERT(IVCVERT) = IPVERT(IVCVERT)
            oDIRVERT(1,IVCVERT) = DIRVERT(1,IVCVERT)
            oDIRVERT(2,IVCVERT) = DIRVERT(2,IVCVERT)
            oDIRVERT(3,IVCVERT) = DIRVERT(3,IVCVERT)
            oABSPVERT(IVCVERT) = ABSPVERT(IVCVERT)
            oABSTPVERT(IVCVERT) = ABSTPVERT(IVCVERT)
            oIVERTI(IVCVERT) = IVERTI(IVCVERT)
            oIVERTF(IVCVERT) = IVERTF(IVCVERT)
 60      CONTINUE

         oFSIPROB = real(EVPIPROB())

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
