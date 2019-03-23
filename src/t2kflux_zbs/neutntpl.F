C     
C     TEST PROGRAM TO READ VECTOR
C
      IMPLICIT NONE


#include "nework.h"
#include "vcwork.h"
#include "vcvrtx.h"
#include "fsihist.h"
#include "posinnuc.h"

C#include "skhead.h"
C#include "skparm.h"
C#include "sktq.h"

      INTEGER*4 LUNI
      PARAMETER(LUNI=10)

      REAL POS(3)

      INTEGER*4  LRECL
      PARAMETER (LRECL=1024)

      INTEGER*4  ISIZE
      PARAMETER (ISIZE=5000000)

      REAL*4     H(ISIZE)
      COMMON/PAWC/H

      INTEGER*4  IERR,I,N1,N2,ICYCLE,LBANK

C     Neut PARAMETER COMMON
      INTEGER*4   IEVENT,MODEN,NPARNEUT,IPNEUT(100)
      REAL*4      ABSPNEUT(100),PNEUT(3,100)
      COMMON/COM1/IEVENT,MODEN,NPARNEUT,
     $            IPNEUT,ABSPNEUT,PNEUT
C
C     FSIHIST PARAMETER COMMON
      INTEGER NVERTR, IFLGVERTR, NVCVERTR, IPVERTR, 
     $        IVERTIR, IVERTFR, IBOUNDSK
      REAL    POSVERTR, DIRVERTR, ABSPVERTR, ABSTPVERTR,FSIPROBSK
      COMMON /COM3/ NVERTR,POSVERTR(3,100),
     $    IFLGVERTR(100), NVCVERTR, DIRVERTR(3,300),
     $    ABSPVERTR(300), ABSTPVERTR(300),
     $    IPVERTR(300),IVERTIR(300),IVERTFR(300),
     $     FSIPROBSK,IBOUNDSK

C     Vector Parameter COMMON
C
      INTEGER*4   NPARVEC,IORGVEC(100),IPVEC(100),ICHVEC(100),
     $            IFLVEC(100)
      REAL*4      ABSPVEC(100),PVEC(3,100),POSV(3)
      COMMON/COM2/NPARVEC,IPVEC,IORGVEC,ICHVEC,IFLVEC,
     $            ABSPVEC,PVEC,POSV
C  
C     Flux info. common
C     
#ifdef FLUX_10A
#include "uhdef.fh"
#else
#ifdef FLUX_10C
#include "uhdef_10c.fh"
#else
#ifdef FLUX_11A
#include "uhdef_11a.fh"
#else
#include "beamntpl.h"
#endif      
#endif      
#endif

C
C     Target type common
C       COMMON /NEUTTARGET/
#include "necard.h"

C 
C     Cross-section common
C      
#include "neutcrs.h"
      
      integer*4 rlu
      External  rlu
      H(1)= rlu(i)

C
C --  INITIALIZATION
C
      CALL KZINIT
      CALL HLIMIT(-1*ISIZE)
      CALL SKOPENF(LUNI,0,'Z',IERR)
C      SK_GEOMETRY=2
C

      CALL HROPEN(85,'CWNT','genvec.nt','N',LRECL,IERR)

      IF (IERR.ne.0) THEN
         WRITE(6,*) 'HROPEN:FAILED TO OPEN FILE genvec.nt'
         STOP
      ENDIF

      CALL HBNT(10,'NEUTVEC',' ')
      CALL HBNAME(10,'NEUTINFO',IEVENT,
     $'NEV,MODE,NUMNU[0,100],IPNU(NUMNU),ABSPNU(NUMNU),PNU(3,NUMNU)')

      CALL HBNAME(10,'FSIHIST',NVERTR,
     $'NVERT[0,100],POSVERT(3,NVERT),IFLGVERT(NVERT),NVCVERT[0,300],
     $DIRVERT(3,NVCVERT),ABSPVERT(NVCVERT),ABSTPVERT(NVCVERT),
     $IPVERT(NVCVERT),IVERTI(NVCVERT),IVERTF(NVCVERT),FSIPROB,IBOUND')

      CALL HBNAME(10,'VECINFO',NPARVEC,
     $'NPAR[0,100],IPV(NPAR),IORGV(NPAR),ICRNV(NPAR),IFLGV(NPAR),
     $ABSPV(NPAR),PMOMV(3,NPAR),POS(3)')

      CALL HBNAME(10,'CRSINFO',CRSENERGY,
     $'CRSENERGY,TOTCRSNE,DIFCRSNE(8),
     $CRSX,CRSY,CRSZ,CRSPHI,CRSQ2')


#ifdef FLUX_10A
      call hbname(10,'nuinfo',enufd,nuFDtag2)
#else
#if defined(FLUX_10C)||defined(FLUX_11A)
      call hbname(10,'nuinfo',enufd,nuFDtag)
#else
      call hbname(10,'nuinfo',enufd,'Enu,ppid:I,modef:I,ppi,xpi(3),
     $npi(3):R,cospibm,norm:R,ppi0,xpi0(3),npi0(3):R,cospi0bm,rnu,
     $xnu,ynu,nnu(3):R,idFD:I')
#endif
#endif

C      call hbname(10,'photons',nhits,'NHITS:I,POTOT:R')

      CALL HBNAME(10,'TARGINFO',NUMBNDN,
     $'NUMBNDN,NUMBNDP,NUMFREP,NUMATOM')

      CALL HBNAME(10,'NEUTCRS',CRSENERGY,
     $'CRSENERGY,TOTCRSNE,DIFCRSNE(8)')

C
C --  READ VECTOR
C

      DO 10 I=1,50000
C         CALL SKOPTN("31,30,29")
C         CALL SKREAD(LUNI,1000,2000,10,10)
         CALL KZREAD(LUNI,IERR)
         IF (IERR.eq.2) GOTO 1000
         IF (IERR.eq.1) GOTO 2000

         CALL KZBLOC('NEUT',LBANK)
         IF (LBANK.eq.0) GOTO 2000
         IF (IERR.EQ.2) GOTO 1000
         IF (IERR.EQ.1) GOTO 2000
         IF (MOD(I,1000).eq.0) then
             write(6,*)  ' *** ',i,' ***'
         endif
         IEVENT=I

         CALL NERDNEBK(POS)

         CALL VCCLCM
         CALL MCRDHD

         DO 60 N1=1,100
            IPNEUT(N1)=0
            ABSPNEUT(N1)=0.
            IPVEC(N1)=0
            IORGVEC(N1)=0
            ICHVEC(N1)=0
            ABSPVEC(N1)=0.

            DO 70 N2=1,3
               PNEUT(N2,N1)=0.
               PVEC(N2,N1) =0.
 70         CONTINUE
 60      CONTINUE

         MODEN=MODENE         
         NPARNEUT = NUMNE
C         write(*,*) "MODEN=",MODEN," / NPAR=",NPAR

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


         DO 80 N1=1,50
            IFLGVERTR(N1)=0
            DO 90 N2=1,3
               POSVERTR(N2,N1) = 0.
 90         CONTINUE
 80      CONTINUE

         DO 110 N1=1,300
            ABSPVERTR(N1)=0
            ABSTPVERTR(N1)=0
            IPVERTR(N1)=0
            IVERTIR(N1)=0
            IVERTFR(N1)=0
            DO 120 N2=1,3
               DIRVERTR(N2,N1) = 0.
 120        CONTINUE
 110     CONTINUE
         FSIPROB = 0.

         CALL NERDFSIBK
         NVERTR=NVERT
         DO 130 N1=1,NVERT
            IFLGVERTR(N1)=IFLGVERT(N1)
            DO 140 N2=1,3
               POSVERTR(N2,N1) = POSVERT(N2,N1)
 140        CONTINUE
 130     CONTINUE

         NVCVERTR=NVCVERT
         DO 150 N1=1,NVCVERT
            ABSPVERTR(N1) = ABSPVERT(N1) 
            ABSTPVERTR(N1)= ABSTPVERT(N1)
            IPVERTR(N1)   = IPVERT(N1)   
            IVERTIR(N1)   = IVERTI(N1)   
            IVERTFR(N1)   = IVERTF(N1)   
            DO 160 N2=1,3
               DIRVERTR(N2,N1) =  DIRVERT(N2,N1)
 160        CONTINUE
 150     CONTINUE



         CALL VCRDVCCM
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

         CALL KZBLOC('JNDFXVEC',LBANK)         
         IF (LBANK.EQ.0) GOTO 100

         call nerdjndfx
         call nerdnetarg
         call nerdcrsbk
         IBOUNDSK = IBOUND
         FSIPROBSK = FSIPROB


 100     continue

         CALL HFNT(10)
         CALL KZECLR
 10   CONTINUE
      GOTO 3000

 1000 WRITE(6,*) ' READ ERROR : ',LUNI
      GOTO 3000

 2000 WRITE(6,*) ' END OF FILE : ',LUNI

 3000 CONTINUE

      CALL HROUT(10,ICYCLE,' ')
      CALL HREND('CWNT')

C      CALL SKCLOSEF(LUNI)

      STOP
      END
