C     
C     TEST PROGRAM TO READ VECTOR
C
      IMPLICIT NONE


#include "skhead.h"
#include "nework.h"
#include "necard.h"
#include "nrcard.h"
#include "nefillver.h"
#include "neutmodel.h"
#include "neutparams.h"
#include "fsihist.h"
#include "neutcrs.h"
#include "vcwork.h"
#include "vcvrtx.h"
#include "posinnuc.h"
#include "rscons.h"
#include "posinnuc.h"


      INTEGER*4 LUNI
      PARAMETER(LUNI=10)

      REAL POS(3)

      INTEGER*4  IERR,I,N1,N2,ICYCLE,LBANK

      
      REAL*4 H

      integer*4 rlu
      External  rlu
      H = rlu(i)

C
C --  INITIALIZATION
C
      CALL KZINIT
      CALL SKOPENF(LUNI,0,'Z',IERR)

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
C         IF (MOD(I,1000).eq.0) then
             write(6,*)  ' *** ',i,' ***'
C         endif

         CALL NERDNEBK(POS)

         CALL VCCLCM
         CALL MCRDHD
         CALL NERDFSIBK
         CALL VCRDVCCM
         CALL NERDNETARG

         write(*,*) corev,nucev,nuccv
         write(*,*) SK_GEOMETRY
         write(*,*) "1--",NEFRMFLG,NEPAUFLG,NENEFO16,NEMODFLG,
     $        PFSURF,PFMAX,VNUINI,VNUFIN,IFORMLEN,NUCRESCAT,XNUCFACT
         write(*,*) "2--",MDLQE,MDLSPI,MDLCOH,MDLDIS,MDLQEAF,XMAQE,XMASPI,
     $        MODELDIS,MODELCOH
         write(*,*) "3--",NUMBNDN,NUMBNDP,NUMFREP,NUMATOM,IBOUND 
         write(*,*) "4--",CRSENERGY,TOTCRSNE
         write(*,*) "5--",MDLQE,XMAQE,MDLQEAF,XMAQE,XMVQE,KAPP,PFSF,
     $        SCCFV,SCCFA,FPQE,IRADCORR
         write(*,*) "6--",MDLSPI,NEIFF,NENRTYPE,XMASPI,XMVSPI,RCA5ISPI,
     $        RBGSCLSPI,XMARSRES,XMVRSRES,IPILESSDCY,RPILESSDCY
         write(*,*) "7--",MDLCOH,XMACOH,RAD0NU,fA1COH,fB1COH
         write(*,*) "8--",MDLDIS ,NEPDF  ,NEBODEK
         write(*,*) "9--",PFSURF,PFMAX,VNUINI,VNUFIN,IFORMLEN,FZMU2,
     $        NEFRMFLG,NEPAUFLG,NENEFO16,NENEFMODL,NENEFMODH,
     $        NENEFKINH,NEABSPIEMIT,FEFQE,FEFQEH,FEFINEL,
     $        FEFABS,FEFCOH,FEFCX,FEFCXH,FEFQEHF,FEFCOHF,
     $        FEFCXHF,FEFCOUL,FEFALL

         write(*,*) NUCRESCAT,XNUCFACT,NUCRESFLG

 100     continue

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
