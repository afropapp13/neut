************************************************************************
*     ------------------------------
      FUNCTION FNELSPF(ENEU,PF,IPAR)
*     ------------------------------
*
*     (Purpose)
*       Calculate elastic cross section
*
*     (Input)
*       ENEU   : NEUTRINO ENERGY ( GEV )
*       PF(3)  : NECLEON FERMI MOMENTUM ( GEV/C )      
*       IPAR   : NEUTRINO TYPE
*                  12 : NEU-E
*                 -12 : NEU-E-BAR            
*                  14 : NEU-MU
*                 -14 : NEU-MU-BAR            
*                  16 : NEU-TAU
*                 -16 : NEU-TAU-BAR
*     (Output)
*       FNELSPF : ELASTIC CROSS SECTION ( 10**-38 CM**2 )
*
*     (Creation Date and Author)
*       1984.01.08 ; M.NAKAHATA
*       1987.06.26 ; N.SATO     FOR TAU
*       1995.02.06 ; K. KANEYUKI  FOR S.K.
*                    ELSPF -> FNELSPF
*       1995.05.11 ; K.KANEYUKI  GPINDP -> DGAUSS ( INTEGRATION )
*
*       1998.02.23 ; K.Ishihara set PFERMI 217 MeV
*
*       2001.04.24 ; Y.Hayato   Add New COMMON NEINTG3
*                               Remove OLD COMMON NEINTEG2      
*                               USE DNELSEL instead of DNELS
*     
************************************************************************
      IMPLICIT NONE

#include "neutparams.h"

      REAL*4    FNELSPF

      REAL*4    ENEU,PF(3)
      INTEGER*4 IPAR

      REAL*8 Q2MINX, Q2MAXX, EPSIN

      REAL*8 DNELSNL,DGAUSS
C      REAL*8 DNELS,DGAUSS
C      REAL*4 FNQ2MAX,FNQ2MIN,FNNUCL
      REAL*4 FNQ2MAX,FNQ2MIN,FNNPOT
      EXTERNAL DGAUSS,  FNQ2MAX, FNQ2MIN,FNNUCL, DNELSNL

C      COMMON /NEINTEG2/ ENSTOP,AM,AMLEP,JMAX,DPHI,DNEU(3),PNEU(3),
C     &     EV(3),GM,PFERMI,IPR


      INTEGER*4 IPR
      REAL*4 PINNUC(2),APINNUC,EINNUC,EINNEU,PTOTLI(2),PTOTLIS,
     $       ENSTOP,FNNUCLPF,AMIN,AMOUT,AMLEP,PFERMI,
     $       ESUM,AMLEP2,AMOUT2

      COMMON /NEINTG3/ PINNUC,APINNUC,EINNUC,EINNEU,PTOTLI,PTOTLIS,
     $                 ENSTOP,IPR,FNNUCLPF,AMIN,AMOUT,AMLEP,PFERMI,
     $                 ESUM,AMLEP2,AMOUT2

      REAL*4 PNEU(3),BETA,GM,EV(3)
C      REAL*4 base,sina,cosa
      integer*4 iplep
C
C *** Initialize ***      
C
C 98/2/223 K.Ishihara
C      PFERMI=0.250E+0
C      PFERMI=0.217+0
C      JMAX=10
      FNELSPF=0.E+0
      EINNEU=ENEU
      IPR=IPAR
      IF (IPAR.gt.0) then
         CALL MCMASS(2112,AMIN)
         CALL MCMASS(2212,AMOUT)
      ELSE
         CALL MCMASS(2212,AMIN)
         CALL MCMASS(2112,AMOUT)
      ENDIF         
      AMIN =AMIN *1.E-3
      AMOUT=AMOUT*1.E-3

      AMOUT2=AMOUT**2
      
C      write(*,'(A30,F12.5,A1,F12.5,A1,F12.5,A1)') 
C     $     'FNELSPF-Check: Init. PF :(',
C     $     PF(1),',',PF(2),',',PF(3),')'

C      BASE=PF(2)**2+PF(3)**2
C      IF (BASE.gt.0.0001) then
C         SINA=sqrt(PF(3)**2/BASE)
C         COSA=sqrt(PF(2)**2/BASE)
C         PINNUC(2)=COSA*PF(2)-SINA*PF(3)
C      else
C         PINNUC(2)=0.
C      endif
C- should be 0      TMPZ=SINA*PF(2)+COSA*PF(3)
C      write(*,'(A30,F12.5)') 
C     $     'FNELSPF-Check: Z= 0? :',SINA*PF(2)+COSA*PF(3)


      PINNUC(1)=PF(1)
      PINNUC(2)=PF(2)
      IF (abs(PF(3)).gt.0.00001) THEN
         write(*,*) 'FNELSPF: PF(3) SHOULD be 0!'
         stop
      endif

      APINNUC=SQRT(PF(1)**2+PF(2)**2+PF(3)**2)
      EINNUC=SQRT(APINNUC**2+AMIN**2)
C      FNNUCLPF=FNNUCL(APINNUC)
      FNNUCLPF=FNNPOT(1,APINNUC)

      PTOTLI(1)=PINNUC(1)+EINNEU
      PTOTLI(2)=PINNUC(2)

      PTOTLIS=(PTOTLI(1)**2+PTOTLI(2)**2)

C     
C     -- CALCUALTE NEUTRINO ENERGY IN NUCELON STOP SYSTEM
C     
      PNEU(1)=ENEU
      PNEU(2)=0.
      PNEU(3)=0.

      BETA=APINNUC/EINNUC
      GM=1./SQRT(1.-BETA**2)
      IF( APINNUC.NE.0. ) THEN
         EV(1)=-PF(1)/APINNUC
         EV(2)=-PF(2)/APINNUC
         EV(3)=-PF(3)/APINNUC
      ELSE
         EV(1)=1.
         EV(2)=0.
         EV(3)=0.
      ENDIF
      
      CALL MCVECBST(PNEU,0.,EV,GM)
      ENSTOP=SQRT(PNEU(1)**2+PNEU(2)**2+PNEU(3)**2)

      IF(ENSTOP.LT.0.001 .AND. IABS(IPAR).EQ.12)RETURN
      IF(ENSTOP.LT.0.112 .AND. IABS(IPAR).EQ.14)RETURN
      IF(ENSTOP.LT.3.48109.AND. IABS(IPAR).EQ.16)RETURN

C      DNEU(1) = PNEU(1)/ENSTOP
C      DNEU(2) = PNEU(2)/ENSTOP
C      DNEU(3) = PNEU(3)/ENSTOP

C***************************************************************
C     INTEGRATE BY Q**2 (AT THIS TIME CONSIDER PAULI PRINCIPLE)
C***************************************************************
      EV(1) = -EV(1)
      EV(2) = -EV(2)
      EV(3) = -EV(3)

      Q2MAXX = FNQ2MAX(ENSTOP,IPAR)
      Q2MINX = FNQ2MIN(ENSTOP,IPAR)
C      DPHI = 6.283185/FLOAT(JMAX)

      IPLEP=ABS(IPAR)-1
      CALL MCMASS(IPLEP,AMLEP)
      AMLEP=AMLEP*1.0E-3
      AMLEP2=AMLEP**2

      ESUM = EINNEU+EINNUC+FNNUCLPF

C
C****************************************
C     START DGAUSS (CERN lib)
C****************************************
C
C     DOUBLE PRECISION FUNCTION DGAUSS(FUNC,A,B,EPSIN)
C
C     PARAMETERS
C
C     FUNC    = FUNCTION ROUTINE FOR THE FUNCTION FUNC(X).TO BE DE-
C               CLARED EXTERNAL IN THE CALLING ROUTINE
C     A       = LOWER BOUNDARY
C     B       = UPPER BOUNDARY
C     EPSIN   = ACCURACY REQUIRED FOR THE APPROXINATION
C
      EPSIN = 1.0D-7
C      FNELSPF = DGAUSS(DNELS,Q2MINX,Q2MAXX,EPSIN)
      FNELSPF = DGAUSS(DNELSNL,Q2MINX,Q2MAXX,EPSIN)

      RETURN
      END

