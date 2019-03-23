************************************************************************
*     -------------------------
      REAL*8 FUNCTION DNELS(Q2)
*     -------------------------
*
*     (Purpose)
*       Calculate differential elastic cross section d(sigma)/d(q2)
*       with pauli principle      
*
*     (Input)
*       Q2     : Q**2
*
*     (Output)
*       DNELS  : differential elastic cross section
*
*     (Creation Date and Author)
*       1987.06.26 ; N.SATO  PRD.
*       1988.10.07 ; T.KAJITA
*                    EFFECT OF NUCLEON BINDING ENERGY WAS ADDITIONALLY
*                    CONSIDERED
*       1995.02.06 ; K. KANEYUKI  FOR S.K.
*                    FUNELS -> DNELS
*     
************************************************************************

      IMPLICIT NONE
      
      REAL*8 Q2, EINP

#include "neutparams.h"
      
      REAL*8   DNELSq2
      EXTERNAL DNELSQ2
      real*4 fnnpot
      external fnnpot

      real*4  enstop,am,amlep,dphi,dneu(3),pneu(3),ev(3),gm,pfermi
      integer*4 jmax,ipr

      COMMON /NEINTEG2/ ENSTOP,AM,AMLEP,JMAX,DPHI,DNEU,PNEU,EV,
     &     GM,PFERMI,IPR

      REAL*4 DLEP(3),PNP(3)

      integer kpmax
      parameter(kpmax=5)

      real*4 pfabs,el,pl,cost,phi,pabs,eabs
      integer*4 ic,j,kp,jj

      integer*4 getpid,kill
      external  getpid,kill
      integer*4 ipid,iret

      write(*,*) "OBSOLETE FUNCTION. Please do not use."
      ipid = getpid()
C--- use SIGNAL 3 : SIGQUIT to make core-dump file
      iret = kill(ipid,3)
      stop

C---- 
      PFABS=SQRT( AM**2*(GM**2-1.) )
C     
      DNELS=0.
      EL=ENSTOP-Q2/2./AM
      PL=SQRT(EL**2-AMLEP**2)
      COST=(2.*ENSTOP*EL-AMLEP**2-Q2)/2./ENSTOP/PL
      IF(ABS(COST).GT.1.) RETURN

      IC=0
      DO 200 J=1,JMAX
         PHI=(FLOAT(J)-0.5)*DPHI

         DO 220 KP=1,KPMAX
            IF(ABS(COST).LT.1.) CALL NECHADIR(COST,PHI,DNEU,DLEP)
            DO 210 JJ=1,3
               PNP(JJ)=PNEU(JJ)-PL*DLEP(JJ)
 210        CONTINUE
            CALL MCVECBST(PNP,AM,EV,GM)
            PABS=SQRT(PNP(1)**2+PNP(2)**2+PNP(3)**2)
C            EABS=SQRT(PABS**2+AM**2) + (FNNUCL(PFABS)-FNNUCL(PABS))
C            EABS=SQRT(PABS**2+AM**2) + (-0.015)
            EABS=SQRT(PABS**2+AM**2)
            IF(EABS.LT.(AM+0.000001)) THEN
               PABS=0.001
            ELSE
               PABS=SQRT(EABS**2-AM**2)
            END IF

C-comment out 2002/09/17            IF(KP.GE.KPMAX) GO TO 220
C            EL=ENSTOP-Q2/2./AM + (FNNUCL(PFABS)-FNNUCL(PABS))
C            EL=ENSTOP-Q2/2./AM + (-0.015)
C            EL=ENSTOP-Q2/2./AM + FIXPOT
            EL=ENSTOP-Q2/2./AM + (FNNPOT(1,PFABS)-FNNPOT(2,PABS))
            PL=0.001
            IF(EL.LT.(AMLEP+0.000001).AND.KP.GE.KPMAX ) GO TO 200
            IF(EL.LT.(AMLEP+0.000001)) GO TO 220
            PL=SQRT(EL**2-AMLEP**2)
            COST=(2.*ENSTOP*EL-AMLEP**2-Q2)/2./ENSTOP/PL
            IF(ABS(COST).GT.1.) THEN
               IF(KP.GE.KPMAX) GO TO 200
               IF(KP.LT.KPMAX) GO TO 220
            END IF
 220     CONTINUE

         IF(PABS.LT.PFSURF)GO TO 200
         IC=IC+1
 200  CONTINUE

      EINP = ENSTOP
      DNELS=DNELSQ2(EINP,IPR,Q2)*FLOAT(IC)/FLOAT(JMAX)
C      write(6,*) dnels,einp,ipr,q2,ic,jmax
      
      RETURN
      END
