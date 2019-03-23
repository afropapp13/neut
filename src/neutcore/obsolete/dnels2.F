************************************************************************
*     -------------------------
      REAL*8 FUNCTION DNELS2(Q2)
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
      
      integer*4 iflag
      
#include "neutparams.h"

      parameter(iflag=0)
      REAL*8 Q2, EINP, DNELSQ2
      EXTERNAL DNELSQ2

      real*4  enstop,am,amlep,dphi,dneu(3),pneu(3),ev(3),gm,pfermi
      integer*4 jmax,ipr

      COMMON /NEINTEG2/ ENSTOP,AM,AMLEP,JMAX,DPHI,DNEU,PNEU,EV,
     &     GM,PFERMI,IPAR

      REAL*4  DLEP(3),PNP(3)
      REAL    EVDUM(3)
C---- 

      integer*4 getpid,kill
      external  getpid,kill
      integer*4 ipid,iret

      write(*,*) "OBSOLETE FUNCTION. Please do not use."
      ipid = getpid()
C--- use SIGNAL 3 : SIGQUIT to make core-dump file
      iret = kill(ipid,3)
      stop

      PFABS=SQRT( AM**2*(GM**2-1.) )
C     

      if (iflag.eq.1) then
         DNELS2=0.
         EL=ENSTOP-Q2/2./AM
         PL=SQRT(EL**2-AMLEP**2)
         COST=(2.*ENSTOP*EL-AMLEP**2-Q2)/2./ENSTOP/PL
         IF(ABS(COST).GT.1.) RETURN
         DO 10 KP=1,3
            IF(ABS(COST).LT.1.) CALL NECHADIR(COST,PHI,DNEU,DLEP)
            DO 20 JJ=1,3
               PNP(JJ)=PNEU(JJ)-PL*DLEP(JJ)
 20         CONTINUE
            PNPABS=SQRT(PNP(1)**2+PNP(2)**2+PNP(3)**2)
            IF (PNPABS.LT.1.0E-3) GOTO 40
            DO 30 JJ=1,3
               EVDUM(JJ)=-PNP(JJ)/PNPABS
 30         CONTINUE
            CALL MCVECBST(PNP,AM,EVDUM,GM)
            PABS=SQRT(PNP(1)**2+PNP(2)**2+PNP(3)**2)
C            EABS=SQRT(PABS**2+AM**2) + (FNNUCL(PFABS)-FNNUCL(PABS))
            EABS=SQRT(PABS**2+AM**2) + (FNNPOT(1,PFABS)-FNNPOT(2,PABS))
            IF(EABS.LT.AM) THEN
               PABS=0.001
            ELSE
               PABS=SQRT(EABS**2-AM**2)
            END IF

            IF(KP.GE.3) GO TO 10
C            EL=ENSTOP-Q2/2./AM + (FNNUCL(PFABS)-FNNUCL(PABS))
            EL=ENSTOP-Q2/2./AM + (FNNPOT(1,PFABS)-FNNPOT(2,PABS))
            PL=0.001
            IF(EL.LT.AMLEP.AND.KP.GE.3 ) GO TO 40
            IF(EL.LT.AMLEP ) GO TO 10
            PL=SQRT(EL**2-AMLEP**2)
            COST=(2.*ENSTOP*EL-AMLEP**2-Q2)/2./ENSTOP/PL
            IF(ABS(COST).GT.1.) THEN
               IF(KP.GE.3) GO TO 40
               IF(KP.LT.3) GO TO 10
            END IF
 10      CONTINUE
         IF (PABS.GT.PFERMI) THEN
            EINP = ENSTOP
            DNELS2=DNELSQ2(EINP,IPAR,Q2)
C     write(6,*) ' Q**2 is enough large ',q2,einp,dnels2
            RETURN
         ENDIF
 40      CONTINUE
      endif
      
      DNELS2=0.
      EL=ENSTOP-Q2/2./AM
      PL=SQRT(EL**2-AMLEP**2)
      COST=(2.*ENSTOP*EL-AMLEP**2-Q2)/2./ENSTOP/PL
      IF(ABS(COST).GT.1.) RETURN

      IC=0
      DO 200 J=1,JMAX
         PHI=(FLOAT(J)-0.5)*DPHI

         DO 220 KP=1,3
            IF(ABS(COST).LT.1.) CALL NECHADIR(COST,PHI,DNEU,DLEP)
            DO 210 JJ=1,3
               PNP(JJ)=PNEU(JJ)-PL*DLEP(JJ)
 210        CONTINUE
            CALL MCVECBST(PNP,AM,EV,GM)
            PABS=SQRT(PNP(1)**2+PNP(2)**2+PNP(3)**2)
C            EABS=SQRT(PABS**2+AM**2) + (FNNUCL(PFABS)-FNNUCL(PABS))
            EABS=SQRT(PABS**2+AM**2) + (FNNPOT(1,PFABS)-FNNPOT(2,PABS))
            IF(EABS.LT.AM) THEN
               PABS=0.001
            ELSE
               PABS=SQRT(EABS**2-AM**2)
            END IF

            IF(KP.GE.3) GO TO 220
C            EL=ENSTOP-Q2/2./AM + (FNNUCL(PFABS)-FNNUCL(PABS))
            EL=ENSTOP-Q2/2./AM + (FNNPOT(1,PFABS)-FNNPOT(2,PABS))
            PL=0.001
            IF(EL.LT.AMLEP.AND.KP.GE.3 ) GO TO 200
            IF(EL.LT.AMLEP ) GO TO 220
            PL=SQRT(EL**2-AMLEP**2)
            COST=(2.*ENSTOP*EL-AMLEP**2-Q2)/2./ENSTOP/PL
            IF(ABS(COST).GT.1.) THEN
               IF(KP.GE.3) GO TO 200
               IF(KP.LT.3) GO TO 220
            END IF
 220     CONTINUE

         IF(PABS.LT.PFERMI)GO TO 200
         IC=IC+1
 200  CONTINUE

      EINP = ENSTOP
      DNELS2=DNELSQ2(EINP,IPAR,Q2)*FLOAT(IC)/FLOAT(JMAX)

      RETURN
      END
