**********************************************************************
*     ------------------
      SUBROUTINE NEFFECT
*     ------------------
*     
*     ( purpose )
*     CONSIDER NUCLEAR EFFECT FOR EACH PION & ETA
*     
*     ( input ) 
*     COMMON VCWORK, VCVRTX
*     
*     ( output ) 
*     COMMON VCWORK, VCVRTX
*     
*     ( creation date and author )
*     1983.12.05 ; M.NAKAHATA
*     1988.06.27 ; T.KAJITA 
*     NUCLEAR EFFECT FOR PAIMOM>800MEV/C IS ADDED
*     1995.03.03 ; K.KANEYUKI
*     NEFFCT -> NEFFECT
*     1997.12.01 ; J.KAMEDA 
*     consider ETA
*     1998.03.05 ; J.KAMEDA 
*     consider  K
*     1998.03.12 ; J.KAMEDA
*     bug fix
*     2002.08.06 ; M.Ishitsuka 
*     consider  omega (only decay)
*     2002.08.27 ; M.Ishitsuka
*     Formation Zones effect is added
*     2003.04.04 ; M.Ishitsuka
*     consider  omega (including nuclear effect)
*     2006.02.23 ; Y.Hayato
*     1) add implicit none
*     2) define common /piinth/   in efpiinth.h
*     3) define common /omegaout/ in efomega.h
*     
*     2006.07.17 ; Y.Hayato
*     particle production   : IFLGVC = 7
*     inelastic scattering  : IFLGVC = 8
*     
*     2007.11.05 ; G.Mitsuka
*     APN,AP are set in necard
*
*     2009.10.xx ; Y.Hayato
*     1) Fix interaction code ( if more than 1 particle interacted,
*                               flags were not correctly set. )
*     2) Avoid storing more than 100 particles. ( array corruption )
*
*     2010.05    ; P. de Perio
*     Modified to work with nucleon ejection
*     Changed EFTRACE arguments NZIG -> NEL,NINEL
*         forward (elastic) scattering : IFLGVC = 9
*
*     2010.09    ; P. de Perio
*     Initialize FSI history vars (NVERT,NVCVERT)
*     
*     2017.11.24 ; Y. Hayato
*     Not to simulate re-interaction of K with momentum larger than 10 GeV.
*      
**********************************************************************
      IMPLICIT NONE

#include "nework.h"
#include "necard.h"
#include "neutparams.h"
#include "vcwork.h"
#include "vcvrtx.h"
#include "efpion.h"

#include "efpiinth.h"
#include "efomega.h"

#include "posinnuc.h"

      REAL*4 X(3),PIN(3),POUT(3),PO1(3),PO2(3)
      REAL*4 PING(3),POUTG(3)
      INTEGER KIN,KFI,NFLAG

      INTEGER NPARTL_PRE

      real*4    RLU
      external  RLU

      REAL*4 APN
C     DATA APN/16./
      REAL*4 AP
C     DATA AP/8./

C     FOR FORMATION ZONES
C     INTEGER IDUM
C     REAL FORMLENGTH, RNFORMLENGTH, X_FORM(3)
C     REAL MU2, FM2GEV
C     DATA MU2/0.08/, FM2GEV/5.07/

      INTEGER*4 I,J
      integer*4 NVC2,NABS,NABS2,NZIG,NSCAT,NDIS,kin1,kin2,ko,NEL,NINEL
      real*4    pabs,dl

C     Clear FSIHIST common block
C     outside of this function since SKDETSIM uses it for possible 
C     multiple secondary interactions 
C      CALL EFCLFSI

C     Define APN,AP
      APN = NUMATOM
      AP  = NUMBNDP
C     
C     ++ ADDED BY M.N ON 23-OCT-88
C     
      NPARTL=0
      NOMEGATL=0
C     
      IF(NVC.LE.4)RETURN

C     
C     -- NEW   INO16  -------------------------  08-OCT-88 -----------
C     IF(INO16.EQ.0)RETURN
C     -- NEW   IBOUND  -------------------------  12-JUNE-06 -----------
      IF(IBOUND.EQ.0)RETURN
C     
C     -- OLD
C     IF(IP(2).EQ.200)GO TO 100
C     RAN=RN(DUM)
C     -- IF FREE PROTON TARGET,NOT CONSIDER NUCLEAR EFFECT
C     IF(RAN.LT.0.2)RETURN
C     ----------------------------------------------------------------
C     
C     
 100  CONTINUE

#ifdef COMPAT_NEUT45
C     C 2006-06-15
C     C 
C     C Position should be set well before this part.
C     C Only for the test
C     C
C     
C     -- SET INTERACTION POSITION
C     
      CALL EFRANBLL(APN,X)
C     
C     -- SAVE INTEACTION POSITION
C     
C     INO16=1
      IBOUND=1

C     CALL UCOPY(X,XIN,3)

C     C for the backward compatibility test
C     C--- should be set well before this
      DO 200 I=1,NVC
         DO 250 J=1,3
            POSNUC(J,I)=X(J)
 250     continue
 200  continue

#else 
C     
C--   Check the positions were correctly set or not.
C     
      DO 200 I=1,NVC
         DO 250 J=1,3
            IF ((IBOUND.eq.-1).or.(POSNUC(J,I).eq.-9999.)) then
               write(*,*) 
     $              'neffect: Interaction positions were not set.'
               write(*,*) 
     $              'neffect: Interaction mode=',modene
               write(*,*) 
     $              'neffect: particle code=',IPVC(I)
               stop
            endif
 250     continue
 200  continue

#endif

      NPARTL=0

      NVC2=NVC
      DO 1000 I=5,NVC

C     NPARTL_PRE = 0
         NPARTL_PRE = NPARTL

C     -YH 18-Feb-2004
C     Avoid tracking "erased" particle.

 959     format(A3,I7,A6,I2,A3,F10.3,A2,F9.3,A1,F9.3,A1,F9.3,A1)

         IF (ICRNVC(I).ne.1) then

C     IF (IPVC(I).ne.200000) then
C     write(*,*) "Erased interaction found"
C     do 955 IJK=1,NVC
C     write(*,959)
C     $                 'IP',IPVC(IJK),' ICRN=',ICRNVC(IJK),
C     $                 ' P=',
C     $                 sqrt(PVC(1,IJK)**2+PVC(2,IJK)**2+PVC(3,IJK)**2),
C     $                 ':(',PVC(1,IJK),',',
C     $                 PVC(2,IJK),',',PVC(3,IJK),')'
C     955           continue
C     endif

            goto 1000
         endif

CCCC  ------- for PION ----------
         IF (ABS(IPVC(I)).EQ.211 .OR. IPVC(I).EQ.111) THEN
            CALL UCOPY(PVC(1,I),PIN,3)
            PABS=SQRT(PIN(1)**2+PIN(2)**2+PIN(3)**2)
CCCC  IF(PABS.GT.850.)GO TO 1000
CCCC  T.KAJITA  27-JUN-88
            KIN=IPVC(I)
            LOOPPI=I
C     C FORMATION ZONES
C     FORMLENGTH=PABS/1000./MU2/FM2GEV
C     RNFORMLENGTH=-1.*FORMLENGTH*LOG(RLU(IDUM))
C     DO J=1,3
C     X_FORM(J)=X(J)+(PIN(J)/PABS)*RNFORMLENGTH
C     END DO
C     CALL EFTRACE(APN,AP,KIN,X_FORM,PIN,NABS,NZIG,KFI,POUT)

            
            DO J=1,3
               X(J)=POSNUC(J,I)
            ENDDO

            CALL EFTRACE(APN,AP,KIN,X,PIN,NABS,NEL,NINEL,KFI,POUT)
            NZIG = NEL+NINEL
C     
C     --        ABSORPTION
C     
C            IF (KFI.EQ.100000) THEN
            IF (NABS.EQ.0) THEN
C     IF (NPARTL.GT.0) THEN
               IF (NPARTL.GT.NPARTL_PRE) THEN
                  IFLGVC(I)=7
               else
                  IFLGVC(I)=3
               endif
               ICRNVC(I)=0
               CALL UCOPY(POSVC,POSFVC(1,I),3)
               IVTFVC(I)=1
C     
C     --        CHARGE EXCHANGE And scattering
C     
            ELSE IF ((KIN.NE.KFI).or.(NZIG.gt.0)) THEN
               IF (KIN.NE.KFI) THEN
                  IFLGVC(I)=4
               ELSE
C     IF (NPARTL.GT.0) THEN
                  IF (NPARTL.GT.NPARTL_PRE) THEN
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
               IF (NVC2.LT.MAXVC) THEN
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
C--   escaped from nucleus
                     POSNUC(J,NVC2)=9999.
                  enddo
               endif
C     
C     --        OTHER CASE (Momentum can change if Coulomb correction is on)
C     
            ELSE IF (FEFCOUL.eq.0) then
               DO J=1,3
                  if (POUT(J).ne.PVC(J,I)) then
                     write(*,*) 'NEFFECT: Momentum changed without 
     &                           interaction'
                     write(*,*) 'NEFFECT: Mode=',modene
                     stop
                  endif
               enddo
C     CALL UCOPY(POUT,PVC(1,I),3)
            ENDIF
         ENDIF
CCCC  ---------------     for ETA    -------------------
CCCC  --------------------------------------------------

         IF (IPVC(I).EQ.221) THEN

            CALL UCOPY(PVC(1,I),PIN,3)
C     C FORMATION ZONES
C     PABS=SQRT(PIN(1)**2+PIN(2)**2+PIN(3)**2)
C     FORMLENGTH=PABS/1000./MU2/FM2GEV
C     RNFORMLENGTH=-1.*FORMLENGTH*LOG(RLU(IDUM))
C     DO J=1,3
C     X_FORM(J)=X(J)+(PIN(J)/PABS)*RNFORMLENGTH
C     END DO
C     CALL EFTRCETA(X_FORM,PIN,NSCAT,NDIS,PO1,PO2,KIN1,KIN2)

            DO J=1,3
               X(J)=POSNUC(J,I)
            ENDDO
            CALL EFTRCETA(X,PIN,NSCAT,NDIS,NABS,NABS2,PO1,PO2,KIN1,KIN2)

            IF (NDIS.EQ.0) then
               if (NSCAT.GT.0) THEN
                  IF (NVC2.LT.MAXVC) THEN
                     NVC2 = NVC2+1
                     IFLGVC(I)=8
                     ICRNVC(I)=0

                     IPVC(NVC2)=IPVC(I)
                     IORGVC(NVC2)=I
                     IFLGVC(NVC2)=0
                     ICRNVC(NVC2)=1
                     IVTIVC(NVC2)=1
                     CALL UCOPY(PO1,PVC(1,NVC2),3)
                     CALL MCMASS(IPVC(NVC2),AMASVC(NVC2))
                     DO J=1,3
C--   escaped from nucleus
                        POSNUC(J,NVC2)=9999.
                     enddo
                  else
                     CALL UCOPY(PO1,PVC(1,I),3)                     
                  endif
               endif
            endif 
CCCC  ------------    change to pi    -----------------
            IF (NDIS.ge.1) then
               IFLGVC(I)=3
               ICRNVC(I)=0
               CALL UCOPY(POSVC,POSFVC(1,I),3)
               IVTFVC(I)=1
               IF (NVC2.LT.MAXVC) THEN
                  NVC2 = NVC2+1
                  IPVC(NVC2)=KIN1
                  CALL UCOPY(PO1,PVC(1,NVC2),3)
                  IORGVC(NVC2)=I
                  IFLGVC(NVC2)=0
                  ICRNVC(NVC2)=1
                  CALL UCOPY(POSVC,POSIVC(1,NVC2),3)
                  IVTIVC(NVC2)=1
                  CALL MCMASS(IPVC(NVC2),AMASVC(NVC2))
                  DO J=1,3
C--   escaped from nucleus
                     POSNUC(J,NVC2)=9999.
                  enddo
                  
c     c   -------------   absorption -------------      
                  IF (NABS.eq.0) THEN
                     IFLGVC(NVC2)=3
                     ICRNVC(NVC2)=0
                     CALL UCOPY(POSVC,POSFVC(1,NVC2),3)
                  ENDIF
c     c   -----------------------------------------
C     
                  IF (NDIS.ge.2) then
                     IF (NVC2.LT.MAXVC) THEN
                        NVC2 = NVC2+1
                        IPVC(NVC2)=KIN2
                        CALL UCOPY(PO2,PVC(1,NVC2),3)
                        IORGVC(NVC2)=I
                        IFLGVC(NVC2)=0
                        ICRNVC(NVC2)=1
                        CALL UCOPY(POSVC,POSIVC(1,NVC2),3)
                        IVTIVC(NVC2)=1
                        CALL MCMASS(IPVC(NVC2),AMASVC(NVC2))
                        DO J=1,3
                           POSNUC(J,NVC2)=9999.
                        enddo
c     c   -------------   absorption -------------      
                        IF (NABS2.eq.0) THEN
                           IFLGVC(NVC2)=3
                           ICRNVC(NVC2)=0
                           CALL UCOPY(POSVC,POSFVC(1,NVC2),3)
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
cccc  -----------------------------------------------

cccc  --------------  for K ---------------------
cccc  ------      add 1998.03.05      ----------- 
cccc  --  K+:321 K-:-321 K0:311 KL:130 KS:310  --

         IF ((ABS(IPVC(I)).EQ.321).or.(IPVC(I).EQ.130).or.
     &        (IPVC(I).eq.310)) THEN
            CALL UCOPY(PVC(1,I),PIN,3)
            KIN = IPVC(I)
C     C FORMATION ZONES
            PABS=SQRT(PIN(1)**2+PIN(2)**2+PIN(3)**2)
C     FORMLENGTH=PABS/1000./MU2/FM2GEV
C     RNFORMLENGTH=-1.*FORMLENGTH*LOG(RLU(IDUM))
C     DO J=1,3
C     X_FORM(J)=X(J)+(PIN(J)/PABS)*RNFORMLENGTH
C     END DO

C-- skip K with momentum > 10GeV ( 2017/11/24 )
            if (PABS.gt.10000) then
               goto 1000
            endif

            DO J=1,3
               X(J)=POSNUC(J,I)
C-- add to convert from MeV to GeV(2020/03/15)
               PING(J)=PIN(J)/1000. 
            ENDDO

C            CALL EFTRCKA(PIN,KIN,POUT,KO,X) 

C-- call with GeV(2020/03/15)
            CALL EFTRCKA(PING,KIN,POUTG,KO,X) 

C-- add to convert back from GeV to MeV (2020/03/15)
            DO J=1,3
               POUT(J)=POUTG(J)*1000.
            ENDDO

            IF (KIN.EQ.KO) then
               if ((( abs(pin(1)-pout(1)) ).gt.0.01 ).or.
     $             (( abs(pin(2)-pout(2)) ).gt.0.01 ).or.
     $             (( abs(pin(3)-pout(3)) ).gt.0.01 )) THEN
                  IF (NVC2.LT.MAXVC) THEN
                     NVC2 = NVC2+1
                     IFLGVC(I)=8
                     ICRNVC(I)=0
                     CALL UCOPY(POSVC,POSFVC(1,I),3)                     
                     IPVC(NVC2)=IPVC(I)
                     IORGVC(NVC2)=I
                     IFLGVC(NVC2)=0
                     ICRNVC(NVC2)=1
                     IVTIVC(NVC2)=1
                     CALL UCOPY(POUT,PVC(1,NVC2),3)
                     CALL UCOPY(POSVC,POSIVC(1,NVC2),3)                     
                     CALL MCMASS(IPVC(NVC2),AMASVC(NVC2))
                     DO J=1,3
                        POSNUC(J,NVC2)=9999.
                     enddo
                  else
                     CALL UCOPY(POUT,PVC(1,I),3)                     
                  endif
               endif
            endif
ccc   --------     charge exchange ----------------
            IF (KIN.NE.KO) then

               IFLGVC(I)=4
               ICRNVC(I)=0
               CALL UCOPY(POSVC,POSFVC(1,I),3)
               IVTFVC(I)=1
               IF (NVC2.LT.MAXVC) THEN
                  NVC2=NVC2+1
                  IPVC(NVC2)=KO
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
            ENDIF
         ENDIF

cccc  ---------------   for omega  --------------
         IF (IPVC(I).EQ.223) THEN
            CALL UCOPY(PVC(1,I),PIN,3)
            KIN=IPVC(I)
            LOOPOMEGA=I
C     C FORMATION ZONES
C     PABS=SQRT(PIN(1)**2+PIN(2)**2+PIN(3)**2)
C     FORMLENGTH=PABS/1000./MU2/FM2GEV
C     RNFORMLENGTH=-1.*FORMLENGTH*LOG(RLU(IDUM))
C     DO J=1,3
C     X_FORM(J)=X(J)+(PIN(J)/PABS)*RNFORMLENGTH
C     END DO
C     CALL EFTRCOMEGA(X_FORM,PIN,KIN,NFLAG)

            DO J=1,3
               X(J)=POSNUC(J,I)
            ENDDO
            CALL EFTRCOMEGA(X,PIN,KIN,NFLAG)

C     
C     --        DECAY
C     
            IF (NFLAG.EQ.1) THEN
               IFLGVC(I)=3
               ICRNVC(I)=0
               CALL UCOPY(POSVC,POSFVC(1,I),3)
               IVTFVC(I)=1
            ELSE
               WRITE(6,782)
 782           FORMAT('  ********** ERROR IN OMEGA TRACE ******')
            ENDIF
         ENDIF

cccc  -------------------------------------------

 1000 CONTINUE
      NVC=NVC2
      IF(NPARTL.EQ.0) GOTO 1101
C     
C     ++ ADDED BY T.KAJITA 27-JUN-88
C     -- RE-CONSIDER PION INTERACTION FOR PIONS PRODUCED IN PAI-N INT.
C     IN 16O
C     
      NVC2=NVC+NPARTL
      IF (NVC2.GT.MAXVC) NVC2=MAXVC
      DO 1100 I=1,NPARTL
C     NPARTL_PRE = 0
         NPARTL_PRE = NPARTL
         IF ( (NVC+I).GT.MAXVC ) GOTO 1100
         IPVC(NVC+I)=IPINT(I)
         CALL UCOPY(PINT(1,I),PVC(1,NVC+I),3)
         IORGVC(NVC+I)=IPIORI(I)
         IFLGVC(NVC+I)=0
         ICRNVC(NVC+I)=1
         CALL UCOPY(POSVC,POSIVC(1,NVC+I),3)
         IVTIVC(NVC+I)=1
         CALL MCMASS(IPVC(NVC+I),AMASVC(NVC+I))
         CALL UCOPY(PINT(1,I),PIN(1),3)
         CALL UCOPY(XINT(1,I),X,3)
         DO J=1,3
            POSNUC(J,NVC+I)=X(J)
         enddo
         PABS=SQRT(PIN(1)**2+PIN(2)**2+PIN(3)**2)
         IF (PABS.GT.850.) GOTO 1100
         KIN=IPINT(I)
C     C FORMATION ZONES
C     FORMLENGTH=PABS/1000./MU2/FM2GEV
C     RNFORMLENGTH=-1.*FORMLENGTH*LOG(RLU(IDUM))
C     DO J=1,3
C     X_FORM(J)=X(J)+(PIN(J)/PABS)*RNFORMLENGTH
C     END DO
C     CALL EFTRACE(APN,AP,KIN,X_FORM,PIN,NABS,NZIG,KFI,POUT)

         CALL EFTRACE(APN,AP,KIN,X,PIN,NABS,NEL,NINEL,KFI,POUT)
         NZIG = NEL+NINEL
C     
C     --     ABSORPTION
C     
C         IF (KFI.EQ.100000) THEN
         IF (NABS.EQ.0) THEN
C     IF (NPARTL.GT.0) THEN
            IF (NPARTL.GT.NPARTL_PRE) THEN
               IFLGVC(NVC+I)=7
            else
               IFLGVC(NVC+I)=3
            endif
            ICRNVC(NVC+I)=0
            CALL UCOPY(POSVC,POSFVC(1,NVC+I),3)
            IVTFVC(NVC+I)=1
C     
C     --     CHARGE EXCHANGE  And scattering
C     
         ELSE IF ((KIN.NE.KFI).or.(NZIG.gt.0)) THEN
            IF (KIN.NE.KFI) THEN
               IFLGVC(NVC+I)=4
            ELSE
C     IF (NPARTL.GT.0) THEN
               IF (NPARTL.GT.NPARTL_PRE) THEN
                  IFLGVC(NVC+I)=7
               ELSE IF (NINEL.gt.0) THEN
                  IFLGVC(NVC+I)=8
               ELSE 
                  IFLGVC(NVC+I)=9
               ENDIF
            ENDIF
            ICRNVC(NVC+I)=0
            CALL UCOPY(POSVC,POSFVC(1,NVC+I),3)
            IVTFVC(NVC+I)=1
            IF (NVC2.LT.MAXVC) THEN
               NVC2=NVC2+1
               IPVC(NVC2)=KFI
               CALL UCOPY(POUT,PVC(1,NVC2),3)
               IORGVC(NVC2)=NVC+I
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
C     --     OTHER CASE (Momentum can change if Coulomb correction is on)
C     
         ELSE IF (FEFCOUL.eq.0) then
C     
            DO J=1,3
               if (POUT(J).ne.PVC(J,NVC+I)) then
                  write(*,*) 
     $                 'NEFFECT: Momentum changed without interaction'
                  write(*,*) 
     $                 'NEFFECT: Mode=',modene
                  stop
               endif
            enddo
C     CALL UCOPY(POUT,PVC(1,NVC+I),3)
         ENDIF
 1100 CONTINUE
      NVC=NVC2
      IF (NVC.GT.MAXVC) NVC=MAXVC
C     
 1101 CONTINUE
C     
cccc  -------------------------------------------
      IF(NOMEGATL.EQ.0) GOTO 1201
C     
C     ++ ADDED BY M.ISHITSUKA 
C     
      NVC2=NVC+NOMEGATL
      IF (NVC2.GT.MAXVC) NVC2=MAXVC
      DO 1200 I=1,NOMEGATL
         IF ((NVC+I).GT.MAXVC) GOTO 1200
         IPVC(NVC+I)=IPINTOMEGA(I)
         CALL UCOPY(PINTOMEGA(1,I),PVC(1,NVC+I),3)
         IORGVC(NVC+I)=IPIORIOMEGA(I)
         IFLGVC(NVC+I)=0
         ICRNVC(NVC+I)=1
         CALL UCOPY(POSVC,POSIVC(1,NVC+I),3)
         IVTIVC(NVC+I)=1   
         CALL MCMASS(IPVC(NVC+I),AMASVC(NVC+I))
         CALL UCOPY(PINTOMEGA(1,I),PIN(1),3)
         CALL UCOPY(XINTOMEGA(1,I),X,3)
         PABS=SQRT(PIN(1)**2+PIN(2)**2+PIN(3)**2)

         DO J=1,3
            POSNUC(J,NVC+I)=X(J)
         ENDDO

         IF (PABS.GT.850.) GOTO 1200
         KIN=IPINTOMEGA(I)
C     C FORMATION ZONES
C     FORMLENGTH=PABS/1000./MU2/FM2GEV
C     RNFORMLENGTH=-1.*FORMLENGTH*LOG(RLU(IDUM))
C     DO J=1,3
C     X_FORM(J)=X(J)+(PIN(J)/PABS)*RNFORMLENGTH
C     END DO
C     DL=SQRT(X_FORM(1)**2+X_FORM(2)**2+X_FORM(3)**2)
C     IF (DL.GT.CC2) GOTO 1200
C     CALL EFTRACE(APN,AP,KIN,X_FORM,PIN,NABS,NZIG,KFI,POUT)
         
         CALL EFTRACE(APN,AP,KIN,X,PIN,NABS,NEL,NINEL,KFI,POUT)
         NZIG = NEL+NINEL
C     
C     --     ABSORPTION
C     
C         IF (KFI.EQ.100000) THEN
         IF (NABS.EQ.0) THEN
            IFLGVC(NVC+I)=3
            ICRNVC(NVC+I)=0
            CALL UCOPY(POSVC,POSFVC(1,NVC+I),3)
            IVTFVC(NVC+I)=1
C     
C     --     CHARGE EXCHANGE
C     
C     --         And Inelastic scattering
C     
         ELSE IF ((KIN.NE.KFI).or.(NZIG.gt.0)) THEN
            IF (KIN.NE.KFI) THEN
               IFLGVC(NVC+I)=4
            ELSE
C     IF (NPARTL.GT.0) THEN
               IF (NPARTL.GT.NPARTL_PRE) THEN
                  IFLGVC(NVC+I)=7
               ELSE IF (NINEL.gt.0) THEN
                  IFLGVC(NVC+I)=8
               ELSE 
                  IFLGVC(NVC+I)=9
               ENDIF
            ENDIF
            ICRNVC(NVC+I)=0
            CALL UCOPY(POSVC,POSFVC(1,NVC+I),3)
            IVTFVC(NVC+I)=1
            IF (NVC2.LT.MAXVC) THEN
               NVC2=NVC2+1
               IPVC(NVC2)=KFI
               CALL UCOPY(POUT,PVC(1,NVC2),3)
               IORGVC(NVC2)=NVC+I
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
C     --     OTHER CASE (Momentum can change if Coulomb correction is on)
C     
         ELSE IF (FEFCOUL.eq.0) then
C     
            DO J=1,3
               if (POUT(J).ne.PVC(J,NVC+I)) then
                  write(*,*) 
     $                 'NEFFECT: Momentum changed without interaction'
                  write(*,*) 
     $                 'NEFFECT: Mode=',modene
                  stop
               endif
            enddo
C     CALL UCOPY(POUT,PVC(1,NVC+I),3)
         ENDIF
 1200 CONTINUE
      NVC=NVC2
      IF (NVC.GT.MAXVC) NVC=MAXVC
C     
 1201 CONTINUE
C     
     
      RETURN
      END


