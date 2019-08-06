************************************************************************
*     --------------------------------------
      SUBROUTINE NEMECVCP(IPAR,E,DNEUT,IERR)
*     --------------------------------------
*
*     (Purpose)
*       Vector generation for mec event
*
*     (Input)
*       IPAR   : NEUTRINO TYPE
*                  12 : NEU-E
*                 -12 : NEU-E-BAR
*                  14 : NEU-MU
*                 -14 : NEU-MU-BAR
*                  16 : NEU-TAU
*                 -16 : NEU-TAU-BAR
*       ENE    : NEUTRINO ENERGY ( GEV )
*       DNEUT(3) : DIRECTION OF NEUTRINO
*
*     (Output)
*       IERR    : ERROR CODE
*       COMMON /NEWORK/
*
*     (Creation Date and Author)
*       2013.02.26 ; P.Sinclair (based on NEELSVCP)
*
************************************************************************
      IMPLICIT NONE

#include "nework.h"
#include "necard.h"


      INTEGER*4 ipar
      real*4    e,cost,theta,phi,eta,dum
      real*4    abspnnc1,pinnuc1,einnuc1,MASSIN1,eoutp1,MASSOUT1
      real*4    abspnnc2,pinnuc2,einnuc2,MASSIN2,eoutp2,MASSOUT2
      real*4    einit,p4in,eoutl

      real*4    efinal,p4out,amlep

      integer*4 ierr,icall,ierr2,ierr3,iii,i,j

      real*4    rlu,fnnpot
      external  rlu,fnnpot

      INTEGER COUNTER
      REAL*4  DNEUT(3),PLEP(3),PF1(3),PF2(3),PBUF0(3),PBUF1(3),PBUF2(3)
      REAL*4  ELEP
C IPMODE(6) contains PID of particles involved in interaction
C 1 = incoming nucleon, 2 = target nucleon(1), 3 = target nucleon(2)
C 3 = outgoing lepton, 3 = outgoing nucleon(1), 4 = outgoing nucleon(2)
      INTEGER IPMODE(6)
      INTEGER IORG(6)/ 0, 0, 0, 1, 2, 3/
      INTEGER IFLG(6)/-1,-1,-1, 0, 0, 0/
      INTEGER ICRN(6)/ 0, 0, 0, 1, 1, 1/

C Final nucleon momenta
      REAL*4 PN1(3),PN2(3)
C Energy and 3-momentum transfer from lepton to hadronic system
      REAL*4 Q0, Q3(3)

      IPMODE(1)=IPAR
      IF (IPAR.GT.0) THEN
        IPMODE(4)=IPAR-1
      ELSE
        IPMODE(4)=IPAR+1
      ENDIF

C For the time being, return an error if the incoming neutrino is not a numu
      IF ((ABS(IPAR).NE.14).AND.(ABS(IPAR).NE.12)) GOTO 800

C
C --  SET INITIAL NEUTRINO DIRECTION TO X-AXIS
C
      COST=DNEUT(3)
      THETA=ACOS(COST)
      IF (DNEUT(1).EQ.0. .AND. DNEUT(2).EQ.0.) THEN
         PHI=0.
      ELSE
         PHI=ATAN2(DNEUT(2),DNEUT(1))
      ENDIF

C
C --  Set final lepton momentum, then try to find a pair of initial state nucleon momenta which will allow the event
C
       IERR3=-1
       COUNTER=0
       DO WHILE ((IERR3.ne.0).and.(COUNTER.lt.20))
         IERR2=0
         CALL NEMECLPV(E,ELEP,PLEP,IPMODE,IERR2)
         IF (IERR2.NE.0) GOTO 900
         Q0 = E-ELEP
         Q3(1) = E-PLEP(1)
         Q3(2) = 0.-PLEP(2)
         Q3(3) = 0.-PLEP(3)
         IERR3=0
         CALL NEMECHAD(Q0,Q3,PF1,PF2,IPMODE,PN1,PN2,IERR3)
         COUNTER=COUNTER+1
       END DO
       IF(IERR3.ne.0) GOTO 910

C
C     -- BACK TO INITIAL COORDINATE
C     - RANDOM ROTATION AROUND X-AXIS
C
            ETA=RLU(DUM)*6.283185
            DO 20 III=1,5
               IF (III.EQ.1) THEN
                  PBUF0(1)=PLEP(1)
                  PBUF0(2)=PLEP(2)
                  PBUF0(3)=PLEP(3)
               ELSE IF (III.eq.2) THEN
                  PBUF0(1)=PF1(1)
                  PBUF0(2)=PF1(2)
                  PBUF0(3)=PF1(3)
               ELSE IF (III.eq.3) THEN
                  PBUF0(1)=PF2(1)
                  PBUF0(2)=PF2(2)
                  PBUF0(3)=PF2(3)
               ELSE IF (III.eq.4) THEN
                  PBUF0(1)=PN1(1)
                  PBUF0(2)=PN1(2)
                  PBUF0(3)=PN1(3)
               ELSE
                  PBUF0(1)=PN2(1)
                  PBUF0(2)=PN2(2)
                  PBUF0(3)=PN2(3)
               ENDIF
               PBUF1(1)=PBUF0(1)
               PBUF1(2)=COS(ETA)*PBUF0(2)-SIN(ETA)*PBUF0(3)
               PBUF1(3)=SIN(ETA)*PBUF0(2)+COS(ETA)*PBUF0(3)
C
C     - THETA ROTATION
C
               PBUF2(1)=SIN(THETA)*PBUF1(1)-COS(THETA)*PBUF1(3)
               PBUF2(2)=PBUF1(2)
               PBUF2(3)=COS(THETA)*PBUF1(1)+SIN(THETA)*PBUF1(3)
C
C     - PHI ROTATION
C
               PBUF0(1)=COS(PHI)*PBUF2(1)-SIN(PHI)*PBUF2(2)
               PBUF0(2)=SIN(PHI)*PBUF2(1)+COS(PHI)*PBUF2(2)
               PBUF0(3)=PBUF2(3)
               IF (III.EQ.1) THEN
                  PLEP(1)=PBUF0(1)
                  PLEP(2)=PBUF0(2)
                  PLEP(3)=PBUF0(3)
               ELSE IF (III.eq.2) THEN
                  PF1(1)=PBUF0(1)
                  PF1(2)=PBUF0(2)
                  PF1(3)=PBUF0(3)
               ELSE IF (III.eq.3) THEN
                  PF2(1)=PBUF0(1)
                  PF2(2)=PBUF0(2)
                  PF2(3)=PBUF0(3)
               ELSE IF (III.eq.4) THEN
                  PN1(1)=PBUF0(1)
                  PN1(2)=PBUF0(2)
                  PN1(3)=PBUF0(3)
               ELSE
                  PN2(1)=PBUF0(1)
                  PN2(2)=PBUF0(2)
                  PN2(3)=PBUF0(3)
               ENDIF
 20         CONTINUE
C
C     -- NOW STORE VECTOR
C
            IF ( IPAR.GT.0 ) MODENE= 2
            IF ( IPAR.LT.0 ) MODENE=-2
            NUMNE = 6
            DO 30 I=1,6
               IPNE(I)=IPMODE(I)
               IORGNE(I)=IORG(I)
               IFLGNE(I)=IFLG(I)
               ICRNNE(I)=ICRN(I)
               DO 40 J=1,3
                  IF (I.EQ.1) THEN
                    PNE(J,I)=E*DNEUT(J)
                  ELSE IF (I.EQ.2) THEN
                    PNE(J,I)=PF1(J)
                  ELSE IF (I.EQ.3) THEN
                    PNE(J,I)=PF2(J)
                  ELSE IF (I.EQ.4) THEN
                    PNE(J,I)=PLEP(J)
                  ELSE IF (I.EQ.5) THEN
                    PNE(J,I)=PN1(J)
                  ELSE IF (I.EQ.6) THEN
                    PNE(J,I)=PN2(J)
                  ENDIF
 40            CONTINUE
 30         CONTINUE

            CALL MCMASS(IPNE(4),AMLEP)
            AMLEP=AMLEP*1.0E-3
            CALL MCMASS(IPNE(2),MASSIN1)
            MASSIN1=MASSIN1*1.0E-3
            CALL MCMASS(IPNE(3),MASSIN2)
            MASSIN2=MASSIN2*1.0E-3
            CALL MCMASS(IPNE(5),MASSOUT1)
            MASSOUT1=MASSOUT1*1.0E-3
            CALL MCMASS(IPNE(6),MASSOUT2)
            MASSOUT2=MASSOUT2*1.0E-3

            ABSPNNC1 =sqrt(PNE(1,5)**2+PNE(2,5)**2+PNE(3,5)**2)
            ABSPNNC2 =sqrt(PNE(1,6)**2+PNE(2,6)**2+PNE(3,6)**2)

C
C --        NORMAL ENDING
C
            if (QUIET.eq.0) then
            write(*,'(A10,F12.5,A1,F12.5,A1,F12.5,A1,F12.5)')
     $           'INLEP =(',PNE(1,1),',',PNE(2,1),',',PNE(3,1),')'
            write(*,'(A10,F12.5,A1,F12.5,A1,F12.5,A1,F12.5)')
     $           'INNUC1 =(',PNE(1,2),',',PNE(2,2),',',PNE(3,2),')'
            write(*,'(A10,F12.5,A1,F12.5,A1,F12.5,A1,F12.5)')
     $           'INNUC2 =(',PNE(1,3),',',PNE(2,3),',',PNE(3,3),')'
            endif
            PINNUC1=sqrt(PNE(1,2)**2+PNE(2,2)**2+PNE(3,2)**2)
            EINNUC1=sqrt(PNE(1,2)**2+PNE(2,2)**2+PNE(3,2)**2+MASSIN1**2)

            PINNUC2=sqrt(PNE(1,3)**2+PNE(2,3)**2+PNE(3,3)**2)
            EINNUC2=sqrt(PNE(1,3)**2+PNE(2,3)**2+PNE(3,3)**2+MASSIN2**2)

C Einit = sum of energy of all initial particles in lab frame
C P4in = (Sum of 4-vectors of all initial particles)^2
            Einit= EINNUC1+EINNUC2+E
            p4in = Einit**2-
     $            ( (PNE(1,1)+PNE(1,2)+PNE(1,3) )**2
     $             +(PNE(2,1)+PNE(2,2)+PNE(2,3) )**2
     $             +(PNE(3,1)+PNE(3,2)+PNE(3,3) )**2 )

            if (QUIET.eq.0) then
            write(*,'(A10,F12.5,A1,F12.5,A1,F12.5,A1,F12.5)')
     $           'OUTLEP=(',PNE(1,4),',',PNE(2,4),',',PNE(3,4),')'
            endif
            EOUTL =sqrt(PNE(1,4)**2+PNE(2,4)**2+PNE(3,4)**2
     $             +AMLEP**2)

            if (QUIET.eq.0) then
            write(*,'(A10,F12.5,A1,F12.5,A1,F12.5,A1,F12.5)')
     $           'OUTNUC1=(',PNE(1,5),',',PNE(2,5),',',PNE(3,5),')'
            endif
            EOUTP1 =sqrt(PNE(1,5)**2+PNE(2,5)**2+PNE(3,5)**2
     $             +MASSOUT1**2)

            if (QUIET.eq.0) then
            write(*,'(A10,F12.5,A1,F12.5,A1,F12.5,A1,F12.5)')
     $           'OUTNUC2=(',PNE(1,6),',',PNE(2,6),',',PNE(3,6),')'
            endif
            EOUTP2 =sqrt(PNE(1,6)**2+PNE(2,6)**2+PNE(3,6)**2
     $             +MASSOUT2**2)

            EFINAL=EOUTL+EOUTP1+EOUTP2
            p4out = (EFINAL)**2-
     $            ( (PNE(1,4)+PNE(1,5)+PNE(1,6) )**2
     $             +(PNE(2,4)+PNE(2,5)+PNE(2,6) )**2
     $             +(PNE(3,4)+PNE(3,5)+PNE(3,6) )**2 )

            if (QUIET.eq.0) then
            write(*,'(A10,F12.5,A10,F12.5)')
     $           'Ein : ',Einit*1000., 'Eout : ',EFINAL*1000.

            write(*,'(A10,F12.5,A10,F12.5,A10,F12.5)')
     $         'Ediff : ',(Einit-Efinal)*1000.,
     $         'VnuclD1:' ,(FNNPOT(2,ABSPNNC1)-FNNPOT(1,PINNUC1))*1000.,
     $         'VnuclD2:' ,(FNNPOT(2,ABSPNNC2)-FNNPOT(1,PINNUC2))*1000.

            write(*,'(A10,F12.5,A10,F12.5)')
     $           'P4in : ',P4in, 'P4out : ',P4out
            endif

            IERR=0
            RETURN
C
C --   ERROR RETURN
C
 900  WRITE(6,901)
 901  FORMAT(' *** ERROR AT NEMECVCP (NO GOOD LEPTON KINEMATICS) ***')
      IERR=1
      RETURN

 800  WRITE(6,801)
 801  FORMAT(' *** ERROR AT NEMECVCP (NU WAS NOT NUMU) ***')
      IERR=1
      RETURN

 910  WRITE(6,911)
 911  FORMAT(' *** ERROR AT NEMECVCP (NO GOOD HADRON KINEMATICS) ***')
      IERR=1
      RETURN

      END
