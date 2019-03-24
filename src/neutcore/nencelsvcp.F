************************************************************************
*     --------------------------------------
      SUBROUTINE NEELSVCP(IPAR,E,DNEUT,IERR)
*     --------------------------------------
*
*     (Purpose)
*       Vector generation for elastic event
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
*       1984.01.11 ; M.Nakahata
*       1987.06.23 ; N.SATO     INCLUDE TAU
*       1995.02.17 ; K. KANEYUKI  FOR S.K.
*                    ELSVCP -> NEELSVCP
*       2001.03.03 ; Y.Hayato  USE NEELSLPV (instead of NEELSLVC)
*     
************************************************************************
      IMPLICIT NONE

#include "nework.h"
#include "necard.h"

      INTEGER*4 ipar
      real*4    e,cost,theta,phi,eta,dum
      real*4    abspnnc,pinnuc,einnuc,am,einit,p4in,eoutl
      real*4    eoutp,efinal,p4out,amlep

      integer*4 ierr,ii,icall,ierr2,iii,i,j
      
C      real*4    rlu,fnnucl
C      external  rlu,fnnucl
      real*4    rlu,fnnpot
      external  rlu,fnnpot

      REAL*4  DNEUT(3),PLEP(3),PF(3),PBUF0(3),PBUF1(3),PBUF2(3)
      INTEGER IPMODE(4,6)/ 12, 2112, 11, 2212, 14, 2112, 13, 2212,
     &                     16, 2112, 15, 2212,-12, 2212,-11, 2112,
     &                    -14, 2212,-13, 2112,-16, 2212,-15, 2112/
      INTEGER IORG(4)/ 0, 0, 1, 2/
      INTEGER IFLG(4)/-1,-1, 0, 0/
      INTEGER ICRN(4)/ 0, 0, 1, 1/

      REAL*4 PNOUT(3)

      IERR = 0
      
      IF (IPAR.EQ. 12) II=1

      IF (IPAR.EQ. 14) II=2
      IF (IPAR.EQ. 16) II=3
      IF (IPAR.EQ.-12) II=4
      IF (IPAR.EQ.-14) II=5
      IF (IPAR.EQ.-16) II=6



C      
C --  SET NEUTRINO DIRECTION
C
C      CALL NERNNEUT(E,DNEUT,IPAR)
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
C --  SET NUCLEON MOMENTUM AND DIRECTION
C      
      DO 10 ICALL=0,100
C         CALL NEELSSPF(E,IPAR,PF,ICALL,IERR)
C         IF (IERR.NE.0) RETURN
C
C --  SET FINAL LEPTON MOMENTUM
C      
C         CALL NEELSLVC(E,PF,IPAR,PLEP,IERR2)
         CALL NEELSLPV(E,PF,IPAR,PLEP,PNOUT,IERR2)
         IF (IERR2.EQ.0) THEN
C
C     -- BACK TO INITIAL COORDINATE
C     - RANDOM ROTATION AROUND X-AXIS
C
            ETA=RLU(DUM)*6.283185
            DO 20 III=1,3
               IF (III.EQ.1) THEN
                  PBUF0(1)=PLEP(1)
                  PBUF0(2)=PLEP(2)
                  PBUF0(3)=PLEP(3)
               ELSE IF (III.eq.2) THEN
                  PBUF0(1)=PF(1)
                  PBUF0(2)=PF(2)
                  PBUF0(3)=PF(3)
               ELSE
                  PBUF0(1)=PNOUT(1)
                  PBUF0(2)=PNOUT(2)
                  PBUF0(3)=PNOUT(3)
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
                  PF(1)=PBUF0(1)
                  PF(2)=PBUF0(2)
                  PF(3)=PBUF0(3)
               ELSE
                  PNOUT(1)=PBUF0(1)
                  PNOUT(2)=PBUF0(2)
                  PNOUT(3)=PBUF0(3)
               ENDIF
 20         CONTINUE
C     
C     -- NOW STORE VECTOR
C     
            IF ( IPAR.GT.0 ) MODENE= 1
            IF ( IPAR.LT.0 ) MODENE=-1
            NUMNE = 4
            DO 30 I=1,4
               IPNE(I)=IPMODE(I,II)
               IORGNE(I)=IORG(I)
               IFLGNE(I)=IFLG(I)
               ICRNNE(I)=ICRN(I)
               DO 40 J=1,3
                  IF (I.EQ.1) PNE(J,I)=E*DNEUT(J)
                  IF (I.EQ.2) PNE(J,I)=PF(J)
                  IF (I.EQ.3) PNE(J,I)=PLEP(J)
                  IF (I.EQ.4) THEN
                     PNE(J,I)=PNE(J,1)+PNE(J,2)-PNE(J,3)
                     PNE(J,I)=PNOUT(J)
                  ENDIF
 40            CONTINUE
 30         CONTINUE

            CALL MCMASS(IPNE(3),AMLEP)
            AMLEP=AMLEP*1.0E-3
            CALL MCMASS(IPNE(4),AM)
            AM=AM*1.0E-3

            ABSPNNC =sqrt(PNE(1,4)**2+PNE(2,4)**2+PNE(3,4)**2)

C            write(*,*) 'PNABS=',PNABS,'(q^2 only=',ABSPNNC,')'
C            DO 60 I=1,3
C               PNE(I,4)=PNE(I,4)*(PNABS/ABSPNNC)
C 60         continue
C            write(*,*) 'Corr :',
C     $           sqrt(PNE(1,4)**2+PNE(2,4)**2+PNE(3,4)**2)
C
C --        NORMAL ENDING            
C            
            if (QUIET.eq.0) then
               write(*,'(A10,F12.5,A1,F12.5,A1,F12.5,A1,F12.5)') 
     $              'INLEP =(',PNE(1,1),',',PNE(2,1),',',PNE(3,1),')'
               write(*,'(A10,F12.5,A1,F12.5,A1,F12.5,A1,F12.5)') 
     $              'INNUC =(',PNE(1,2),',',PNE(2,2),',',PNE(3,2),')'
            endif
            
            PINNUC=sqrt(PNE(1,2)**2+PNE(2,2)**2+PNE(3,2)**2)
            EINNUC=sqrt(PNE(1,2)**2+PNE(2,2)**2+PNE(3,2)**2+AM**2)
            Einit= EINNUC+E
            p4in =(EINNUC+E)**2-
     $            ( (PNE(1,2)+PNE(1,1))**2+(PNE(2,2)+PNE(2,1))**2
     $             +(PNE(3,2)+PNE(3,1))**2)

            if (QUIET.eq.0) then
               write(*,'(A10,F12.5,A1,F12.5,A1,F12.5,A1,F12.5)') 
     $              'OUTLEP=(',PNE(1,3),',',PNE(2,3),',',PNE(3,3),')'
            endif
            EOUTL =sqrt(PNE(1,3)**2+PNE(2,3)**2+PNE(3,3)**2
     $             +AMLEP**2)

            if (QUIET.eq.0) then
               write(*,'(A10,F12.5,A1,F12.5,A1,F12.5,A1,F12.5)') 
     $              'OUTNUC=(',PNE(1,4),',',PNE(2,4),',',PNE(3,4),')'
            endif
            EOUTP =sqrt(PNE(1,4)**2+PNE(2,4)**2+PNE(3,4)**2
     $           +AM**2)

            EFINAL=EOUTL+EOUTP
            p4out = (EOUTL+EOUTP)**2-
     $              ( (PNE(1,3)+PNE(1,4))**2+(PNE(2,3)+PNE(2,4))**2
     $               +(PNE(3,3)+PNE(3,4))**2)

            if (QUIET.eq.0) then
               write(*,'(A10,F12.5,A10,F12.5)') 
     $              'Ein : ',Einit*1000., 'Eout : ',EFINAL*1000.
               
               write(*,'(A10,F12.5,A10,F12.5,A10,F12.5)') 
     $              'Ediff : ',(Einit-Efinal)*1000., 
     $              'VnuclD:' ,(FNNPOT(2,ABSPNNC)-FNNPOT(1,PINNUC))*1000.
C     $           'VnuclD:' ,(FNNUCL(ABSPNNC)-FNNUCL(PINNUC))*1000.
C     $           'Vnucl :' ,(FNNUCL(ABSPNNC))*1000.,
               
               write(*,'(A10,F12.5,A10,F12.5)') 
     $              'P4in : ',P4in, 'P4out : ',P4out
            endif

            RETURN
         ENDIF
 10   CONTINUE
C     
C --   ERROR RETURN
C     
      WRITE(6,900)
 900  FORMAT(' *** ERROR AT ELSVCP (TOO MANY ELSLVC CALL) ***')
      IERR=1
      RETURN
      END
