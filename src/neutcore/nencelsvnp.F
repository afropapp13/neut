************************************************************************
*     --------------------------------------
      SUBROUTINE NENCELSVNP(IPAR,E,MODE,DNEUT,IERR)
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
*       MODE   : Interaction mode
*       ENE    : NEUTRINO ENERGY ( GEV )
*       DNEUT(3) : DIRECTION OF NEUTRINO
*     
*     (Output)
*       IERR    : ERROR CODE      
*       COMMON /NEWORK/
*
*     (Creation Date and Author)
*       2011.03.15 ; Y.Hayato
*     
************************************************************************
      IMPLICIT NONE

#include "nework.h"
#include "neutmodel.h"
#include "necard.h"

      INTEGER*4 ipar,mode
      real*4    e,cost,theta,phi,eta,dum
      real*4    abspnnc,pinnuc,einnuc,am,einit,p4in,eoutl
      real*4    eoutp,efinal,p4out

      integer*4 ierr,ii,icall,ierr2,iii,i,j
      integer*4 itarg
      
C      real*4    rlu,fnnucl
C      external  rlu,fnnucl
      real*4    rlu,fnnpot
      external  rlu,fnnpot

      REAL*4  DNEUT(3),PLEP(3),PF(3),PBUF0(3),PBUF1(3),PBUF2(3)
      INTEGER IORG(4)/ 0, 0, 1, 2/
      INTEGER IFLG(4)/-1,-1, 0, 0/
      INTEGER ICRN(4)/ 0, 0, 1, 1/

      REAL*4 PNOUT(3)
      IERR = 0
C
      IF (abs(mode).eq.51) then
         ITARG = 2212
      ELSE IF (abs(mode).eq.52) then
         ITARG = 2112
      else
         write(*,*) 'NENCELSVNP: Invalid mode',MODE
         STOP
      ENDIF         
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
C
C --  SET FINAL LEPTON MOMENTUM
C      
         CALL NENCELSLPV(E,PF,IPAR,MODE,PLEP,PNOUT,IERR2)
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
            MODENE= MODE
            NUMNE = 4
            
            IPNE(1)=IPAR
            IPNE(2)=ITARG
            IPNE(3)=IPAR
            IPNE(4)=ITARG

            DO 30 I=1,4
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
            EOUTL =sqrt(PNE(1,3)**2+PNE(2,3)**2+PNE(3,3)**2)

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
