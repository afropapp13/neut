************************************************************************
*     -------------------------
      REAL*8 FUNCTION DNELSNL(Q2)
*     -------------------------
*     
*     (Purpose)
*     Calculate differential elastic cross section d(sigma)/d(q2)
*     with pauli principle      
*     
*     (Input)
*     Q2     : Q**2
*     
*     (Output)
*     DNELSNL: differential elastic cross section(new vers.)
*     
*     (Creation Date and Author)
*     1987.06.26 ; N.SATO  PRD.
*     1988.10.07 ; T.KAJITA
*     EFFECT OF NUCLEON BINDING ENERGY WAS ADDITIONALLY
*     CONSIDERED
*     1995.02.06 ; K. KANEYUKI  FOR S.K.
*     FUNELS -> DNELS
*     2000.09.05 ; K.Kaneyuki COMMENT      IF(KP.GE.3) GO TO 220
*     2000.09.14 ; K.Kaneyuki FNNUCL(PFABS) -> FNNUCLPF
*     2001.04.24 ; Y.Hayato (Use LAB Frame)
*     In this code, PNEU and DNEU are in LAB frame
*     
************************************************************************
      IMPLICIT NONE
      
      REAL*8   Q2

#include "neutparams.h"

      REAL*8   DNELSQ2
      EXTERNAL DNELSQ2
C      real*4   FNNUCL
C      EXTERNAL FNNUCL
      real*4   FNNPOT
      EXTERNAL FNNPOT

      REAL*8   EINP

      INTEGER*4 IPR
      REAL*4 PINNUC(2),APINNUC,EINNUC,EINNEU,PTOTLI(2),PTOTLIS,
     $     ENSTOP,FNNUCLPF,AMIN,AMOUT,AMLEP,PFERMI,
     $     ESUM,AMLEP2,AMOUT2
      COMMON /NEINTG3/ PINNUC,APINNUC,EINNUC,EINNEU,PTOTLI,PTOTLIS,
     $     ENSTOP,IPR,FNNUCLPF,AMIN,AMOUT,AMLEP,PFERMI,
     $     ESUM,AMLEP2,AMOUT2

      REAL*4    PNFABS,pnfabs2
      REAL*8    ptotliy2

      INTEGER ILOOPMX,IPHIMX,ILOOP,IPHI,ISIGN
      PARAMETER (ILOOPMX=4)
      PARAMETER (IPHIMX=10)
      REAL*8  DPHI

      integer*4 incctr, itotctr
      real*8    dpassctr

      real*8    cosphi, cosphi2, dsign

      real*8    repa,repb,repc,repd,repf,repg,reph
      real*8    repadve, repadve2
      real*8    elep,  enf
      REAL*4    PLEPF(3)
      
      real*8    dfact

      integer*4 ierrctr
      logical first /.true./

      if (first) then
         DPHI=3.1415926535/dble(IPHIMX*2)
         first = .false.
      endif

      DNELSNL = 0.0D0
      DFACT   = 0.0D0
      INCCTR  = 0
      ITOTCTR = 0
      DPASSCTR= 0.0D0

      DO 1250 IPHI=0,IPHIMX-1
         IF (IPHI.ne.0) THEN
            COSPHI =COS(DPHI*dble(IPHI))
         else
            COSPHI =0.99999999
         ENDIF

         COSPHI2=COSPHI**2

         REPA=-1.*q2-AMLEP2

         REPADVE=REPA/EINNEU
         REPADVE2=(REPA/EINNEU)**2
         
         DSIGN=-1.

         DO 1200 ISIGN=0,1
            DSIGN=DSIGN*(-1.)
            
            PNFABS=APINNUC

            PTOTLIY2=(PTOTLI(2)**2)*COSPHI2

            ierrctr=0
            DO 1510 ILOOP=1,ILOOPMX

C               REPB=ESUM-FNNUCL(PNFABS)
               REPB=ESUM-FNNPOT(2,PNFABS)
#ifdef DEBUG
C               IF (abs(REPB-(EINNEU+EINNUC+FNNUCLPF-FNNUCL(PNFABS)))
               IF (abs(REPB-(EINNEU+EINNUC+FNNUCLPF-FNNPOT(2,PNFABS)))
     $              .gt.1.E-6) THEN
                  write(*,*) "REPB has prescision problem."
               endif
#endif

               REPC=REPB**2-PTOTLIS+AMLEP2-AMOUT2+(PTOTLI(1)*REPADVE)
#ifdef DEBUG
               if (abs( repc - (REPB**2
     $              -( (EINNEU+PINNUC(1))**2+PINNUC(2)**2)
     $              +AMLEP*AMLEP
     $              -AMOUT*AMOUT
     $              +(PTOTLI(1)*REPA/EINNEU)))
     $              .gt.1.E-6) THEN
                 write(*,*) "REPC has precision problem."
               endif
#endif

               REPD=2*(REPB-PTOTLI(1))
#ifdef DEBUG
               if (abs(REPD-2*(REPB-EINNEU-PINNUC(1)))
     $              .gt.1.e-6)  THEN
                  write(*,*) "REPD has precision problem."                  
               endif
#endif

               REPF=REPC**2+(PTOTLIY2*REPADVE2)+(4*AMLEP2*PTOTLIY2)
#ifdef DEBUG
               if (abs(REPF -( REPC**2
     $              +( PINNUC(2)*PINNUC(2)
     $              *COSPHI*COSPHI
     $              *REPA*REPA
     $              /EINNEU/EINNEU)
     $              +(4*AMLEP*AMLEP*PINNUC(2)*PINNUC(2)*COSPHI*COSPHI)))
     $              .gt.1.e-6)  THEN
                  write(*,*) "REPF has precision problem."
               endif
#endif

               REPG=(4*REPADVE*PTOTLIY2)-(2*REPC*REPD)
#ifdef DEBUG
              if (abs(REPG-
     $              ((4*REPA*PINNUC(2)*PINNUC(2)*COSPHI*COSPHI/EINNEU)
     $              -(2*REPC*REPD)))
     $              .gt.1.e-6)  THEN
                 write(*,*) "REPG has precision problem."
              endif
#endif
               
               REPH=(REPG**2)-(4*REPD**2*REPF)

#ifdef DEBUG
               write(*,'(A,F12.7,A,F12.7,A,F12.7,A,F12.7)') 
     $              'A=',REPA,', B=',REPB,', C=',REPC,', D=',REPD
               write(*,'(A,F12.7,A,F12.7,A,F12.7)') 
     $              'F=',REPF,', G=',REPG,', H=',REPH
#endif

               IF (ABS(REPH).lt.1.D-8) REPH=0.D0
               
               IF (REPH.LT.0.D0) THEN
C                  write(*,*) 'value in SQRT lt.0(REPH=',REPH,')'
                  IF (ILOOP.eq.1)       goto 1200
                  IF (ILOOP.eq.ILOOPMX) THEN
C                     write(*,*) 'value in SQRT lt.0(REPH=',REPH,')'
                     ierrctr=1
                     goto 1190
                  ENDIF

                  REPH=0.
                  
               endif

 1515          ELEP=(-1.*REPG+DSIGN*dsqrt(REPH))/(2*REPD**2)

               IF (ELEP.LT.AMLEP) THEN
                  IF (ILOOP.eq.1)       goto 1200
                  if (ILOOP.eq.ILOOPMX) THEN
#ifdef DEBUG
                     write(*,*) 
     $                    'Energy of Final LEPTON lt.MLEP',
     $                    ' (ELEP=',ELEP,')'
#endif
                     ierrctr=1
                     goto 1190
                  ENDIF
                  ENF=REPB-AMLEP
                  IF (ENF.LT.AMOUT) THEN
                     ENF=AMOUT+0.0001
                  ENDIF
                  goto 1575
               endif
               
               PLEPF(1)=REPADVE/2.+ELEP
               if ((elep**2-plepf(1)**2-amlep2).le.0.E0) THEN
                  IF (ILOOP.eq.1)       goto 1200
                  if (ILOOP.eq.ILOOPMX) THEN
#ifdef DEBUG
                     write(*,*) 
     $                    'Energy of Final LEPTON lt.MLEP',
     $                    ' (ELEP=',ELEP,')'
#endif
                     ierrctr=1
                     goto 1190
                  ENDIF
                  ENF=REPB-AMLEP
                  IF (ENF.LT.AMOUT) THEN
                     ENF=AMOUT+0.0001
                  ENDIF
                  goto 1575
               ENDIF

               IF (abs(PTOTLI(2)).gt.0.00001) THEN
                  PLEPF(2)=(REPD*ELEP-REPC)/(2*PTOTLI(2))
#ifdef DEBUG
                  if (     abs( plepf(2)
     $                         -sqrt(elep**2-plepf(1)**2-AMLEP2)
     $                          *cosphi*rsign)
     $                 .gt.0.005) then
                     Write(*,*) "INCONSISTENT?"
                     write(*,*) "PLEPF(2) =",PLEPF(2)
                     write(*,*) "PLEPF(2)'=",
     $                    sqrt(elep**2-plepf(1)**2-AMLEP2)*cosphi
                     write(*,*) "sqrt(elep**2-plepf(1)**2-AMLEP2)=",
     $                    sqrt(elep**2-plepf(1)**2-AMLEP2)
                  endif
#endif
               else
                  PLEPF(2)= (-REPADVE2/4.-REPADVE*ELEP-AMLEP2)*COSPHI2
                  IF (PLEPF(2).gt.0.) THEN
                     PLEPF(2)=sqrt(PLEPF(2))*DSIGN
                  else
                     PLEPF(2)=0.
                  ENDIF
               endif

               IF (     PLEPF(1)**2+PLEPF(2)**2+AMLEP2-elep*elep
     $              .gt.0.0001)   THEN
                  IF (ILOOP.eq.1)       goto 1200
                  if (ILOOP.eq.ILOOPMX) THEN
#ifdef DEBUG
                    write(*,*) 
     $                    'kinematical boundary of lepton failed'
#endif
                     ierrctr=1
                     goto 1190
                  ENDIF
                  ENF=REPB-ELEP                  
                  IF (ENF.LT.AMOUT) THEN
                     ENF=AMOUT+0.0001
                  ENDIF
                  goto 1575

#ifdef DEBUG
               else
                  PLEPF232=-1.*( REPADVE2/4.
     $                 +REPADVE*ELEP
     $                 +AMLEP2)
                  
                  IF (PLEPF232.le.0.) THEN
                     IF (ILOOP.eq.1) goto 1200
                     IF (ILOOP.eq.ILOOPMX) then
                        write(*,*) 
     $                       'kinematical boundary of lepton failed'
                        ierrctr=1
                        goto 1190
                     endif
                     PLEPF232=0.
                  endif
                  
                  IF (    abs(sqrt(PLEPF232)*COSPHI*DSIGN-PLEPF(2))
     $                 .gt.0.010) THEN
                     
                     write(*,*) "PTOTLI(2)=",PTOTLI(2),'cosphi=',cosphi
                     write(*,'(A10,F12.8,A,F12.8,A,F12.8,A)') 
     $                    'PLEPF(2)=(',PLEPF(2),',',
     $                    sqrt(PLEPF232)*COSPHI*DSIGN,'(plepf23^2)=',
     $                    plepf232,')'
                  endif
#endif
                  
               endif

               ENF=REPB-ELEP

               IF (ENF**2.le.AMOUT2) THEN
#ifdef DEBUG
                  write(*,*) 
     $                 'Insufficient energy of final state nucleon',
     $                 'ENF=',ENF,')'
#endif

                  IF (ILOOP.eq.1)       goto 1200
                  if (ILOOP.eq.ILOOPMX) THEN
#ifdef DEBUG
                     write(*,*) 
     $                    'Insufficient energy of final state nucleon',
     $                    'ENF=',ENF,')'
#endif
                     ierrctr=1
                     goto 1190
                  endif

                  ENF=AMOUT+0.0001
                  goto 1575
               endif


 1575          PNFABS2=ENF**2-AMOUT2
               PNFABS = sqrt(PNFABS2)

#ifdef DEBUG
              write(*,'(A10,F12.8,A,F12.8)') 'q2(lep)=(',
     $              -1*( (ELEP-ENEU)**2-(PLEPF(1)-ENEU)**2
     $                  -(PLEPF(2)-0.)**2-(PLEPF(3)-0.)**2)
     $              ,')       : should be ',q2
     
              write(*,'(A10,F12.8,A10,F12.8)') 
     $              'Ein : ',Einit, 'Eout : ',EFINAL
     
              write(*,'(A10,F12.3,A10,F12.3)') 
     $             'Ediff : ',(Einit-Efinal)*1000., 
     $             'Vnucl: ' ,(FNNPOT(2,PNFABS)-FNNUCLPF)*1000.
C     $             'Vnucl: ' ,(FNNUCL(PNFABS)-FNNUCLPF)*1000.
#endif

 1510       CONTINUE
            
            IF(PNFABS.LT.PFSURF) GOTO 1190

C           IF ((IPHI.eq.0).or.(IPHI.eq.10)) THEN
C               IPASSCTR=IPASSCTR+1
#ifndef RECOIL
            DPASSCTR=DPASSCTR+1.D0
#else
            DPASSCTR= DPASSCTR
     $               + (2*AMIN*(ELEP**2-AMLEP2))
     $                /( 2*AMIN*(ELEP**2-AMLEP2)
     $                  -(2*EINNEU*AMLEP2)-(ELEP*REPA))
#endif
C            else
C              IPASSCTR=IPASSCTR+2               
C            ENDIF

#ifdef DEBUG
C           write(*,'(A10,F12.8,A,F12.8,A,F12.8,A,F12.8)') 
C     $           'OUTLEP=(',PLEPF(1),',',PLEPF(2),',',PLEPF(3),
C     $           ')/Eoutlep=',elep
C           write(*,'(A10,F12.8,A,F12.8,A,F12.8,A,F12.8)')
C     $           'OUTNUC=(',PNF(1),  ',',PNF(2),  ',',PNF(3),
C     $           ')/Eoutnuc=',enf
C           write(*,'(A10,F12.8,A,F12.8)') 'q2(lep)=(',
C     $           -1*( (ELEP-ENEU)**2-(PLEPF(1)-ENEU)**2
C     $           -(PLEPF(2)-0.)**2-(PLEPF(3)-0.)**2)
C     $           ,')       : should be ',q2
C     
C           write(*,'(A10,F12.8,A10,F12.8)') 
C     $           'Ein : ',Einit, 'Eout : ',EFINAL
C     
C           write(*,'(A10,F12.3,A10,F12.3)') 
C     $           'Ediff : ',(Einit-Efinal)*1000., 
C     $           'Vnucl: ' ,(FNNPOT(2,PNFABS)-FNNUCLPF)*1000.
CC     $           'Vnucl: ' ,(FNNUCL(PNFABS)-FNNUCLPF)*1000.
C     
C           write(*,'(A10,F12.3,A10)') 
C     $           'Econserv: ',
CC     $           (EINIT+FNNUCLPF-(EFINAL+FNNUCL(PNFABS)))*1000,
C     $           (EINIT+FNNUCLPF-(EFINAL+FNNPOT(2,PNFABS)))*1000,
C     $           'MeV'
C     
CC           IF ((abs(EINIT+FNNUCLPF-(EFINAL+FNNUCL(PNFABS)))*1000
C           IF ((abs(EINIT+FNNUCLPF-(EFINAL+FNNNPOT(2,PNFABS)))*1000
C     $           .GT.2)) THEN
C              write(*,'(A,I8,A)')
C     $                 'EDIFF > 2.0MeV (COUNT=',ICOUNT,')'
C              ICOUNT = ICOUNT+1
C           ENDIF
#endif

C     1190       IF ((IPHI.eq.0).or.(IPHI.eq.10)) THEN
 1190       CONTINUE
            ITOTCTR=ITOTCTR+1
C          else
C             ITOTCTR=ITOTCTR+2               
C          ENDIF
            

#ifdef DEBUG
            if (ierrctr.eq.0) then
               IF (   (sqrt(PLEPF(1)**2+PLEPF(2)**2+AMLEP2)+0.001)
     $              .lt.ELEP)   THEN
                  
                  PLEPF(3)=sqrt( (ELEP**2-PLEPF(1)**2-AMLEP2)
     $                          *(1-COSPHI2))
               
                  PNF(1)=PTOTLI(1)-PLEPF(1)
                  PNF(2)=PTOTLI(2)-PLEPF(2)
                  PNF(3)=0        -PLEPF(3)
               
                  ENF2  = PNF(1)**2+PNF(2)**2+PNF(3)**2+AMOUT2
                  
                  IF (abs(sqrt(ENF2)-ENF).gt.0.001) THEN
                  
                     write(*,'(A10,F12.8,A,F12.8,A,F12.8,A,F12.8)') 
     $                    'OUTLEP=(',PLEPF(1),',',PLEPF(2),',',
     $                    PLEPF(3),')/Eoutlep=',elep
                  
                     write(*,'(A10,F12.8,A,F12.8,A,F12.8,A)')
     $                    'TOTIN=(',PTOTLI(1),  ',',PTOTLI(2),',',
     $                    0,')'
                  
                     write(*,'(A10,F12.8,A,F12.8,A,F12.8,A,F12.8,A)')
     $                    'OUTNUC=(',PNF(1),  ',',PNF(2),  ',',PNF(3),
     $                    ')/Poutnuc=',sqrt(ENF2-AMOUT2),')'
                  
                     write(*,'(A10,F12.8,A,F12.8,A)')
     $                    'Eoutnuc=(' , sqrt(ENF2) , ',' , ENF , ')'
                  endif
               endif
            endif
#endif
 1200    continue

 1250 continue

      if (ITOTCTR.eq.0)  then
         return
      endif
      EINP = ENSTOP
      DFACT= DPASSCTR/dble(ITOTCTR)

      DNELSNL=DNELSQ2(EINP,IPR,Q2)*DFACT
      
      RETURN
      END
      
