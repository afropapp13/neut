************************************************************************
*     -------------------------------------------
      SUBROUTINE NEELSLPV(ENEU,PF,IPAR,PLEP,PNUC,IERR)
*     -------------------------------------------
*
*     (Purpose)
*       set final lepton momentum vector for elastic scattering
*
*     (Input)
*       ENE    : NEUTRINO ENERGY ( GEV )
*       PF(3)  : FERMI MOMENTUM ( GEV/C )      
*       IPAR   : NEUTRINO TYPE
*                  12 : NEU-E
*                 -12 : NEU-E-BAR            
*                  14 : NEU-MU
*                 -14 : NEU-MU-BAR            
*                  16 : NEU-TAU
*                 -16 : NEU-TAU-BAR
*     
*     (Output)
*       PLEP(3) : FINAL LEPTON MOMENTUM ( GEV/C )
*       IERR    : ERROR CODE      
*
*     (Creation Date and Author)
*       1984.01.08 ; M.Nakahata
*       1987.06.23 ; N.SATO     INCLUDE TAU
*       1988.11.28 ; T.KAJITA   NUCLEAR VINDING ENERGY IS CONSIDERD
*       1995.01.30 ; K. KANEYUKI  FOR S.K.
*                    ELSLVC -> NEELSLVC
*       1996.01.11 ; K.Kaneyuki  NEPAUFLG=1 -> No Pauli Effect
*     
*       1998.02.23 ; K.Ishihara set PFERMI 217 MeV
*       2000.09.05 ; K.Kaneyuki COMMENT OUT   IF(KP.GE.3) GO TO 220
*       2000.09.14 ; K.Kaneyuki FNNUCL(PFABS) -> FNNUCLPF
*       2001.03.03 ; Y.Hayato   create neelslvp.F to return P_nucl
*       2013.11.19 ; A.Redij RPA correction added to CCQE cross-section
*
************************************************************************
      IMPLICIT NONE
#include "necard.h"
#include "neutparams.h"
      REAL    ENEU,AMIN,AMOUT
      INTEGER IPAR,IERR,iplep,i,ilast,loop

      SAVE    ICOUNT,IEVENT
      INTEGER ICOUNT,IEVENT
      
      REAL    pfabs,e,beta,gm,enstop,q2maxx,dq2,amlep,bmax,blast,q2
      real    el,pl,cost,amp,dum,hit,sigma
C      real    phi,pabs,eabs
      
C      real    fnq2max,fnnucl,rlu
C      external fnq2max,fnnucl,rlu
      real    fnq2max,fnnpot,rlu, fnels_rpascl
      external fnq2max,fnnpot,rlu, fnels_rpascl 

      REAL PF(3),PNEU(3),EV(3),DNEU(3),PLEP(3),ENF
C      real      DLEP(3),PNP(3)
      REAL    pneulab(3),pnulabi(3),ptotli(3),plepf(3),pnf(3)
      REAL    einnuc,einit,elep,efinal,ptotlis,pnfabs,pnfabs2
      REAL    phiang,cosphi,ptotliy2
      REAL*8  repa,repb,repc,repd,repf,repg,reph
      REAL*8  DELEP,DTMPY,DENEU,DAMLEP
      REAL*8  DPLEPF(3)
      
      integer*4  isub,iloop
      REAL PNUC(3)
      REAL XMAP(51,11)
      REAL*8 DNELSQ2, DBLE0, DBLQ2
      EXTERNAL DNELSQ2
      REAL FNNUCLPF

      real*8 DSGN
      DATA    ICOUNT/0/
      DATA    IEVENT/1/

C     
      IERR=0
      DSGN=0.D0
      DENEU=ENEU
      
      IF (IPAR.gt.0) then
         CALL MCMASS(2112,AMIN)
         CALL MCMASS(2212,AMOUT)
      ELSE
         CALL MCMASS(2212,AMIN)
         CALL MCMASS(2112,AMOUT)
      ENDIF         

      AMIN=AMIN*1.E-3
      AMOUT=AMOUT*1.E-3

C      write(*,*) '---------------------------------------'
C      write(*,'(A,I8,A,I8)') 'Event #',ievent,'  INPAR:',ipar
C      ievent = ievent+1
C
C -- CALCUALTE NEUTRINO ENERGY IN NUCELON STOP SYSTEM
C

C      DO 100 LOOP=0,1000000
C--reduce loop to 1000times.
      DO 100 LOOP=0,1000

         PNEU(1)=ENEU
         PNEU(2)=0.
         PNEU(3)=0.

C-- set pf (not considering pauli blocking or potential)
         call neelsspff(eneu,ipar,pf,LOOP,ierr)
         if (ierr.ne.0) goto 9000

         PFABS=SQRT(PF(1)**2+PF(2)**2+PF(3)**2)
         E=SQRT(PFABS**2+AMIN**2)
         BETA=PFABS/E
         GM=1./SQRT(1.-BETA**2)
         IF (PFABS.NE.0.) THEN
            EV(1)=-PF(1)/PFABS
            EV(2)=-PF(2)/PFABS
            EV(3)=-PF(3)/PFABS
         ELSE
            EV(1)=1.
            EV(2)=0.
            EV(3)=0.
         ENDIF
         CALL MCVECBST(PNEU,0.,EV,GM)
         ENSTOP=SQRT(PNEU(1)**2+PNEU(2)**2+PNEU(3)**2)
         IF(ENSTOP.LT.0.001 .AND. IABS(IPAR).EQ.12)GO TO 9000
         IF(ENSTOP.LT.0.112 .AND. IABS(IPAR).EQ.14)GO TO 9000
         IF(ENSTOP.LT.3.48109 .AND. IABS(IPAR).EQ.16)GO TO 9000
         DNEU(1)=PNEU(1)/ENSTOP
         DNEU(2)=PNEU(2)/ENSTOP
         DNEU(3)=PNEU(3)/ENSTOP
C     
C     -- MAKE (Q**2 -- SCATTERING ANGLE) MAP
C     
         CALL VZERO(XMAP,561)

         EV(1)=-EV(1)
         EV(2)=-EV(2)
         EV(3)=-EV(3)
         Q2MAXX=FNQ2MAX(ENSTOP,IPAR)
         DQ2=Q2MAXX/50.
         IPLEP=IABS(IPAR)-1
         CALL MCMASS(IPLEP,AMLEP)
         AMLEP=AMLEP*1.0E-3
         DAMLEP=AMLEP
         BMAX=0.
         BLAST=0.E+0
         FNNUCLPF=FNNPOT(1,PFABS)

C For RPA correction to CCQE cross-section, QE amplitude is scaled by the ratio of (rpa corrected CCQE/CCQE)cross-section
C FNELSPAU_RPASCL(Enu,Q2,IPAR) returns this ratio given, neutrino Energy, Q**2 and neutrino type for repective target element.
C The function will return 1 when RPA effect is disabled through neut.card.
C
CCCCCC Here the imeplemetation is wrong as the E is in lab frame and Q2 in nucleon frame.
C
C     
C     --  SEARCH PEAK ROUGHLY
C     
         DO 10 I=0,50 
            Q2=DQ2*FLOAT(I)
            EL=ENSTOP-Q2/2./AMIN
            PL=SQRT(EL**2-AMLEP**2)
            COST=(2.*ENSTOP*EL-AMLEP**2-Q2)/2./ENSTOP/PL
            IF(ABS(COST).GT.1.)GO TO 10
            DBLE0=ENSTOP
            DBLQ2=Q2
            AMP=DNELSQ2(DBLE0,IPAR,DBLQ2)*FNELS_RPASCL(ENEU,Q2,IPAR)
            IF (BMAX.LT.AMP) BMAX=AMP
            IF (AMP.LT.BLAST) GOTO 20
            BLAST=AMP
 10      CONTINUE
 20      CONTINUE
         ILAST=I-1
C     
C     --  SEARCH PEAK CLOSELY
C     
         Q2=DQ2*FLOAT(ILAST)+DQ2*0.02
         EL=ENSTOP-Q2/2./AMIN
         PL=SQRT(EL**2-AMLEP**2)
         COST=(2.*ENSTOP*EL-AMLEP**2-Q2)/2./ENSTOP/PL
         IF(ABS(COST).GT.1.) THEN
            AMP=0.0E+0
         ELSE
            DBLE0=ENSTOP
            DBLQ2=Q2
            AMP=DNELSQ2(DBLE0,IPAR,DBLQ2)*FNELS_RPASCL(ENEU,Q2,IPAR)
         ENDIF
         IF (AMP.GT.BLAST) THEN
C     
C     --     FORWARD
C     
            BMAX=AMP
            BLAST=AMP
            DO 30 I=2,50
               Q2=DQ2*FLOAT(ILAST)+DQ2*0.02*FLOAT(I)
               EL=ENSTOP-Q2/2./AMIN
               PL=SQRT(EL**2-AMLEP**2)
               COST=(2.*ENSTOP*EL-AMLEP**2-Q2)/2./ENSTOP/PL
               IF(ABS(COST).GT.1.)GO TO 30
               DBLE0=ENSTOP
               DBLQ2=Q2
               AMP=DNELSQ2(DBLE0,IPAR,DBLQ2)*FNELS_RPASCL(ENEU,Q2,IPAR)
               IF (BMAX.LT.AMP) BMAX=AMP
               IF (AMP.LT.BLAST) GOTO 50
               BLAST=AMP
 30         CONTINUE
         ELSE IF (ILAST.NE.0) THEN
C     
C     --     BACKWORD
C     
            DO 40 I=1,50
               Q2=DQ2*FLOAT(ILAST)-DQ2*0.02*FLOAT(I)
               EL=ENSTOP-Q2/2./AMIN
               PL=SQRT(EL**2-AMLEP**2)
               COST=(2.*ENSTOP*EL-AMLEP**2-Q2)/2./ENSTOP/PL
               IF(ABS(COST).GT.1.)GO TO 40
               DBLE0=ENSTOP
               DBLQ2=Q2
               AMP=DNELSQ2(DBLE0,IPAR,DBLQ2)*FNELS_RPASCL(ENEU,Q2,IPAR)
               IF (BMAX.LT.AMP) BMAX=AMP
               IF (AMP.LT.BLAST) GOTO 50
               BLAST=AMP
 40         CONTINUE
         ENDIF
 50      CONTINUE
C     
C     --  THROW DICE
C     
         EINNUC=sqrt(PFABS**2+AMIN**2)
         Einit=EINNUC + ENEU 

#ifdef DEBUG
         write(*,'(A10,F12.8,A,F12.8,A,F12.8,A,F12.8)')
     $        'INLEP =(',ENEU,   ',',0,      ',',0    ,  
     $        ')/Einlep =',eneu
         write(*,'(A10,F12.8,A,F12.8,A,F12.8,A,F12.8)') 
     $        'INNUC =(',PF(1),  ',',PF(2),  ',',PF(3),  
     $        ')/Einnuc =',einnuc
         write(*,'(A10,F12.8,A1)')
     $        'q2=(',q2,')'
#endif

C     DO 100 LOOP=1,1000000
         Q2=RLU(DUM)*Q2MAXX
         DBLQ2=Q2
         HIT=RLU(DUM)*BMAX

         
         SIGMA=DNELSQ2(DBLE0,IPAR,DBLQ2)*FNELS_RPASCL(ENEU,Q2,IPAR)
         IF (SIGMA.GT.HIT) THEN
            
            PHIANG=6.283*RLU(DUM)
            cosphi=COS(PHIANG)
            if (COSPHI.gt.0) THEN
               DSGN=1.D0
            ELSE
               DSGN=-1.D0
            ENDIF

            
            PNEULAB(1)=ENEU
            PNEULAB(2)=0.
            PNEULAB(3)=0.
            DO 1500 ISUB=1,3
               PNULABI(ISUB)=PF(ISUB)
 1500       continue

C---  for the comparison --- 
C 1200       EL=ENSTOP-Q2/2./AMIN
C            PL=SQRT(EL**2-AMLEP**2)
C            COST=(2.*ENSTOP*EL-AMLEP**2-Q2)/2./ENSTOP/PL
C     
C            DLEP(1)=COST*DNEU(1)-sqrt(1.-COST**2)*DNEU(2)
C            DLEP(2)=sqrt(1.-COST**2)*DNEU(1)+COST*DNEU(2)
C            DLEP(3)=0.
C     
C            do 1210 JJ=1,3
C               PLEP(JJ)=PL*DLEP(JJ)
C 1210       continue
C     
C            CALL MCVECBST(PLEP,AMLEP,EV,GM)
C            el = sqrt(plep(1)**2+plep(2)**2+plep(3)**2+AMLEP**2)
C            write(*,'(A10,F12.8,A,F12.8,A,F12.8,A,F12.8)') 
C    $           'OUTLEPX=(',PLEP(1),',',PLEP(2),',',
C    $           PLEP(3),')/Eoutlep=',el
C            write(*,'(A10,F12.8)')  'q^2=',
C    $           ( (el-eneu)**2
C    $            -(plep(1)-eneu)**2-plep(2)**2-plep(3)**2)


            PNULABI(1)=PF(1)
            PNULABI(2)=PF(2)
            PNULABI(3)=0.
            
            PTOTLI(1)=PNULABI(1)+ENEU
            PTOTLI(2)=PNULABI(2)
            PTOTLI(3)=PNULABI(3)
            PTOTLIS=(PTOTLI(1)**2+PTOTLI(2)**2+PTOTLI(3)**2)

            PTOTLIY2=(PTOTLI(2)**2)*COSPHI**2

            PNFABS=PFABS
            DO 1510 ILOOP=1,5

               REPA=-1.*q2-AMLEP**2
C               REPB=dble(ENEU)+dble(einnuc)
C     $              +(FNNUCLPF-FNNUCL(PNFABS))
               REPB=dble(ENEU)+dble(einnuc)
     $              +(FNNUCLPF-FNNPOT(2,PNFABS))
               REPC=REPB**2-dble(PTOTLIS)+dble(AMLEP)**2-dble(AMOUT)**2
     $              +dble(PTOTLI(1))*REPA/dble(ENEU)
               REPD=2*(REPB-PTOTLI(1))
               REPF= REPC**2+(PTOTLIY2)*(REPA**2)/(ENEU**2)
     $              +4*(AMLEP**2)*(PTOTLIY2)
               REPG=4*REPA*(PTOTLIY2)/ENEU-2*REPC*REPD
               REPH=REPG**2-4*REPD**2*REPF
C               write(*,'(A,F12.7,A,F12.7,A,F12.7,A,F12.7)') 
C     $              'A=',REPA,', B=',REPB,', C=',REPC,', D=',REPD
C               write(*,'(A,F12.7,A,F12.7,A,F12.7)') 
C     $              'F=',REPF,', G=',REPG,', H=',REPH

               IF (ABS(REPH).lt.1.D-8) REPH=0.D0

               IF (REPH.LT.0) THEN
                  if (QUIET.eq.0) then
                     write(*,*) 'neelslpv:ILOOP=',ILOOP,
     $                    'value in SQRT lt.0',
     $                    '(REPH=',REPH,')'
                  endif
                  IF ((ILOOP.eq.1).or.(ILOOP.eq.5)) then 
                     goto 100
                  endif
                  reph=0.
               endif

 1515          DELEP=(-1.*REPG+DSGN*dsqrt(REPH))/(2*REPD**2)

               ELEP = REAL(DELEP)
               IF (ELEP.LT.AMLEP) THEN
                  write(*,*) 
     $                 'neelslpv:ILOOP=',ILOOP,
     $                 'Energy of Final LEPTON lt.MLEP',
     $                 '(ELEP=',ELEP,')'
                  IF ((ILOOP.eq.1).or.(ILOOP.eq.5)) then 
                     goto 100
                  endif
C                  plepf(1)=0.
C                  plepf(2)=0.
C                  plepf(3)=0.
C                  goto 1575
                  ELEP=AMLEP
                  ENF=REPB-ELEP
                  PNFABS2=ENF**2-AMOUT**2
                  IF (PNFABS2.lt.0) PNFABS2=0.
                  goto 1585
               endif
               
               DPLEPF(1)=-q2/(2*DENEU)-(AMLEP**2)/(2*DENEU)+DELEP
               DTMPY = DELEP**2-AMLEP**2-DPLEPF(1)**2

               IF (DTMPY.LT.0) THEN
                  if (QUIET.eq.0) then
                     write(*,*) 
     $                    'neelslpv:Insufficient energy of LEPTON',
     $                    '(P_Y^2=',DTMPY,')'
                  endif 
                  IF ((ILOOP.eq.1).or.(ILOOP.eq.5)) then 
                     goto 100
                  endif
C                  plepf(1)=0.
C                  plepf(2)=0.
C                  plepf(3)=0.
C                  goto 1575
                  ENF=REPB-ELEP
                  PNFABS2=ENF**2-AMOUT**2
                  IF (PNFABS2.lt.0) PNFABS2=0.
                  goto 1585
               endif

               PLEPF(1)=DPLEPF(1)

               if (abs(PTOTLI(2)).gt.0.000001) then
                  PLEPF(2)=(REPD*ELEP-REPC)/(2*PTOTLI(2))
               else
                  PLEPF(2)= (-((REPA/ENEU)**2)/4.
     $                       -(REPA/ENEU)*ELEP
     $                       -AMLEP**2)*(COSPHI**2)
C                  PLEPF(2) = (DELEP**2-AMLEP**2-DPLEPF(1)**2)
C     $                      *(COSPHI**2)
                  IF (PLEPF(2).gt.0.) THEN
                     PLEPF(2)=sqrt(PLEPF(2))*DSGN
                  else
                     PLEPF(2)=0.
                  ENDIF
               endif
               IF (     PLEPF(1)**2+PLEPF(2)**2+AMLEP**2-elep*elep
     $              .gt.0.0001)   THEN
                  IF ((ILOOP.eq.1).or.(ILOOP.eq.5)) then 
                     goto 100
                  endif
C                  plepf(1)=0.
C                  plepf(2)=0.
C                  plepf(3)=0.
C                  goto 1575
                  ENF=REPB-ELEP
                  PNFABS2=ENF**2-AMOUT**2
                  IF (PNFABS2.lt.0) PNFABS2=0.
                  goto 1585
               endif

               PLEPF(3)=real(dsqrt(DTMPY))*SIN(PHIANG)
               
 1575          PNF(1)=ENEU+PF(1)-PLEPF(1)
               PNF(2)=PF(2)-PLEPF(2)
               PNF(3)=PF(3)-PLEPF(3)
               PNFABS2=(PNF(1)**2+PNF(2)**2+PNF(3)**2)

               ENF    = sqrt(PNFABS2+AMOUT**2)
               EFINAL = ELEP + ENF
 1585          PNFABS = sqrt(PNFABS2)

C              write(*,'(A10,F12.8,A,F12.8,A,F12.8,A,F12.8)') 
C    $              'OUTLEP=(',PLEPF(1),',',PLEPF(2),',',PLEPF(3),
C    $              ')/Eoutlep=',elep
C              write(*,'(A10,F12.8,A,F12.8,A,F12.8,A,F12.8)')
C    $              'OUTNUC=(',PNF(1),  ',',PNF(2),  ',',PNF(3),
C    $              ')/Eoutnuc=',enf

C              write(*,'(A10,F12.8,A,F12.8)') 'q2(lep)=(',
C    $              -1*( (ELEP-ENEU)**2-(PLEPF(1)-ENEU)**2
C    $                  -(PLEPF(2)-0.)**2-(PLEPF(3)-0.)**2)
C    $              ,')       : should be ',q2
C     
C              write(*,'(A10,F12.8,A10,F12.8)') 
C    $              'Ein : ',Einit, 'Eout : ',EFINAL
C     
C              write(*,'(A10,F12.3,A10,F12.3)') 
C     $             'Ediff : ',(Einit-Efinal)*1000., 
C     $             'Vnucl: ' ,(FNNPOT(2,PNFABS)-FNNUCLPF)*1000.
CC     $             'Vnucl: ' ,(FNNUCL(PNFABS)-FNNUCLPF)*1000.


 1510       CONTINUE
            
            IF(NEPAUFLG.EQ.0 .AND. PNFABS.LT.PFSURF) GOTO 100

            IF (PLEPF(1)/sqrt(ELEP**2-AMLEP**2).le.-0.998) THEN
               goto 100
            endif

#ifdef DEBUG
            write(*,'(A10,F12.8,A,F12.8,A,F12.8,A,F12.8)') 
     $           'OUTLEP=(',PLEPF(1),',',PLEPF(2),',',PLEPF(3),
     $           ')/Eoutlep=',elep
            write(*,'(A10,F12.8,A,F12.8,A,F12.8,A,F12.8)')
     $           'OUTNUC=(',PNF(1),  ',',PNF(2),  ',',PNF(3),
     $           ')/Eoutnuc=',enf
            write(*,'(A10,F12.8,A,F12.8)') 'q2(lep)=(',
     $           -1*( (ELEP-ENEU)**2-(PLEPF(1)-ENEU)**2
     $           -(PLEPF(2)-0.)**2-(PLEPF(3)-0.)**2)
     $           ,')       : should be ',q2
            
            write(*,'(A10,F12.8,A10,F12.8)') 
     $           'Ein : ',Einit, 'Eout : ',EFINAL
            
            write(*,'(A10,F12.3,A10,F12.3)') 
     $           'Ediff : ',(Einit-Efinal)*1000., 
     $           'Vnucl: ' ,(FNNPOT(2,PNFABS)-FNNUCLPF)*1000.
C     $           'Vnucl: ' ,(FNNUCL(PNFABS)-FNNUCLPF)*1000.
            
            write(*,'(A10,F12.3,A10)') 
     $           'Econserv: ',
     $           (EINIT+FNNUCLPF-(EFINAL+FNNPOT(2,PNFABS)))*1000,
     $           'MeV'
C     $           (EINIT+FNNUCLPF-(EFINAL+FNNUCL(PNFABS)))*1000,
C     $           'MeV'
            
#endif
C            IF ((abs(EINIT+FNNUCLPF-(EFINAL+FNNUCL(PNFABS)))*1000
            IF ((abs(EINIT+FNNUCLPF-(EFINAL+FNNPOT(2,PNFABS)))*1000
     $           .GT.2)) THEN
               write(*,'(A,I8,A)')
     $              'EDIFF > 2.0MeV (COUNT=',ICOUNT,')'
               ICOUNT = ICOUNT+1
            ENDIF

            PLEP(1)=PLEPF(1)
            PLEP(2)=PLEPF(2)
            PLEP(3)=PLEPF(3)

            PNUC(1)=PNF(1)
            PNUC(2)=PNF(2)
            PNUC(3)=PNF(3)



C           write(*,*) '---------------------------------------'
C           write(*,'(A10,F12.8,A,F12.8,A,F12.8,A,F12.8)') 
C    $           'OUTLEP=(',PLEP(1),',',PLEP(2),',',PLEP(3),
C    $           ')/Eoutlep=',elep
C           write(*,'(A10,F12.8,A,F12.8,A,F12.8,A,F12.8)')
C    $           'OUTNUC=(',PNUC(1),  ',',PNUC(2),  ',',PNUC(3),
C    $           ')/Eoutnuc=',enf


C            write(*,*) '---------------------------------------'

            RETURN

         ENDIF
 100  if (QUIET.eq.0) write(*,*) 'neelslpv: retry..'
      CONTINUE
      WRITE(6,*) ' IN NEELSLPV ( TOO MANY TRY ) '
C     -Add 97/08/19 _ Y.H. for errortrap
      IERR=1
      RETURN
C     
C     ++ ERROR RETURN
C     
 9000 WRITE(6,900)
 900  FORMAT(' ***ERROR IN NEELSLPV(NOT ENOUGH NEUTRINO ENERGY)***')
      IERR=1 
      RETURN
      END
