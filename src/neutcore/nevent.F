**********************************************************************
*     ----------------------------
      SUBROUTINE NEVENT(IPAR,PMOM,IRET)
*     ----------------------------
*     
*     ( purpose )
*       VECTOR GENERATION ( NEUTRINO INTERACTION )
*
*     ( input ) 
*       IPAR    : PARICLE CODE
*       PMOM(3) : MOMENTUM 
*
*     ( output ) 
*       COMMON NEWORK
*       IRET    : ERROR status
C                  No errors       : 0      
C                  Error condition : 1
*
*     ( creation date and author )
*       1985.01.12 ; First version by M.Nakahata
*       1987.07.29 ; N.Sato  include tau
*       1988.10.08 ; T.Kajita ( N.C. ELASTIC INT. )
*       1995.02.25 ; K.Kaneyuki for S.K.
*                    NEVENT -> NEVENT      
*       1997.08.19 ; Y.Hayato  / debug 
*                    (RE-try if kinematical error happened)
*       1997.12.01 ; J.Kameda (single eta production)
*       1997.12.01 ; J.Kameda (single  K  production)
*
*       2006.06.12 ; Y.Hayato  Use IBOUND ; Remove INO16 & IFREE
*                    Add IMPLICIT NONE
*
*       2006.12.31 ; Y.Hayato  Add IRET to indicate errors
C                              and not to change energy of neutrino      
*       2007.08.23 ; G.Mitsuka add delta->gamma decay
*       2013.02.28 ; P.Sinclair Add multi-nucleon mode
*
**********************************************************************
      IMPLICIT NONE

#include "neutparams.h"
#include "nework.h"
#include "vcwork.h"
#include "necard.h"
#include "neutmodel.h"
      
      integer*4 IPAR
      REAL PMOM(3)

      integer*4 nemodsel
      external  nemodsel

      REAL DIRNEU(3), PMOMG(3), PMOMM(3), xpos(3), tmppos(3)
C     
#include "posinnuc.h"
C---  IBOUND = 1 : Bound nucleon
C---  IBOUNT = 0 : Free proton
C     
      REAL*4 APN
      DATA APN/16./

      real*4 E,eneu

      integer*4 mode,modabs
      integer*4 ierr,i,j

      integer*4 iret

      integer*4 iskip

      integer*4 fshadpos
c -- dummy variable used in spectral function elastic events - set to true for event generation
      logical first
      first = .TRUE.
C
      iret = 1
C     
C --  SET NEUTRINO ENERGY
C
 1    continue
      PMOMG(1)=PMOM(1)*1.0E-3
      PMOMG(2)=PMOM(2)*1.0E-3
      PMOMG(3)=PMOM(3)*1.0E-3
      E = SQRT(PMOMG(1)**2+PMOMG(2)**2+PMOMG(3)**2)
      IF (E.GT.0.001) THEN
         DIRNEU(1)=PMOMG(1)/E
         DIRNEU(2)=PMOMG(2)/E
         DIRNEU(3)=PMOMG(3)/E
      ELSE
         DIRNEU(1)=1.0
         DIRNEU(2)=0.0
         DIRNEU(3)=0.0
      ENDIF
      
      if (QUIET.le.1) then
         write(*,*) 'Enu: ',E
      endif

C      DO 10 I=1,1000
C     
C --  SELECT INTERACTION MODE
C     
         CALL NECLRCRS

         MODE=NEMODSEL(IPAR,E)
         MODENE=MODE
         MODABS=IABS(MODE)

C     
C --     GENERATE VECTOR
C     
         CALL NECLRWRK
         IF (MODE.EQ.0) RETURN
         IF (MODABS.EQ.1) THEN
C
C --     ELASTIC PROCESS
C            
C --     Elastic process *with free proton*
C            IF (IFREE.NE.0) THEN
            IF (IBOUND.EQ.0) THEN
C-DEBUG/Y.H/96-7-26               CALL NEELSVCT(IPAR,E,MODE,DIRNEU)
               CALL NEELSVCT(IPAR,E,DIRNEU)
               GOTO 1000
            ELSE
C
C --      ELASTIC PROCESS IN O16
C

C --      Use spectral function
           IF (mod(MDLQE, 10000)/1000 .EQ. 2) then
              call ne_1p1h_n(ipar,mode,e,dirneu,ierr)
           else IF (mod(MDLQE, 1000) / 100 .EQ. 4) then
              CALL SFEVENT(IPAR,E,DIRNEU,IERR,.TRUE.,mode,first)
           ELSE IF ( (mod(MDLQE, 1000) / 100 .EQ. 6).or.
     $               (mod(MDLQE, 1000) / 100 .EQ. 7) ) then
              CALL EFFSFEVENT(IPAR,E,DIRNEU,IERR,.TRUE.,mode,first)
           ELSE
              CALL NEELSVCP(IPAR,E,DIRNEU,IERR)
           ENDIF

C
C --           SAVE INTERACTION POSITION IN NUCEI
C               CALL EFRANBLL(APN,TMPPOS)
C               INO16=1
C               IBOUND=1

               IF(IERR.EQ.0)GO TO 1000
               WRITE(6,620)MODE,E
 620           FORMAT(' *** WARNING *** KINEMATICAL ERROR',
     &              ' AT ELSVCP (MODE,E) ',I3,G15.7)
               GO TO 10
            ENDIF
         ELSE IF (MODABS.EQ.2) THEN

C --     Error if proton is free
            IF (IBOUND.EQ.0) THEN
               WRITE(6,640)MODE,E
 640           FORMAT(' *** WARNING *** CANT HAVE MEC ON FREE PROTONS',
     &              I3,G15.7)
               GO TO 10
            ELSE
C
C --      MEC PROCESS
C               
               if (MDL2P2H.eq.2) then 
                  call ne_2p2h_ht(ipar,0,E,dirneu,ierr) 
               else if (MDL2P2H.eq.1) THEN
                  CALL NEMECVCP(IPAR,E,DIRNEU,IERR)
               else
                  write(*,*) 'Unknown MDL2P2H value ',MDL2P2H
                  stop
               endif
               
               IF(IERR.EQ.0)GO TO 1000
               WRITE(6,630)MODE,E
 630           FORMAT(' *** WARNING *** KINEMATICAL ERROR',
     &              ' AT MECVCP (MODE,E) ',I3,G15.7)
               GO TO 10
            ENDIF


         ELSE IF (( 10.LT.MODABS .AND. MODABS.LT.20 ).OR.
     &            ( 30.LT.MODABS .AND. MODABS.LT.40 )) THEN
C     
C --        SINGLE PION PROCESS
C
            IF ( MODABS.EQ.16 .OR. MODABS.EQ.36) THEN
C
C  -- COHERENT PI PRODUCTION - T.KAJITA AND K.OKUMURA 06-SEP-94---
C
               CALL NECOHVCT(IPAR,E,MODE,DIRNEU,IERR)
               IF (IERR.ne.0) then
                  WRITE(6,650)MODE,E
 650              FORMAT(' *** WARNING *** KINEMATICAL ERROR',
     &                 ' AT NECOHVCT (MODE,E) ',I3,G15.7)
                  goto 10
               endif
               GOTO 1000
C     
C     --           DIFFRACTIVE PION PROCESS(resonance -> nucleon + gamma)
C     
            ELSEIF ( MODABS.EQ.15 .OR. MODABS.EQ.35 ) THEN
               IF (IBOUND.EQ.0) THEN
                  CALL NEDIFVCT(IPAR,E,MODE,DIRNEU)
                  GOTO 1000
               ELSE
                  GOTO 10
               ENDIF
C     
C --           SINGLE PION PROCESS(resonance -> nucleon + gamma)
C     
            ELSEIF ( MODABS.EQ.17 .OR. MODABS.EQ.38 .OR. MODABS.EQ.39 ) THEN
               CALL RSSGVCT(IPAR,MODE,E,DIRNEU,IERR)
               IF(IERR.EQ.0)GO TO 1000
               WRITE(6,601)MODE,E
 601           FORMAT(' *** WARNING *** KINEMATICAL ERROR',
     &              ' AT SGVCT (MODE,E) ',I3,G15.7)
               GO TO 10
C     
C --           SINGLE PION PROCESS
C
            ELSE
               CALL RSSPIVCT(IPAR,MODE,E,DIRNEU,IERR)
               IF(IERR.EQ.0)GO TO 1000
               WRITE(6,600)MODE,E
 600           FORMAT(' *** WARNING *** KINEMATICAL ERROR',
     &              ' AT SPIVCT (MODE,E) ',I3,G15.7)
               GO TO 10
            ENDIF
C            
C --        MULTI PION PROCESS
C
C         ELSE IF( MODABS .EQ. 21 .OR. MODABS .EQ. 41) THEN
         ELSE IF( MODABS .EQ. 21 .OR. MODABS .EQ. 41
     &           .or. MODABS .EQ. 26 .OR. MODABS .EQ. 46) THEN
            CALL NEAPIVCT(IPAR,MODE,E,DIRNEU,IERR)
            IF(IERR.EQ.0)GO TO 1000
            WRITE(6,610)MODE,E
 610        FORMAT(' *** WARNING *** KINEMATICAL ERROR',
     &           'AT APIVCT (MODE,E) ',I3,G15.7)
            GO TO 10
C     
C --        N.C. ELASTIC PROCESS         08-OCT-88   T.KAJITA
C            
         ELSE IF( MODABS .EQ. 51 .OR. MODABS.EQ.52 ) THEN
            IF (IBOUND.NE.1) THEN
               CALL NEELSVNT(IPAR,E,MODE,DIRNEU)
               GO TO 1000
            ELSE
C
C --           N.C. ELASTIC PROCESS IN O16

               if (mod(MDLQE, 1000) / 100 .EQ. 4) then
c               use spectral function               
                  CALL SFEVENT(IPAR,E,DIRNEU,IERR,.TRUE.,mode,first)
               else if (mod(MDLQE, 1000) / 100 .EQ. 6) then
                  CALL EFFSFEVENT(IPAR,E,DIRNEU,IERR,.TRUE.,mode,first)
               else
                  CALL NEELSVNP(IPAR,E,DIRNEU,IERR,MODE)
               endif
C
C --           SAVE INTERACTION POSITION IN NUCEI
C               
C               CALL EFRANBLL(APN,TMPPOS)
C               INO16=1
C               IBOUND=1

               IF(IERR.EQ.0)GO TO 1000
               WRITE(6,621)MODE,E
 621           FORMAT(' *** WARNING *** KINEMATICAL ERROR ',
     &              'AT ELSVNP (MODE,E) ',I3,G15.7)
               GO TO 10
            ENDIF
         ELSE IF( MODABS .EQ. 22 .OR. MODABS.EQ.42 .OR. 
     &       MODABS .EQ. 43) THEN
CC ----------      SINGLE ETA PRODUCTION --1997.12.01.J.KAMEDA

               CALL RSSETVCT(IPAR,MODE,E,DIRNEU,IERR)
               IF(IERR.EQ.0)GO TO 1000
               WRITE(6,603)MODE,E
 603           FORMAT(' *** WARNING *** KINEMATICAL ERROR',
     &              'AT RSSETVCT (MODE,E) ',I3,G15.7)
               GO TO 10

         ELSE IF( MODABS .EQ. 23 .OR. MODABS.EQ.44 .OR. 
     &       MODABS .EQ. 45) THEN
CC ----------      SINGLE  K  PRODUCTION --1998.02.25.J.KAMEDA

               CALL RSSKVCT(IPAR,MODE,E,DIRNEU,IERR)
               IF(IERR.EQ.0)GO TO 1000
               WRITE(6,604)MODE,E
 604           FORMAT(' *** WARNING *** KINEMATICAL ERROR',
     &              'AT RSSKVCT (MODE,E) ',I3,G15.7)
               GO TO 10
             
         ENDIF
 10   CONTINUE

      WRITE(6,*) '*** WARNING *** TOO MANY ERR AT NEVENT'
      WRITE(6,*) 'Need to change the energy and direction of neutrino'

      IRET = 1
      RETURN

 1000 CONTINUE

C-- set interaction position here (2006/06/12)

      if (ibound.eq.-1) then
         write(*,*) 'nevent: IBOUND flag is still -1..'
         write(*,*) 'nevent: mode=',modene
         stop
      endif         

      if (IBOUND.eq.0) then
         TMPPOS(1)=9999.
         TMPPOS(2)=9999.
         TMPPOS(3)=9999.
         do 1150 J=1,NUMNE
            do 1050 I=1,3
               POSNUC(I,J)=TMPPOS(I)
 1050       continue
 1150    continue
      endif

      if (IBOUND.eq.1) then

C-- set initial interactions for leptons and nucleon
C                    1: neutrino
C                    2: nucleon
C                    3: lepton
C
C Position of nucleons different for MEC events. Set fshadpos to loop from correct position
C Position of lepton and initial state nucleon is set in nemechad for MEC events, skip to 900 

         IF(MODABS.eq.2) THEN
           fshadpos=5
           goto 900
         ELSE 
           fshadpos=4
         ENDIF

CC-- if position has been set already,
         if ( ( posnuc(1,1) .ge. -9998 ) .or.
     $        ( posnuc(2,1) .ge. -9998 ) .or.
     $        ( posnuc(3,1) .ge. -9998 )) THEN
            xpos(1)=posnuc(1,1)
            xpos(2)=posnuc(2,1)
            xpos(3)=posnuc(3,1)
            goto 900
         endif

         DO 1010 J=1,3
            PMOMM(J) = 0.
 1010    continue

         CALL NERANBLL(APN,PMOMM,TMPPOS,XPOS)
         do 1030 I=1,3
            do 1020 J=1,3
               POSNUC(J,I)=XPOS(J)
 1020       continue
C            write(*,'(I3,A3,I5,A5,F8.3,A1,F8.3,A1,F8.3,A1)')
C     $           I,'IP=',ipne(i),'POS=(',posnuc(1,I),',',
C     $                                   posnuc(2,I),',',
C     $                                   posnuc(3,I),')'
 1030    continue
         
C-- set initial interactions position for the other particles
C                    4: nucleon
C                    5: pion(?)
C

 900     do 1040 I=fshadpos,numne
            DO 1070 J=1,3
               PMOMM(J) = PNE(J,I)
 1070       continue
            ISKIP = 0

            if  ((IFORMLEN.eq.110).and.
     $           ( (abs(mode).eq.1).or.
     $             (abs(mode).eq.2).or.
     $             (abs(mode).eq.51).or.
     $             (abs(mode).eq.52))) then
               ISKIP = 1
            endif

            if (IFORMLEN.eq.100) THEN
               if ( (abs(mode).eq.21).or.
     $              (abs(mode).eq.26).or.
     $              (abs(mode).eq.41).or.
     $              (abs(mode).eq.46)) then
                  ISKIP = 0
               ELSE
                  ISKIP = 1
               endif
            endif
            
            if (ISKIP.eq.1) then
               DO 1110 J=1,3
                  POSNUC(J,I) = XPOS(J)
 1110          continue
            else
               CALL NERANBLL(APN,PMOMM,XPOS,TMPPOS)
               do 1060 J=1,3
                  POSNUC(J,I) = TMPPOS(J)
 1060          continue
            endif
            
C           write(*,'(I3,A3,I5,A10,F8.3,A1,F8.3,A1,F8.3,A1)')
C     $           I,'IP=',ipne(i),'POS_NEV=(',posnuc(1,I),',',
C     $                                   posnuc(2,I),',',
C     $                                   posnuc(3,I),')'
 1040    continue
         
      endif

      if (QUIET.le.1) then
         write(*,*) 'Mode: ',MODENE
      endif

C-- no errors
      IRET = 0

      RETURN
      END
