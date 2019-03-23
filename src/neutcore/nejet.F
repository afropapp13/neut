***********************************************************************
      SUBROUTINE NEJET(W, MODE, IERR)
*                                                                     *
*     Multi hadron production using PYTHIA                            *
*                                                                     *
*     Input                                                           *
*       COMMON blocks                                                 *
*         NEWORK                                                      *
*       Arguments                                                     *
*         W: Invariant mass of hadron system [GeV]                    *
*         MODE: Interaction mode                                      *
*                 +-21: CC-Npi, +-41: NC-Npi                          *
*                                                                     *
*     Output                                                          *
*       COMMON blocks                                                 *
*         NEWORK, LUJETS                                              *
*       Arguments                                                     *
*         IERR: Status                                                *
*                 1: OK, 0: NG                                        *
*                                                                     *
*     (creation date and author)                                      *
*  1998.02.?? ; 'tatakidai' by M.Eto                                  *
*  1998.03.13 ;  K.Kobayashi      some bugs fix                       *
*  2008.02.18 ;  G.Mitsuka        Give up if #of try is over 10000
***********************************************************************

      IMPLICIT NONE


#include "nework.h"
#include "lujets.h"
#include "pysubs.h"
#include "pypars.h"
#include "ludat3.h"
#include "necard.h"

#include "vcwork.h"
#include "posinnuc.h"
      INTEGER MODE, IERR
      REAL W 
      INTEGER LUCOMP
      INTEGER INEU, IHAD
C      INTEGER FF,INO16
      INTEGER FF
      REAL COSTH, PFABS, IWAVE, DIRN(3)
      REAL PNEU(3), PHAD(3)
      INTEGER I, J, ITMP1, ITMP2, NCCC, ITYP, IMODE
      REAL DUMMY, PTMP1(3), PTMP2(3), WPY, TOTP(3), ENENE(MAXNE), TOTE
      CHARACTER*8 CHNEU, CHHAD

C      COMMON /SAVXIN/INO16

*     *** FOR PDFLIB(Saji Aug. 02)
      INTEGER  IDBMUP,PDFGUP,PDFSUP
      REAL  EBMUP
      COMMON/HEPRUP/IDBMUP(2), EBMUP(2), PDFGUP(2), PDFSUP(2)
      Integer ntry,ntry1

      ntry = 0
      ntry1 = 0
*-----------------------------------------------------------------*
*     NEWORK -> Auto variable, and initialize some parameters     *
*-----------------------------------------------------------------*

      INEU = IPNE(1)
      IHAD = IPNE(2)

 300  CONTINUE
      DO I=1,3
         PNEU(I) = PNE(I,1)
         PHAD(I) = PNE(I,2)
      ENDDO

      IF (INEU.EQ.12) THEN
         CHNEU = 'nu_e  '
      ELSE IF (INEU.EQ.-12) THEN
         CHNEU = 'nu_e~ '
      ELSE IF (INEU.EQ.14) THEN
         CHNEU = 'nu_mu '
      ELSE IF (INEU.EQ.-14) THEN
         CHNEU = 'nu_mu~'
      ELSE IF (INEU.EQ.16) THEN
         CHNEU = 'nu_tau '
      ELSE IF (INEU.EQ.-16) THEN
         CHNEU = 'nu_tau~ '
      ELSE
         WRITE(0,*) 'nejet: Error: Invalid neutrino type, INEU=', INEU
         IERR = 0
         RETURN
      ENDIF

      IF (IHAD.EQ.2112) THEN
         CHHAD = 'n0    '
      ELSE IF (IHAD.EQ.2212) THEN
         CHHAD = 'p+    '
      ELSE
         WRITE(0,*) 'nejet: Error: Invalid nucleon type, IHAD=', IHAD
         IERR = 0
         RETURN
      ENDIF


*---------------------------*
*     Initialize PYTHIA     *
*---------------------------*

*     Lower edge of allowed sqrt{s} [GeV]   --->  pysubs.h
      CKIN(1)     = 0.001
*     Lower cut-off on p_t [GeV/c]          --->  pysubs.h
      CKIN(5)     = 0.0001
*     Lower CM energy [GeV]                 --->  pypars.h
      PARP(2)     = 0.001
*     Switch to be allowed to decay or not  --->  ludat3.h
      MDCY(LUCOMP(111),1) = 0
      MDCY(LUCOMP(221),1) = 0
      MDCY(LUCOMP(311),1) = 0
      MDCY(LUCOMP(223),1) = 0
      MDCY(LUCOMP(130),1) = 0
      MDCY(LUCOMP(310),1) = 0
      MDCY(LUCOMP(331),1) = 1
**** without tau decay(decay at tauola)
      IF(ITAUFLGCORE.eq.1) THEN
         MDCY(LUCOMP(15),1) = 0
      ENDIF

*-----> USE PDFLIB(1: INTERNAL PYTHIA, 2:USE PDFLIB)
         MSTP(52) = 2
*-----> author gourp(PDFGUP(1) = 5: GRV)
         PDFGUP(1) = 5
*-----> PDF set(PDFSUP(1) = 7: GRV94DI)
C         PDFSUP(1) = 7 ! PDF is set by neut.card
         PDFSUP(1) = NEPDF
*-----> USE GRV94(1000: nucleon)
         MSTP(51) = 1000*PDFGUP(1)+PDFSUP(1)
         
C#ifndef NO_GRV94CORR
*     W^2 range allowed in DIS process(approximately cut)[GeV^2]
*     [D=4.]
C      CKIN(39)    = 0.
C#else
C      CKIN(39)    = 0.4
C#endif
      if(NEBODEK.eq.1) then
         CKIN(39)    = 0.
      elseif(NEBODEK.eq.0) then
         if(NEPDF.eq.7) then
            CKIN(39)    = 0.23
         elseif(NEPDF.eq.12) then
            CKIN(39)    = 0.8
         endif
      endif

      DO I=1,3
         P(1,I) = PNEU(I)
         P(2,I) = PHAD(I)
      ENDDO

      CALL PYINIT('USER', CHNEU, CHHAD, DUMMY)
      
*------------------------------------------*
*     Generate events till W > 2.0 GeV     *
*------------------------------------------*

 100  CONTINUE
      CALL PYEVNT
*      CALL LULIST(1)
      NCCC = -1
      DO I=3,N
*     Check type of lepton in final state      
         IF (NCCC.EQ.-1 .AND. K(I,1).EQ.1) THEN
            ITYP = IABS(K(I,2))
            IF (ITYP.EQ.11 .OR. ITYP.EQ.13 .OR. ITYP.EQ.15) THEN ! CC
               NCCC = 1
               FF = I
            ELSE IF (ITYP.EQ.12 .OR. ITYP.EQ.14 .OR. ITYP.EQ.16) THEN ! NC
               NCCC = 0
               FF = I
            ENDIF
         ENDIF
      ENDDO
*           Check COSTH
      if ( (SQRT(P(1,4)**2-P(1,5)**2)*SQRT(P(FF,4)**2-P(FF,5)**2))
     $    .ge.1.e-7) then

         COSTH = (P(1,1)*P(FF,1)+P(1,2)*P(FF,2)+P(1,3)*P(FF,3))/
     &        (SQRT(P(1,4)**2-P(1,5)**2)*SQRT(P(FF,4)**2-P(FF,5)**2))

         IF(ntry.gt.1000) goto 400
         IF(COSTH.LT.-0.75) THEN
            ntry1=ntry1+1
            If (mod(ntry1,100).eq.0) then
               print *, 'nejet: retry ',ntry1,' times to get COSTH>-0.75 '
               print *, '     : E,mode=',ineu,mode
               print*,'Pnu=',pne(1,1),pne(2,1),pne(3,1)
               print*,'Pp =',pne(1,2),pne(2,2),pne(3,2)
C               IF (INO16.EQ.1) THEN
               IF (IBOUND.EQ.1) THEN
C                  CALL RNFERM(PFABS,IWAVE)
                  CALL nefrmmom(PFABS,IWAVE)
                  PFABS=PFABS*1.0E-3
               ENDIF
               CALL RNDIR(DIRN)
               PNE(1,2)=DIRN(1)*PFABS
               PNE(2,2)=DIRN(2)*PFABS
               PNE(3,2)=DIRN(3)*PFABS
               GOTO 300
            endif
            GOTO 100
         ENDIF
      endif

*        Check interaction mode     
 400  IMODE = IABS(MODE)
      IF (((IMODE.EQ.26) .AND. NCCC.EQ.0)
     +     .OR. ((IMODE.EQ.46) .AND. NCCC.EQ.1)) THEN
         GOTO 100
      ENDIF


*--------------------------*
*     LUJETS -> NEWORK     *
*--------------------------*

      NUMNE = 2     ! There have been 2 particles in initial state
      DO I=1,N 
         IF (K(I,1).EQ.1) THEN          ! If final state
            NUMNE = NUMNE + 1
            IPNE(NUMNE) = K(I,2)
            DO J=1,3
               PNE(J,NUMNE) = P(I,J)
            ENDDO
            ENENE(NUMNE) = P(I,4)
         ENDIF
      ENDDO


*---------------------------------------------------------------*
*     Move nucleon in final state into IPNE(4) and PNE(I,4)     *
*---------------------------------------------------------------*

      IF (IPNE(4).LT.1000) THEN          ! If not baryon 
         DO I=5,NUMNE
            IF (IPNE(I).GT.1000) THEN
               ITMP1 = IPNE(4)
               ITMP2 = IPNE(I)
               DO J=1,3
                  PTMP1(J) = PNE(J,4)
                  PTMP2(J) = PNE(J,I)
               ENDDO
               IPNE(4) = ITMP2
               IPNE(I) = ITMP1
               DO J=1,3
                  PNE(J,4) = PTMP2(J)
                  PNE(J,I) = PTMP1(J)
               ENDDO
               GOTO 200
            ENDIF
         ENDDO
      ENDIF


 200  CONTINUE

      DO I=4,NUMNE
         ICRNNE(I) = 1 
         IF((IPNE(I).EQ.12).OR.(IPNE(I).EQ.14).OR.(IPNE(I).EQ.16)) THEN
            ICRNNE(I) = 0
         ENDIF
      ENDDO
      IORGNE(3) = 1
      DO I=4,NUMNE
         IORGNE(I) = 2
      ENDDO

*                   CHECK W>2.0GEV

      DO J=1,3
         TOTP(J) = 0.
      ENDDO
      TOTE = 0.

      DO I=4,NUMNE
         DO J=1,3
            TOTP(J) = TOTP(J) + PNE(J,I)
         ENDDO
         TOTE = TOTE + ENENE(I)
      ENDDO
      WPY = SQRT(TOTE**2-TOTP(1)**2-TOTP(2)**2-TOTP(3)**2)
      IF (WPY.LT.2.0) THEN
         ntry=ntry+1
         If (mod(ntry,100).eq.0) then
            print *, 'nejet: retry ',ntry,' times to get W>2.0GeV.'
            print *, '     : E,mode=',ineu,mode
            print*,'Pnu=',pne(1,1),pne(2,1),pne(3,1)
            print*,'Pp =',pne(1,2),pne(2,2),pne(3,2)
C            IF (INO16.EQ.1) THEN
            IF (IBOUND.EQ.1) THEN
C               CALL RNFERM(PFABS,IWAVE)
               CALL nefrmmom(PFABS,IWAVE)
               PFABS=PFABS*1.0E-3
            ENDIF
            CALL RNDIR(DIRN)
            PNE(1,2)=DIRN(1)*PFABS
            PNE(2,2)=DIRN(2)*PFABS
            PNE(3,2)=DIRN(3)*PFABS
            IF (ntry.eq.1000) THEN
		IF(WPY.GT.(2.0-0.05*Real(ntry/1000))) THEN
                   GOTO 600		
		ENDIF	
            ENDIF
            IF (ntry.eq.10000) THEN ! Give up this mode.
               IERR = 2
               GOTO 601
            ENDIF
            GOTO 300
         Endif
         GOTO 100
      ENDIF

  600 CONTINUE

      IERR = 1

  601 CONTINUE

      RETURN
      END




