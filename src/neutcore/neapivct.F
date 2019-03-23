**********************************************************************
*     ------------------------------------------------
      SUBROUTINE NEAPIVCT(IPAR,JMOD0,ENEUT,DIRNEU,IERR)
*     ------------------------------------------------
*     
*     ( purpose )
*       VECTOR GENERATION FOR MULTI PION PRODUCTION
*
*     ( input ) 
*       IPAR    : PARICLE CODE
*       JMOD0   : INTERACTION MODE
*       ENEUT   : ENERGY OF NEUTRINO ( GEV )
*       DIRNEU(3) : NEUTRINO DIRECTION
*
*     ( output ) 
*       IERR    : ERROR CODE
*       COMMON NEWORK
*
*     ( creation date and author )
*       1983.11.27 ; First version by M.Nakahata
*       1985.11.06 ; K.KAJITA SMALL  MODIFICATION
*                       ( FOR MULTI PION    DELTA PI  --> N PI PI)
*                       ( FORWARD BACKWARD PION MULTIPLICITY )
*       1988.09.08 ; T.KAJITA CROSS SECTION AT W=1.3-1.4 WITH N(PAI)=2
*                       IS ADDED.
*       1988.10.08 ; T.KAJITA  SMALL MOD. FLAG INO16 IS SET IN THIS PROGRAM
*       1988.10.18 ; T.KAJITA  FERMI MOTION OF NUCLEON IS INCLUDED
*       1995.02.25 ; K.Kaneyuki for S.K.
*                    APIVCT -> NEAPIVCT
*	1995.10.29 ; K.Knaeyuki BUG FIX 
*       1997.12.12 ; Y.Hayato   To use IFREE to determine FREE or BOUND.
*       1999.08.03 ; M.Shiozawa modify procedure for re-calculate 
*                    the pion multiplicity
*       2006.06.12 ; Y.Hayato   Not to use IFREE nor INO16
*                               use IBOUND                     
*       2007.11.05 ; G.Mitsuka  p/n ratio is calculated 
*                               using target information
*       2007.11.10 ; T.Tanaka   add upmu mode
*       2008.02.18 ; G.Mitsuka  Bug fix
*       2016.03.08 ; C.Bronner  W (x,y) generation now depends on target nucleon
*                               Target nucleon chosen based on comparison of xsec
*
*
**********************************************************************
      IMPLICIT NONE

#include "vcwork.h"
#include "nework.h"
#include "necard.h"
C
C     PNRAT -- P : N TSRGET RATIO
C
      integer*4 ipar,jmod0
      real*4    eneut
      REAL*4    DIRNEU(3)
      integer*4 ierr
      REAL*4 FNAPICRS
      EXTERNAL FNAPICRS
      real*4    rlu,rnmass,fngamfct
      external  rlu,rnmass,fngamfct

      integer*4 NEMULTIP
      external  NEMULTIP

      REAL*4    DIRN(3), PNEU(3)
      REAL*4    DIR(3),  GDIR(3),PHAD(3),PPHA(3,10),PBUF(3)
      REAL*4    EV(3),EVN(3)
      REAL*4    DIRTMP(3), DIRRS(3)
C      COMMON/MASS/AM(9)
      REAL*4    AM(9)

      real*4    E,W,WMEV,X,Y,q2


#include "posinnuc.h"
C      COMMON/FREPRO/IFREE
C      COMMON /SAVXIN/INO16,XIN(3)

      real*4 PMASS
      DATA PMASS/0.938/
      
      real*4 PNRAT(3)
C      DATA PNRAT/0.3846 , 0.7143 , 0.5/
      real*4 PXSEC, NXSEC, R1320P, R1320N
      integer*4 MAXERR
      DATA MAXERR/10/

      real*4 FRFREE
c     FRFREE is derived by NUMBNDP and NUMFREP
c      DATA FRFREE/0.2/

      INTEGER IERJET
      Integer idummy

      integer*4 MULP,MULPD,MUL,NFPI,NBPI
      integer*4 IEC,IWAVE,itype,iwmode,IPN,IPP,ipbuff
      real*4    ratio,ran,dum,prin

      integer*4 I,J,K,ierrr,ierrp,ichk

      real*4    pfabs,pabs,enucl,enstop,beta,betars,gm,gam,gmh
      real*4    elep,amlep,plep,cost,pphad,totmas,ptotrs
      real*4    bm

      real*4 NUCMASS,W2Max

C                FRFREE --- FREE PROTON/TOTAL PROTON     08-OCT-88
C
      E=ENEUT
      IERR=0
      IEC=0

C     Determine p/n fraction with cross section ratio is taken into account
C     Note:9/13 and 9/7 is given by cross section ratio(p/total)
C     (2007 Nov G.Mitsuka)
C      PNRAT(1)=9./13.*float(NUMBNDP+NUMFREP)/float(NUMBNDN+NUMBNDP+NUMFREP)
C      PNRAT(2)= 9./7.*float(NUMBNDP+NUMFREP)/float(NUMBNDN+NUMBNDP+NUMFREP)

C CB Use ratio of cross sections instead of number of valence quarks
C CB     PNRAT(1)=float(NUMBNDP+NUMFREP)
C CB     $         /float(2*NUMBNDN+NUMBNDP+NUMFREP)
C CB      PNRAT(2)=float(2*(NUMBNDP+NUMFREP))
C CB     $     /float(NUMBNDN+2*(NUMBNDP+NUMFREP))

C      Previously assumed iso-scalar nucleus
C      PNRAT(3)=0.5
C CB      PNRAT(3) = float(NUMBNDP+NUMFREP)/float(NUMBNDP+NUMFREP+NUMBNDN)
      
C CB Determine whether interaction is CC or NC      
      ITYPE=1
      IF(IABS(JMOD0).EQ.41.or.IABS(JMOD0).EQ.46)ITYPE=0

C CB Compute cross section on neutron and on proton
      PXSEC=FNAPICRS(E,IPAR,ITYPE,1.,R1320P)
      NXSEC=FNAPICRS(E,IPAR,ITYPE,0,R1320N)
      If (abs(JMOD0).eq.21.or.abs(JMOD0).eq.41) then
         PXSEC=PXSEC*R1320P
         NXSEC=NXSEC*R1320N
      endif
C CB      PRINT *,'PXSEC, NXSEC: ',PXSEC, NXSEC
C     Free hydrogen target should always be a proton
C CB      if (NUMBNDP.eq.0 .and. NUMBNDN.eq.0) PNRAT(3) = 1.0
      if (NUMBNDP.eq.0 .and. NUMBNDN.eq.0) then
         RATIO = 1.0
      else
         RATIO=float(NUMBNDP+NUMFREP)*PXSEC/(float(NUMBNDN)*NXSEC+float(NUMBNDP+NUMFREP)*PXSEC)
      endif
      
C CB      PRINT *,'Enu, ratio: ',E,RATIO
      
      FRFREE  =float(NUMFREP)/float(NUMBNDP+NUMFREP)
#ifdef UPMU_RCK
C CB      PNRAT(1)=0.5
C CB      PNRAT(2)=0.5
C CB      PNRAT(3)=0.5
      RATIO=0.5
      FRFREE  =0.0
#endif

C
C -- SET TARGET PARTICLE AND INTERACTION MODE
C
      MODENE=JMOD0
      IPNE(1)=IPAR
      IFLGNE(1)=-1
C CB      IF(JMOD0.EQ.21)RATIO=PNRAT(1)
C CB      IF(JMOD0.EQ.26)RATIO=PNRAT(1)
C CB      IF(JMOD0.EQ.41)RATIO=PNRAT(3)
C CB      IF(JMOD0.EQ.46)RATIO=PNRAT(3)
C CB      IF(JMOD0.EQ.-21)RATIO=PNRAT(2)
C CB      IF(JMOD0.EQ.-26)RATIO=PNRAT(2)
C CB      IF(JMOD0.EQ.-41)RATIO=PNRAT(3)
C CB      IF(JMOD0.EQ.-46)RATIO=PNRAT(3)
      RAN=RLU(DUM)
      IPNE(2)=2112
      IF(RAN.LT.RATIO)IPNE(2)=2212
      IFLGNE(2)=-1
C     
C -- SET INO16 FLAG     08-OCT-88 T.KAJITA----------------------------
C      INO16=1
      IBOUND = 1
      RAN=RLU(DUM)
C      IF(IPNE(2).EQ.2212 .AND. RAN.LT.FRFREE) INO16=0
      IF(IPNE(2).EQ.2212 .AND. RAN.LT.FRFREE) IBOUND=0
C      IF (INO16.eq.0) IFREE=1
C -------------------------------------------------------------------
C
      IF(JMOD0.EQ.21.or.JMOD0.EQ.26) THEN
         IPNE(3)=IPAR-1
         IORGNE(3)=1
         ICRNNE(3)=1
      ELSE IF (JMOD0.EQ.-21.or.JMOD0.EQ.-26) THEN
         IPNE(3)=IPAR+1
         IORGNE(3)=1
         ICRNNE(3)=1
      ELSE IF(IABS(JMOD0).EQ.41.or.IABS(JMOD0).EQ.46) THEN
         IPNE(3)=IPAR
         IORGNE(3)=1
C -- BUG FIX by K.K. 29-Oct-1995
C	 IFLGNE(3)=1
	 IFLGNE(3)=2
      ENDIF
C
C -- SET NEUTRINO DIRECTION AND TARGET MOMENTUM
C
C      CALL RNNEUT(E,DIR,IPAR)
      DO 100 I=1,3
         PNE(I,1)=E*DIRNEU(I)
         PNE(I,2)=0.
         PNEU(I)=PNE(I,1)
  100 CONTINUE
C
C
C -- CONSIDER FERMI MOTION OF NUCLEON T.KAJITA  18-NOV-88 ------------
 96   PFABS=0.
C      IF (INO16.EQ.1 .AND. ABS(IPAR).NE.16) THEN
      IF (IBOUND.EQ.1 .AND. ABS(IPAR).NE.16) THEN
C         CALL RNFERM(PFABS,IWAVE)
         CALL nefrmmom(PFABS,IWAVE)
         PFABS=PFABS*1.0E-3
      ENDIF
      CALL RNDIR(DIRN)
      PNE(1,2)=DIRN(1)*PFABS
      PNE(2,2)=DIRN(2)*PFABS
      PNE(3,2)=DIRN(3)*PFABS
C
      ENUCL=SQRT(PFABS**2+PMASS**2)
      BETA=PFABS/ENUCL
      GM=1./SQRT(1.-BETA**2)
      IF(PFABS.EQ.0.)GO TO 10
      EV(1)=-DIRN(1)
      EV(2)=-DIRN(2)
      EV(3)=-DIRN(3)
      GO TO 20
   10 EV(1)=1.
      EV(2)=0.
      EV(3)=0.
   20 CALL MCVECBST(PNEU,0.,EV,GM)
      ENSTOP=SQRT(PNEU(1)**2+PNEU(2)**2+PNEU(3)**2)
      DIR(1)=PNEU(1)/ENSTOP
      DIR(2)=PNEU(2)/ENSTOP
      DIR(3)=PNEU(3)/ENSTOP
C      IF(ENSTOP.LT.0.550 .AND. IABS(IPAR).EQ.12 .AND. INO16.EQ.1)
      IF(ENSTOP.LT.0.550 .AND. IABS(IPAR).EQ.12 .AND. IBOUND.EQ.1)
     &                                                         GO TO 96
C      IF(ENSTOP.LT.0.650 .AND. IABS(IPAR).EQ.14 .AND. INO16.EQ.1)
      IF(ENSTOP.LT.0.650 .AND. IABS(IPAR).EQ.14 .AND. IBOUND.EQ.1)
     &                                                         GO TO 96
      E=ENSTOP
C---------------------------------------------------------------------
C
C
C------------------------------------------------ ADD T.KAJITA 06-SEP-88
C -- DETERMINE W<1.4 OR W>1.4 GEV   06-SEP-88  T.KAJTIA
C
C     mode was already selected in nemodsel
C                                   -- IWMODE=0 -> W<2.0GeV
C                                      IWMODE=1 -> W>2.0GeV
      If (abs(JMOD0).eq.21.or.abs(JMOD0).eq.41) then
         IWMODE=0
      Else
         IWMODE=1

         ! To avoid infinite loops
         CALL MCMASSGV(2112,NUCMASS)
         W2Max = NUCMASS**2. + 2.*NUCMASS*ENEUT ! (GeV^2)
         if (W2Max.le.2.0**2) then            
            write(*,*),'Error in neapivct, Mode=',JMOD0,
     $           ' Enu=',ENEUT,' W=',sqrt(W2Max),'< 2.0GeV'
            goto 1000
         endif
      Endif
CMS      DUMMY= FNAPICRS(E,IPAR,1, -1.0, R1314)
CMS      IWMODE=0
CMS      IF(RNDM(DUM).GT.R1314) IWMODE=1
CMSC                                   -- IWMODE=0 -> W<1.4GEV
CMSC                                      IWMODE=1 -> W>1.4GEV
C-----------------------------------------------------------------------
C
C
C -- SET X, Y AND W
C
C CB      Move definition of itype to the beginning as it is need for xsec on p and n
C CB      ITYPE=1
C CB      IF(IABS(JMOD0).EQ.41.or.IABS(JMOD0).EQ.46)ITYPE=0
C
C------------------------------------------------ MOD T.KAJITA 06-SEP-88
CC200 CALL RANW(IPAR,E,X,Y,W,ITYPE,IERRR)    OLD T.KAJITA 06-SEP-88
  200 CONTINUE
      IERRR = 0
      IF(IWMODE.EQ.1) THEN
C     Use PYTHIA if W > 2.0 GeV
C     Added by M.Eto on 14/Feb/1998
         CALL NEJET(W, JMOD0, IERJET)
         IF (IERJET.EQ.2) GOTO 1000 ! Give up this mode
         IF (IERJET.NE.0) THEN
            IERR = 0
         ENDIF
         GOTO 5050
      ELSE
         CALL RNW1320(IPAR,E,X,Y,W,ITYPE,IPNE(2),IERRR)
      END IF
C-----------------------------------------------------------------------
      IF(IERRR.EQ.1)GO TO 1000


      Q2=PMASS**2+2.*PMASS*Y*E-W**2
      ELEP=E-Y*E
      CALL MCMASS(IPNE(3),AMLEP)
      AMLEP=AMLEP*1.0E-3
      PLEP=ELEP**2-AMLEP**2
      IF(PLEP.GE.0.)GO TO 201
      WRITE(6,610)PLEP,E,X,Y,ELEP,W
  610 FORMAT(' *** ERROR IN APIVCT *** PLEP < 0.',G15.7//
     &   '  E,X,Y,ELEP,W ',5G15.7)
      IEC=IEC+1
      IF(IEC.LT.MAXERR)GO TO 200
      IERR=1
      RETURN
  201 PLEP=SQRT(PLEP)
      COST=(2.*E*ELEP-AMLEP**2-Q2)/2./E/PLEP
      IF(ABS(COST).LE.1.)GO TO 210
      WRITE(6,600)COST,E,X,Y,PLEP,W
  600 FORMAT(' *** ERROR IN APIVCT *** !COST! > 1.',G15.7//
     &   '  E,X,Y,PLEP,W ',5G15.7)
      IEC=IEC+1
      IF(IEC.LT.MAXERR)GO TO 200
      IERR=1
      RETURN
  210 CALL RNROT(COST,DIR,GDIR)
C
C -- BOOST TO THE ORIGINAL EXP. SYSTEM  T.KAJITA  18-NOV-88 ------------
      DO 2310 I=1,3
 2310 PBUF(I)=PLEP*GDIR(I)
      CALL MCVECBST(PBUF,AMLEP,DIRN,GM)
      DO 2330 I=1,3
         PNE(I,3)=PBUF(I)
         PHAD(I)=PNE(I,1)+PNE(I,2)-PNE(I,3)
 2330 CONTINUE

      PPHAD=SQRT(PHAD(1)**2+PHAD(2)**2+PHAD(3)**2)
      GMH=SQRT(W**2+PPHAD**2)/W
      DO 230 I=1,3
  230 EV(I)=PHAD(I)/PPHAD
C
C -- SET PION MULTIPLICITY
C
C     set counter = 0
      idummy=0
 240  idummy = idummy+1
      MULP=NEMULTIP(W)
      MUL=MULP+1
C === NEXT LINE WAS ADDED BY M.N ON 9-AUG-88
      IF(MULP.GE.8)MULP=7
      IF(MUL.GE.9)MUL=8
C------------------------------------------------ ADD T.KAJITA 06-SEP-88
CMS      IF(IWMODE.EQ.0 .AND. MULP.NE.2 ) GO TO 200
CC--------------------------------------------MODIFIED by M.Shiozawa 02-MAR-98
C      IF(IWMODE.EQ.0 .AND. MULP.LE.2 ) GO TO 200
C--------------------------------------------BUG FIXED by M.Shiozawa 28-APR-99
C      IF(IWMODE.EQ.0 .AND. MULP.LT.2 ) GO TO 200
CC-----------------------------------------------------------------------
C----------------------------------- W is not recalculated
C----------------------------------- (only pion multiplicity is recalculated)
      IF(IWMODE.EQ.0 .AND. MULP.LT.2 .and.idummy.lt.100 ) GO TO 240
      IF(IWMODE.EQ.0 .AND. MULP.LT.2 ) GO TO 200
CC-----------------------------------------------------------------------
C
C
C
C -- SET KIND OF PARTICLE
C
      CALL NEKINDPA(MUL)
      NUMNE=MUL+3
C
C
 6000 IF(MUL.LE.2) THEN
C
C
C -- SET MOMENTUM OF EACH HADRON (1) -- HAD RON MULTIPLICITY =2
C
         TOTMAS=0.
         DO 301 I=1,MUL
            CALL MCMASS(IPNE(I+3),AM(I))
            AM(I)=AM(I)*1.0E-3
            TOTMAS=TOTMAS+AM(I)
  301    CONTINUE
         IF(TOTMAS.LT.W)GO TO 401
         IEC=IEC+1
         IF(IEC.LT.MAXERR)GO TO 200
         IERR=1
         RETURN
C -- VECTOR GENERATION BY PHASE SPACE
  401    CALL VCPHSVCT(W,AM,MUL,PPHA,IERRP)
         IF(IERRP.EQ.0)GO TO 411
         IERR=1
         RETURN
C -- LORENTZ BOOST
  411    DO 311 I=1,MUL
            DO 321 J=1,3
  321       PBUF(J)=PPHA(J,I)
            BM=AM(I)
            CALL MCVECBST(PBUF,BM,EV,GMH)
            DO 331 J=1,3
  331       PNE(J,I+3)=PBUF(J)
  311    CONTINUE
         GO TO 4000
C
C
C
C -- SET MOMENTUM OF EACH HADRON
C                             T.KAJITA 06-NOV-85 MODIFICATION
C                                 ( INCLUDE DELTA )
C
      ELSE
         CONTINUE
         IPP=0
         DO 299 I=1,MULP
            IF ( IPNE(4).EQ.2212 .AND. IPNE(I+4).EQ. 211) IPP=I
            IF ( IPNE(4).EQ.2112 .AND. IPNE(I+4).EQ.-211) IPP=I
  299    CONTINUE
         IF(IPP.EQ.0) IPP=1
C
         IPBUFF=IPNE(IPP+4)
         DO 297 I=IPP,MULP-1
            IPNE(I+4)=IPNE(I+4+1)
  297    CONTINUE
         IPNE(MULP+4)=IPBUFF
C
         TOTMAS=0.
         WMEV=W*1000.
  789    AM(1) = RNMASS(WMEV,1232.,0.,115.)/1000.
         IF(AM(1).LT.(PMASS+.140)) GO TO 789
         TOTMAS = TOTMAS+AM(1)
         DO 300 I=2,MULP
            CALL MCMASS(IPNE(I+3),AM(I))
            AM(I)=AM(I)*1.0E-3
            TOTMAS=TOTMAS+AM(I)
  300    CONTINUE
         IF(TOTMAS.LT.W)GO TO 400
         IEC=IEC+1
         IF(IEC.LT.MAXERR)GO TO 200
         IERR=1
         RETURN
C
C -- VECTOR GENERATION BY PHASE SPACE
  400    MULPD=MUL-1
C+++++++++++++
         CALL VCPHSVCT(W,AM,MULPD,PPHA,IERRP)
C+++++++++++++
         IF(IERRP.EQ.0)GO TO 410
         IERR=1
         RETURN
C -- LORENTZ BOOST
  410    DO 310 I=1,MULPD
            DO 320 J=1,3
  320       PBUF(J)=PPHA(J,I)
            BM=AM(I)
            CALL MCVECBST(PBUF,BM,EV,GMH)
            DO 330 J=1,3
  330       PNE(J,I+3)=PBUF(J)
  310    CONTINUE
C
C + MAKE A DECAY OF DELTA
         CALL NEDCY2BD (AM(1), PMASS, 0.140, PABS)
         IPN=1
         CALL RNDIR(DIRTMP)
         PBUF(1)=PABS*DIRTMP(1)
         PBUF(2)=PABS*DIRTMP(2)
         PBUF(3)=PABS*DIRTMP(3)
C
         PTOTRS=SQRT(PNE(1,4)**2+PNE(2,4)**2+PNE(3,4)**2)
         DIRRS(1)=PNE(1,4)/PTOTRS
         DIRRS(2)=PNE(2,4)/PTOTRS
         DIRRS(3)=PNE(3,4)/PTOTRS
         BETARS = SQRT ( PTOTRS**2 / ( AM(1)**2 + PTOTRS**2 ) )
C                  BETA OF RESONANCE
         GAM=FNGAMFCT(BETARS)
C                  GAMMA FACTOR
C
 1700    IF(IPN.EQ.1)BM=PMASS
         IF(IPN.EQ.2)BM=0.140
         CALL MCVECBST(PBUF,BM,DIRRS,GAM)
         DO 1330 J=1,3
            IF(IPN.EQ.1)PNE(J,4)=PBUF(J)
            IF(IPN.EQ.2)PNE(J,MULP+4)=PBUF(J)
 1330    CONTINUE
         IF(IPN.EQ.1) THEN
            IPN=2
            PBUF(1)=-PABS*DIRTMP(1)
            PBUF(2)=-PABS*DIRTMP(2)
            PBUF(3)=-PABS*DIRTMP(3)
            GO TO 1700
         ELSE
         END IF
      END IF
C
C
C                                             09-NOV-85 T.KAJITA
C -- CHECK THE FORWARD AND BACKWARD PION MULTIPLICITY
 4000 NFPI=0
      NBPI=0
      CALL MCMASS(IPNE(MUL+3),AM(MUL))
      AM(MUL)=AM(MUL)*1.0E-3
      DO 5000 K=1,MULP
C          -- LORENTZ BOOST
         DO 2320 J=1,3
 2320    PBUF(J)=PNE(J,K+4)
         BM=AM(K+1)
         EVN(1)=-EV(1)
         EVN(2)=-EV(2)
         EVN(3)=-EV(3)
         CALL MCVECBST(PBUF,BM,EVN,GMH)
         PRIN=PBUF(1)*EV(1) + PBUF(2)*EV(2) + PBUF(3)*EV(3)
         IF(PRIN.GE.0.) NFPI=NFPI+1
         IF(PRIN.LT.0.) NBPI=NBPI+1
 5000 CONTINUE
CMS move 5050
CMS nechkfb is only for original multi pi production, not for jet set.
CMS M.Shiozawa  Mar.-5-1998
CMS 5050 CONTINUE
      CALL NECHKFB(W,MULP,NFPI,ICHK)
C                                        ICHK=1 OK,  =0 NO
      IF(ICHK.EQ.0) GO TO 6000
C
 5050 CONTINUE
C
      RETURN
C
C -- ERROR ROUTINE
C
 1000 IERR=1
      RETURN
      END
