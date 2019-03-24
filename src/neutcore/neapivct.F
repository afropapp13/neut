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
*       2016.07.20 ; C.Bronner  KNO multiplicity applies to number of hadrons and not pions
*                               Effect of multiplicity on W generation
*                               Remove delta mass when comparing mass of particles to W
*       2016.07.26 ; C.Bronner  Cleaning/commenting for new release
*
**********************************************************************
      IMPLICIT NONE

#include "vcwork.h"
#include "nework.h"
#include "necard.h"

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

      REAL*4    AM(9)

      real*4    E,W,WMEV,X,Y,q2


#include "posinnuc.h"


      REAL *4 PMASS  
      real*4 xmp,xmn
      PARAMETER (xmp    = 0.93827231)
      PARAMETER (xmn    = 0.93956563)
      
      real*4 PNRAT(3)
      real*4 PXSEC, NXSEC, R1320P, R1320N
      integer*4 MAXERR
      DATA MAXERR/10/

      real*4 FRFREE



      INTEGER IERJET
      Integer idummy

      integer*4 MULP,MULPD,MUL,NFPI,NBPI
      integer*4 IEC,IWAVE,itype,iwmode,IPN,IPP,ipbuff,iwmult
      real*4    ratio,ran,dum,prin

      integer*4 I,J,K,ierrr,ierrp,ichk

      real*4    pfabs,pabs,enucl,enstop,beta,betars,gm,gam,gmh
      real*4    elep,amlep,plep,cost,pphad,totmas,ptotrs
      real*4    bm

      real*4 NUCMASS,W2Max

     
      real*4 adzazd

*************************************************************************
*     Following functions don't  seem to be built in the library if they 
*     are not called somewhere in neutcore, so put a random call to them 
*     here so that nuceff can be built.
*     Should be fixed, but not sure how, it was not needed in 5.3.4 - CB
      CALL FNPKNO(2,1.5)
      adzazd=FNGAMFCT(1.5)
*************************************************************************    

      E=ENEUT
      IERR=0
      IEC=0

*     Determine if interactions occur on proton or neutron
*     For this compute ratio of cross-sections on each nucleon
      
C  Determine whether interaction is CC or NC      
      ITYPE=1
      IF(IABS(JMOD0).EQ.41.or.IABS(JMOD0).EQ.46)ITYPE=0

C  Compute cross section on neutron and on proton
      PXSEC=FNAPICRS(E,IPAR,ITYPE,1.,R1320P)
      NXSEC=FNAPICRS(E,IPAR,ITYPE,0,R1320N)
      If (abs(JMOD0).eq.21.or.abs(JMOD0).eq.41) then   ! --- multi-pi mode
         PXSEC=PXSEC*R1320P
         NXSEC=NXSEC*R1320N
      else   ! --- DIS/PYTHIA mode
         PXSEC=PXSEC*(1.-R1320P)
         NXSEC=NXSEC*(1.-R1320N)
      endif


C     Get the ratio for the target composition
C     Free hydrogen target should always be a proton
      if (NUMBNDP.eq.0 .and. NUMBNDN.eq.0) then
         RATIO = 1.0
      else
         RATIO=float(NUMBNDP+NUMFREP)*PXSEC/(float(NUMBNDN)*NXSEC+
     &         float(NUMBNDP+NUMFREP)*PXSEC)
      endif
      

C     FRFREE --- FREE PROTON/TOTAL PROTON     08-OCT-88 
      FRFREE  =float(NUMFREP)/float(NUMBNDP+NUMFREP)

#ifdef UPMU_RCK
      RATIO=0.5
      FRFREE  =0.0
#endif

C
C -- SET TARGET PARTICLE AND INTERACTION MODE
C
      MODENE=JMOD0
      IPNE(1)=IPAR
      IFLGNE(1)=-1
      RAN=RLU(DUM)
      IPNE(2)=2112
      IF(RAN.LT.RATIO)IPNE(2)=2212
      IFLGNE(2)=-1

C   Assign target mass       
      IF (IPNE(2).EQ.2112) THEN
         PMASS=xmn
      else
         PMASS=xmp
      endif
    
C   Determine if interaction occurs on free or bound proton
      IBOUND = 1
      RAN=RLU(DUM)
      IF(IPNE(2).EQ.2212 .AND. RAN.LT.FRFREE) IBOUND=0

C   Determine outgoing lepton/neutrino
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
	 IFLGNE(3)=2
      ENDIF
C
C -- SET NEUTRINO DIRECTION AND TARGET MOMENTUM
C
      DO 100 I=1,3
         PNE(I,1)=E*DIRNEU(I)
         PNE(I,2)=0.
         PNEU(I)=PNE(I,1)
  100 CONTINUE

C -- CONSIDER FERMI MOTION OF NUCLEON T.KAJITA  18-NOV-88 ------------

 96   PFABS=0.
C    Generate Fermi momentum for target nucleon      
      IF (IBOUND.EQ.1 .AND. ABS(IPAR).NE.16) THEN
         CALL nefrmmom(PFABS,IWAVE)
         PFABS=PFABS*1.0E-3
      ENDIF
      CALL RNDIR(DIRN)
      PNE(1,2)=DIRN(1)*PFABS
      PNE(2,2)=DIRN(2)*PFABS
      PNE(3,2)=DIRN(3)*PFABS
      
C    Neutrino is boosted to the rest frame of the target nucleon 
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

      IF(ENSTOP.LT.0.550 .AND. IABS(IPAR).EQ.12 .AND. IBOUND.EQ.1)
     &                                                         GO TO 96
      IF(ENSTOP.LT.0.650 .AND. IABS(IPAR).EQ.14 .AND. IBOUND.EQ.1)
     &                                                         GO TO 96
      E=ENSTOP


C  Interaction mode 
C     -- IWMODE=0 -> W<2.0GeV  => multi-pi mode
C     -- IWMODE=1 -> W>2.0GeV  => PYTHIA
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

      iwmult=0    !--- counter to avoid regenerating W forever
  250 iwmult=iwmult+1 
      iec=0
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
C  Rest of the program is multi-pi mode
      ELSE
C  Generate X and Y to have 1.3 GeV < W < 2.0 GeV        
         CALL RNW1320(IPAR,E,X,Y,W,ITYPE,IPNE(2),IERRR)
      END IF
C-----------------------------------------------------------------------
      IF(IERRR.EQ.1)GO TO 1000

*   Generate lepton kinematics
*   Regenerate W if impossible values for momentum and cosine
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
C -- SET MULTIPLICITY
C
* Counter to avoid infinite loop
      idummy=0
 240  idummy = idummy+1

* Generate hadron multiplicity
      MUL=NEMULTIP(W,IPAR,IPNE(2))  ! --- number of hadrons
      MULP=MUL-1                    ! --- number of pions
C === NEXT LINE WAS ADDED BY M.N ON 9-AUG-88
      IF(MULP.GE.8)MULP=7
      IF(MUL.GE.9)MUL=8

*  Multi-pi mode is only for events with at least two pions, one
*  particle production is handled by resonant modes.
*  If less than 2 pions, regenerate W, and not just the multiplicity
*  this is to take into account the dependance of multiplicity on W      
      IF(MULP.ge.2 ) GO TO 225
      IF(IWMULT.LT.100)GO TO 250

      IERR=1
      RETURN
CC-----------------------------------------------------------------------

C
C -- Assign particle content
C
 225  CONTINUE

      CALL NEKINDPA(MUL)
      NUMNE=MUL+3

 6000 CONTINUE

C
C -- SET MOMENTUM OF EACH HADRON
C        
 
*     Total mass of the hadronic system
      TOTMAS=0.
      DO 300 I=1,MUL
         CALL MCMASS(IPNE(I+3),AM(I))
         AM(I)=AM(I)*1.0E-3
         TOTMAS=TOTMAS+AM(I)
 300  CONTINUE
C     CB     Check W is enough to create the hadronic system         
      IF(TOTMAS.LT.W)GO TO 400
C     CB     If not, first try to regenerate the multiplicity for this W
      IF(idummy.lt.100)GO TO 240
C     CB     If proper multiplicity cannot be generated, rethrow W
      IEC=IEC+1
      IF(IEC.LT.MAXERR)GO TO 200
      
      IERR=1
      RETURN

400   CONTINUE

C     
C     -- VECTOR GENERATION BY PHASE SPACE
C
 
* Generate by phase space in the rest frame of hadronic system
      CALL VCPHSVCT(W,AM,MUL,PPHA,IERRP)
      IF(IERRP.EQ.0)GO TO 410
      IERR=1
      RETURN

* Boost back to experimental frame     
 410  DO 310 I=1,MUL
         DO 320 J=1,3
            PBUF(J)=PPHA(J,I)
 320     CONTINUE
         BM=AM(I)
         CALL MCVECBST(PBUF,BM,EV,GMH)
         DO 330 J=1,3
            PNE(J,I+3)=PBUF(J)
 330     CONTINUE
 310  CONTINUE

C     
C     09-NOV-85 T.KAJITA
C     -- CHECK THE FORWARD AND BACKWARD PION MULTIPLICITY
 4000 NFPI=0
      NBPI=0
      CALL MCMASS(IPNE(MUL+3),AM(MUL))
      AM(MUL)=AM(MUL)*1.0E-3
      DO 5000 K=1,MULP
C     -- LORENTZ BOOST
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
