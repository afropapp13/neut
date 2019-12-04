**********************************************************************
*     -------------------------
      DOUBLE PRECISION FUNCTION EVDIFCRS()
*     -------------------------
*
*     ( purpose )
*        SAVE INTERACTION CROSS-SECTION
*
*     ( input )
*       COMMON NEWORK
*          modene
*          numne (=4)
*          ipne
*          pne
*       COMMON NECARD
*       COMMON NEUTMODEL
*       COMMON NEUTPARAMS
*
*     ( output )
*       DIFCRSSECT
*
*     ( creation date and author )
*       2010.07.19 ; Taken from nemodsel
*       2010.10.13 ; Integrated into NEUTCORE - P. de Perio
*
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C    MODSEL   : INTERACTION MODE
C
C
C        ######### NEUTRINO MODE #########
C
C            ***** CHARGED CURRENT *****
C
C               -- ELASTIC --
C           1 : NEU,N --> LEPTON-,P
C
C               -- SINGLE PI FROM DELTA RESONANCE --
C          11 : NEU,P --> LEPTON-,P,PI+
C          12 : NEU,N --> LEPTON-,P,PI0
C          13 : NEU,N --> LEPTON-,N,PI+
C
C          16 : NEU,O(16) --> LEPTON-,O(16),PI+
C
C               -- SINGLE GAMMA FROM DELTA RESONANCE --
C          17 : NEU,N --> LEPTON-,P,GAMMA
C
C               -- MULTI PI (1.3 < W < 2.0 GeV) --
C          21 : NEU,(N OR P) --> LEPTON-,(N OR P),MULTI PI
C
C               -- SINGLE ETA FROM DELTA RESONANCE --
C                                     (added 97/12/01 J.Kameda)
C          22 : NEU,N --> LEPTON-,P,ETA0
C
C               -- SINGLE K FROM DELTA RESONANCE --
C                                     (added 98/02/25 J.Kameda)
C          23 : NEU,N --> LEPTON-,LAMBDA,K+
C
C               -- DEEP INELASTIC (2.0 GeV < W , JET set) --
C          26 : NEU,(N OR P) --> LEPTON-,(N OR P),MESONS
C
C            ***** NEUTAL CURRENT *****
C
C               -- SINGLE PI FROM DELTA RESONANCE --
C          31 : NEU,N --> NEU,N,PI0
C          32 : NEU,P --> NEU,P,PI0
C          33 : NEU,N --> NEU,P,PI-
C          34 : NEU,P --> NEU,N,PI+
C
C          36 : NEU,O(16) --> NEU,O(16),PI0
C
C               -- SINGLE GAMMA FROM DELTA RESONANCE --
C          38 : NEU,N --> NEU,N,GAMMA
C          39 : NEU,P --> NEU,P,GAMMA
C
C               -- MULTI PI (1.3 GeV < W < 2.0 GeV) --
C          41 : NEU,(N OR P) --> NEU,(N OR P),MULTI PI
C
C               -- SINGLE ETA FROM DELTA RESONANCE --
C                                     (added 97/12/01 J.Kameda)
C          42 : NEU,N --> NEU,N,ETA0
C          43 : NEU,P --> NEU,P,ETA0
C
C               -- SINGLE  K  FROM DELTA RESONANCE --
C                                     (added 98/02/20 J.Kameda)
C          44 : NEU,N --> NEU,LAMBDA,K0
C          45 : NEU,P --> NEU,LAMBDA,K+
C
C               -- DEEP INELASTIC (2.0 GeV < W , JET set) --
C          46 : NEU,(N OR P) --> NEU,(N OR P),MESONS
C
C               -- ELASTIC --
C          51 : NEU,P --> NEU,P
C          52 : NEU,N --> NEU,N
C
C
C        ######### ANTI NEUTRINO MODE #########
C
C            ***** CHARGED CURRENT *****
C
C               -- ELASTIC --
C          -1 : NEUBAR,P --> LEPTON+,N
C
C               -- SINGLE PI FROM DELTA RESONANCE --
C         -11 : NEUBAR,N --> LEPTON+,N,PI-
C         -12 : NEUBAR,P --> LEPTON+,N,PI0
C         -13 : NEUBAR,P --> LEPTON+,P,PI-
C
C         -16 : NEUBAR,O(16) --> LEPTON+,O(16),PI-
C
C               -- SINGLE GAMMA FROM DELTA RESONANCE --
C         -17 : NEUBAR,P --> LEPTON+,N,GAMMA
C
C               -- MULTI PI (W > 1.4 GEV) --
C         -21 : NEUBAR,(N OR P) --> LEPTON+,(N OR P),MULTI PI
C
C               -- SINGLE ETA FROM DELTA RESONANCE --
C                                     (added 97/12/01 J.Kameda)
C         -22 : NEUBAR,P --> LEPTON+,N,ETA0
C
C               -- SINGLE  K FROM DELTA RESONANCE --
C                                     (added 98/02/25 J.Kameda)
C         -23 : NEUBAR,P --> LEPTON+,LAMBDA,K0
C
C               -- DEEP INELASTIC (2.0 GeV < W , JET set) --
C         -26 : NEUBAR,(N OR P) --> LEPTON+,(N OR P),MESONS
C
C               ** NEUTAL CURRENT **
C
C               -- SINGLE PI FROM DELTA RESONANCE --
C         -31 : NEUBAR,N --> NEUBAR,N,PI0
C         -32 : NEUBAR,P --> NEUBAR,P,PI0
C         -33 : NEUBAR,N --> NEUBAR,P,PI-
C         -34 : NEUBAR,P --> NEUBAR,N,PI+
C
C         -36 : NEUBAR,O(16) --> NEUBAR,O(16),PI0
C
C               -- SINGLE GAMMA FROM DELTA RESONANCE --
C         -38 : NEUBAR,N --> NEUBAR,N,GAMMA
C         -39 : NEUBAR,P --> NEUBAR,P,GAMMA
C
C               -- MULTI PI (W > 1.4 GEV) --
C         -41 : NEUBAR,(N OR P) --> NEUBAR,(N OR P),MULTI PI
C
C               -- SINGLE ETA FROM DELTA RESONANCE --
C                                     (added 97/12/01 J.Kameda)
C         -42 : NEUBAR,N --> NEUBAR,N,ETA0
C         -43 : NEUBAR,P --> NEUBAR,P,ETA0
C
C               -- SINGLE  K  FROM DELTA RESONANCE --
C                                     (added 98/02/20 J.Kameda)
C          -44 : NEUBAR,N --> NEUBAR,LAMBDA,K0
C          -45 : NEUBAR,P --> NEUBAR,LAMBDA,K+
C
C               -- DEEP INELASTIC (2.0 GeV < W , JET set) --
C         -46 : NEUBAR,(N OR P) --> NEUBAR,(N OR P),MESONS
C
C               -- ELASTIC --
C         -51 : NEUBAR,P --> NEUBAR,P
C         -52 : NEUBAR,N --> NEUBAR,N
C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C    IPAR : NEUTRINO TYPE
C        = 12  NUE
C        =-12  NUEBAR
C        = 14  NUMU
C        =-14  NUMUBAR
C    E    : ENERGY OF NEUTRINO (UNIT : GEV)
C
C    MEANING OF DATA
C    RATN   : NEUTRON RATIO IN TARGET
C    RATP   : PROTON RATIO IN TARGET
C    RATFRE : FREE PROTON RATIO IN TARGET
C    IFREE  : =0 BOUND NUCLEON, =1 FREE PROTON
C
C    IBOUND : =0 Free proton  , =1 Bound nucleon
C
      IMPLICIT NONE

#include "nework.h"
#include "rscons.h"
#include "neutcrs.h"
#include "neutmodel.h"
#include "necard.h"
#include "vcwork.h"
#include "posinnuc.h"

      integer*4 ipar,i
      real*4    e
C      integer*4 mdlqe
c      integer*4 IBOUND

C
C -- Preparation
C
      real*4 plep,   elep,   xmlep
      real*4 xmn01,  xmn02,  enuc1, enuc2
      real*4 q2,     coslep
      real*4 w,      x,      xcoh,  y
      real*8 e8,     x8,     y8,    q28

      real*4 dcp,    dc33p
      real*4 dcm,    dc33m

C -- For boosting
      real*4 PNEU(3),PFABS,AMIN,BETA,EV(3),GM

      real*8   qedifcrs
      external qedifcrs
      real*8   dndxdy
      external dndxdy
      real*8   cohma
      external cohma

      double precision LH_GeV
      external LH_GeV

      real*4   dcohnerein
      external dcohnerein

      real*4   dcohneberger
      external dcohneberger

      real*8 dnelsq2
      external dnelsq2
      real*8 dnncelsq2
      external dnncelsq2

C -- Parameters
      integer*4 iorig

C -- Flags
      Integer*4 ccnc

      integer*4 IMOD

C -- for SF nuclei tests
      logical HaveSpecFuncForNuclei, HaveEffSpecFuncForNuclei
      external HaveSpecFuncForNuclei, HaveEffSpecFuncForNuclei
C -- for SF parameters
      double precision k_GeV(4)
      double precision kprime_GeV(4)
      double precision p_GeV(4)
      double precision qTilde_GeV(4)
      logical nc
      integer*4 nucPdg

      evdifcrs = -1

      ipar = ipne(1)

      IF (MODENE.eq.11)  IMOD=1
      IF (MODENE.eq.12)  IMOD=2
      IF (MODENE.eq.13)  IMOD=3
      IF (MODENE.eq.31)  IMOD=6
      IF (MODENE.eq.32)  IMOD=4
      IF (MODENE.eq.33)  IMOD=7
      IF (MODENE.eq.34)  IMOD=5
      IF (MODENE.eq.-11) IMOD=11
      IF (MODENE.eq.-12) IMOD=12
      IF (MODENE.eq.-13) IMOD=13
      IF (MODENE.eq.-31) IMOD=16
      IF (MODENE.eq.-32) IMOD=14
      IF (MODENE.eq.-33) IMOD=17
      IF (MODENE.eq.-34) IMOD=15

      IF (MODENE.eq.22)  IMOD=1
      IF (MODENE.eq.42)  IMOD=11
      IF (MODENE.eq.43)  IMOD=12
      IF (MODENE.eq.-22) IMOD=2
      IF (MODENE.eq.-42) IMOD=13
      IF (MODENE.eq.-43) IMOD=14

      IF (MODENE.eq.23)  IMOD=1
      IF (MODENE.eq.44)  IMOD=11
      IF (MODENE.eq.45)  IMOD=12
      IF (MODENE.eq.-23) IMOD=2
      IF (MODENE.eq.-44) IMOD=13
      IF (MODENE.eq.-45) IMOD=14

      IF (MODENE.eq.17)  IMOD=1
      IF (MODENE.eq.38)  IMOD=11
      IF (MODENE.eq.39)  IMOD=12
      IF (MODENE.eq.-17) IMOD=2
      IF (MODENE.eq.-38) IMOD=13
      IF (MODENE.eq.-39) IMOD=14

      if (abs(ipne(3)).eq.11) then
         xmlep=XME
      else if (abs(ipne(3)).eq.13) then
         xmlep=XMMU
      else if (abs(ipne(3)).eq.15) then
         xmlep=XMTAU
      else
         xmlep=0.
      endif

      if (ipne(2).eq.2112) then
         xmn01 = XMNE
      else
         xmn01 = XMP
      endif
      if (ipne(4).eq.2112) then
         xmn02 = XMNE
      else
         xmn02 = XMP
      endif

      e = sqrt(pne(1,1)**2+pne(2,1)**2+pne(3,1)**2)

      enuc1= sqrt((pne(1,2)**2+pne(2,2)**2+pne(3,2)**2)+xmn01**2)
      plep = sqrt((pne(1,3)**2+pne(2,3)**2+pne(3,3)**2))
      elep = sqrt((pne(1,3)**2+pne(2,3)**2+pne(3,3)**2)+xmlep**2)
      enuc2= sqrt((pne(1,4)**2+pne(2,4)**2+pne(3,4)**2)+xmn02**2)

cRT: variables for SF
      k_GeV(1) = e
      k_GeV(2) = pne(1,1)
      k_GeV(3) = pne(2,1)
      k_GeV(4) = pne(3,1)

      kPrime_GeV(1) = elep
      kPrime_GeV(2) = pne(1,3)
      kPrime_GeV(3) = pne(2,3)
      kPrime_GeV(4) = pne(3,3)

      p_GeV(1) = enuc1
      p_GeV(2) = pne(1,2)
      p_GeV(3) = pne(2,2)
      p_GeV(4) = pne(3,2)

      qTilde_GeV(1) = enuc2-enuc1
      qTilde_GeV(2) = pne(1,1)-pne(1,3)
      qTilde_GeV(3) = pne(2,1)-pne(2,3)
      qTilde_GeV(4) = pne(3,1)-pne(3,3)

      nucPdg = ipne(2)
cRT: end variables for SF

      q2 =   (e - elep)**2
     $     - ( (pne(1,1)-pne(1,3))**2
     $        +(pne(2,1)-pne(2,3))**2
     $        +(pne(3,1)-pne(3,3))**2)

      w = sqrt ( (e + enuc1 - elep)**2
     $          -( (pne(1,1)+pne(1,2)-pne(1,3))**2
     $            +(pne(2,1)+pne(2,2)-pne(2,3))**2
     $            +(pne(3,1)+pne(3,2)-pne(3,3))**2))

      coslep = (  pne(1,1)*pne(1,3)
     $           +pne(2,1)*pne(2,3)
     $           +pne(3,1)*pne(3,3))/(e*plep)

      X = -1. * q2 / (2 * XMN01 * ( E - Elep ))
      Y = (E-Elep)/e
      xcoh = -q2/(2*0.938*E*Y)

      E8=e
      X8=X
      Y8=Y
      Q28=q2

      if (abs(modene).gt.30) then
         IORIG = 1
      else
         IORIG = 0
      endif

      nc = (abs(modene).gt.30)
C
C
C --  FOR NEUTRINO
C
      if (abs(modene).eq.1) then

         ccnc = 1

         if (IBOUND .eq. 1) then
            if ( (mod(mdlqe,1000)/100 .EQ. 4) .AND.
     $           HaveSpecFuncForNuclei() ) then
              evdifcrs = LH_GeV(k_GeV,kPrime_GeV,p_GeV,qTilde_GeV,IPAR,nc,nucPdg, NUMATOM)
            else if (((mod(mdlqe,1000)/100 .EQ. 6 ) .or.
     $                (mod(mdlqe,1000)/100 .EQ. 7 )      ).AND.
     &           HaveEffSpecFuncForNuclei() ) then
               evdifcrs
     $              = LH_GeV(k_GeV,kPrime_GeV,p_GeV,qTilde_GeV,IPAR,nc,nucPdg, NUMATOM)
            else
               evdifcrs
     $              = qedifcrs(ipar, ccnc, e, elep, coslep)
            endif
         else
            evdifcrs = dnelsq2(REAL(e,8),ipar,REAL(q2,8))
         endif
c         if (mdlqe.lt.400) then
c            write(*,*) "model less than 400"
c            evdifcrs
C     $        = qedifcrs(ipar, ccnc, e*1000., elep*1000., coslep)
c     $           = qedifcrs(ipar, ccnc, e, elep, coslep)
C         write(*,*) "New   :",evdifcrs
C         write(*,*) "Old   :",
C     $        ddifcrs(e*1000.,coslep,elep*1000.,ipar)
C         write(*,*) "ratio :",
C     $        evdifcrs/ddifcrs(e*1000.,coslep,elep*1000.,ipar)
c         else
c            write(*,*) "model greater than 400"
c         evdifcrs
c     $        = LH_GeV(k_GeV,kPrime_GeV,p_GeV,qTilde_GeV,IPAR,nc,nucPdg)
c     endif
      else if (    ((abs(modene).ge.11).and.(abs(modene).le.13))
     $         .or.((abs(modene).ge.31).and.(abs(modene).le.34)))
     $        then
         if ( w.le.XMN+XMPI ) then
            WRITE(*,*) "Resetting W: ", w, " -> ", XMN+XMPI
            w = XMN+XMPI
         endif
C     Boost neutrino to nucleon rest frame
        PFABS=SQRT(pne(1,2)**2+pne(2,2)**2+pne(3,2)**2)
        IF (PFABS.NE.0.) THEN
          BETA=PFABS/enuc1
          GM=1./SQRT(1.-BETA**2)
          do i =1,3
              PNEU(i) = pne(i,1)
              EV(i)=-pne(i,2)/PFABS
          enddo
          CALL MCVECBST(PNEU,0.,EV,GM)
          e = SQRT(PNEU(1)**2+PNEU(2)**2+PNEU(3)**2)
        ENDIF

         call rsdcrs(IMOD,IORIG,1, xmlep,e,q2,w,dcp,dc33p);
         call rsdcrs(IMOD,IORIG,-1,xmlep,e,q2,w,dcm,dc33m);
         evdifcrs = dcp + dcm
      else if ((abs(modene).eq.21).or.(abs(modene).eq.26)) then
         evdifcrs = dndxdy(IPAR,IPNE(2),1,E8,X8,Y8)
      else if ((abs(modene).eq.41).or.(abs(modene).eq.46)) then
         evdifcrs = dndxdy(IPAR,IPNE(2),0,E8,X8,Y8)
      else if ((abs(modene).eq.51).or.(abs(modene).eq.52)) then
         ccnc   = 0
         if (IBOUND .eq. 1) then


            if ( (mod(mdlqe,1000)/100 .EQ. 100) .AND.
     $            HaveSpecFuncForNuclei() ) then
               evdifcrs
     $              = LH_GeV(k_GeV,kPrime_GeV,p_GeV,qTilde_GeV,IPAR,nc,nucPdg, NUMATOM)
            else if ( ((mod(mdlqe,1000)/100 .EQ. 6 ).or.
     $                 (mod(mdlqe,1000)/100 .EQ. 7 )    ) .AND.
     &                 HaveEffSpecFuncForNuclei() ) then
               evdifcrs
     $              = LH_GeV(k_GeV,kPrime_GeV,p_GeV,qTilde_GeV,IPAR,nc,nucPdg, NUMATOM)
            else
               evdifcrs
     $              = qedifcrs(ipar, ccnc, e, elep, coslep)
            endif
         else
C RT: this is for free protons; will need to amend when free
C RT: neutrons are available
            write(*,*) "Free proton NCEL"
            evdifcrs = dnncelsq2(e,ipar,2212,q2)
         endif
c         evdifcrs
C     $        = qedifcrs(ipar, ccnc, e*1000., elep*1000., coslep)
c     $        = qedifcrs(ipar, ccnc, e, elep, coslep)
      else if (abs(modene).eq.16.or.abs(modene).eq.36) then
C        write(*,*) 'in evdifcrs, MDLCOH = ',MDLCOH
        if (MDLCOH.eq.0) then
          evdifcrs = dcohnerein(e,ipar,modene,CRSX,CRSY,CRSZ,CRSPHI)
C          write(*,*) 'in RS if statement, e=',e,
C     $ ' ipar=',ipar,' x=',CRSX,' y=',CRSY,' z=',CRSZ,' phi=',CRSPHI
       else if (MDLCOH.eq.2) then
          evdifcrs = dcohneberger(e,ipar,modene,CRSX,CRSY,CRSZ,CRSPHI)
C          write(*,*) 'in BS if statement, evdifcrs = ',evdifcrs,' e=',e,
C     $ ' ipar=',ipar,' x=',CRSX,' y=',CRSX,' z=',CRSX,' phi=',CRSPHI
       endif
      else if (   (abs(modene).eq.22)
     $        .or.(abs(modene).eq.42).or.(abs(modene).eq.43))
     $        then
         if ( w.le.XMN+XMETA ) then
            w = XMN+XMETA
         endif
         call rsdcrse(IMOD,IORIG,1, xmlep,e,q2,w,dcp,dc33p);
         call rsdcrse(IMOD,IORIG,-1,xmlep,e,q2,w,dcm,dc33m);
         evdifcrs = dcp + dcm
        else if (   (abs(modene).eq.23)
     $        .or.((abs(modene).eq.44).or.(abs(modene).eq.45)))
     $        then
         if ( w.le.XMLAMD+XMKAON ) then
            w = XMLAMD+XMKAON
         endif
         call rsdcrsk(IMOD,IORIG,1, xmlep,e,q2,w,dcp,dc33p);
         call rsdcrsk(IMOD,IORIG,-1,xmlep,e,q2,w,dcm,dc33m);
         evdifcrs = dcp + dcm
      else if (    (abs(modene).eq.17)
     $         .or.((abs(modene).eq.38).or.(abs(modene).eq.39)))
     $        then
         call rsdcrsg(IMOD,IORIG,1, xmlep,e,q2,w,dcp,dc33p);
         call rsdcrsg(IMOD,IORIG,-1,xmlep,e,q2,w,dcm,dc33m);
         evdifcrs = dcp + dcm
      else
         write(*,*) 'Unknown interaction mode',modene
C         stop
      endif

c      write(*,*) "evdifcrs: ", evdifcrs

      RETURN
      END
