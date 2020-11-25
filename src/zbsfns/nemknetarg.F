
      SUBROUTINE NEMKNETARG

      IMPLICIT NONE

#include "skhead.h"
#include "nework.h"
#include "necard.h"
#include "nrcard.h"
#include "nefillver.h"
#include "neutmodel.h"
#include "neutparams.h"
#include "neutcrs.h"
#include "vcwork.h"
#include "posinnuc.h"
#include "rscons.h"

      integer*4 NSEGM
      PARAMETER (NSEGM=15)

      integer*4 LSEGMTBL(15)
      DATA LSEGMTBL/ 3, 6,11, 7, 5,28,28,10,10,11,
     $               5, 3, 1,25, 3 /

      integer*4 IDATA(500)
      real*4    RDATA(500)

      EQUIVALENCE(IDATA(1),RDATA(1))

      CHARACTER*32 CNAME

      integer*4 IERR,ISEGM,NDATA,I

      real*4    fntotpau
      external  fntotpau

      CNAME='NEUTINFO'

      ierr = 0
      CALL KZBCR0(CNAME,IERR)
      if (ierr.ne.0) then 
         goto 999
      endif

C--------------------------      
C- EMB
C-    + 1 : Length of this segment
C     + 2 : Number of segments
C-    + 3 : Length of 1st segment  : Versions of libraries
C     + 4 : Length of 2nd segment  : Flux & Geometry
C     + 5 : Length of 3rd segment  : Interaction and nuclear effects
C     + 6 : Length of 4th segment  : Reaction models
C     + 7 : Length of 5th segment  : Target information
C     + 8 : Length of 6th segment  : Cross-section weights ( neutrino )
C     + 9 : Length of 7th segment  : Cross-section weights ( anti neutrino )
C     +10 : Length of 8th segment  : Total cross-section
C     +11 : Length of 9th segment  : CCQE parameters
C     +12 : Length of 10th segment : Resonance production parameters
C     +13 : Length of 11th segment : coherent related parameters 
C     +14 : Length of 12th segment : DIS related parameters
C     +15 : Length of 13th segment : MEC related parameters
C     +16 : Length of 14th segment : general / pion nucl. effect parameters
C     +17 : Length of 15th segment : nucleon nuclear effect related parameters

      ISEGM   = 0
      NDATA   = 17
      IDATA(1)= NSEGM+2
      IDATA(2)= NSEGM
      DO 10 I=1,NSEGM
         IDATA(2+I) = LSEGMTBL(I)
 10   continue

      call kzrep0(CNAME,ISEGM,'I',NDATA,IDATA)

C--------------------------      
C- Segment 1
C    Version of libraries
C     + 1 : NeutCore
C     + 2 : Nuceff
C     + 3 : Nucleon rescattering
C     

      ISEGM    = 1
      NDATA    = LSEGMTBL(ISEGM)
      IDATA(1) = corev
      IDATA(2) = nucev
      IDATA(3) = nuccv

      call kzrep0(CNAME,ISEGM,'I',NDATA,IDATA)

C--------------------------      
C- Segment 2
C    Flux & Geometry
C     + 1 : SK Geometry
C     + 2 : Solar Activity on/off ( for ATMPD )
C     + 3 : Solar Activity param. ( for ATMPD )
C     + 4 : Neut Flux type        ( for ATMPD )
C     + 5 : Years                 ( for ATMPD )
C     + 6 : ITAUFLG               ( for ATMPD )

      ISEGM    = 2
      NDATA    = LSEGMTBL(ISEGM)
      IDATA(1) = SK_GEOMETRY

      DO 20 I=2,NDATA
         IDATA(2) = 0
 20   continue

      call kzrep0(CNAME,ISEGM,'I',NDATA,IDATA)
      
C--------------------------      
C- Segment 3
C    Interaction
C     + 1 : Fermi momentum flag
C     + 2 : Pauli blocking flag
C     + 3 : Nuclear effect flag
C     + 4 : Interaction mode
C     + 5 : Fermi surface momentum
C     + 6 : Maximum Fermi momentum
C     + 7 : Initial nuclear potential        
C     + 8 : Final nuclear potential
C     + 9 : Formation zone flag
C     +10 : Nucleon rescattering flag     
C     +11 : Factor for nucleon rescattering prob.
C      

      ISEGM    = 3
      NDATA    = LSEGMTBL(ISEGM)
      IDATA(1) = NEFRMFLG
      IDATA(2) = NEPAUFLG
      IDATA(3) = NENEFO16
      IDATA(4) = NEMODFLG
      RDATA(5) = PFSURF
      RDATA(6) = PFMAX
      RDATA(7) = VNUINI
      RDATA(8) = VNUFIN
      IDATA(9) = IFORMLEN
      IDATA(10)= NUCRESCAT
      RDATA(11)= XNUCFACT

      call kzrep0(CNAME,ISEGM,'I',NDATA,IDATA)
      
C--------------------------     
C- Segment 4
C       Reaction models
C     + 1 : QE model
C     + 2 : Single pi model
C     + 3 : Coherent pi model
C     + 4 : Deep Inelastic scattering model
C     + 5 : QE Axial form factor model
C     + 6 : QE Axial mass
C     + 7 : Single pi Axial mass
C

      ISEGM    = 4
      NDATA    = LSEGMTBL(ISEGM)
      IDATA(1) = MDLQE
      IDATA(2) = MDLSPI
      IDATA(3) = MDLCOH
      IDATA(4) = MDLDIS
      IDATA(5) = MDLQEAF      
      RDATA(6) = XMAQE
      RDATA(7) = XMASPI

      if (MODELDIS.ne.MDLDIS) then
         write(*,*) 'MODELDIS=',MODELDIS,' but MDLDIS=',MDLDIS
         stop
      endif
      if (MODELCOH.ne.MDLCOH) then
         write(*,*) 'MODELCOH=',MODELCOH,' but MDLCOH=',MDLCOH
         stop
      endif

      call kzrep0(CNAME,ISEGM,'I',NDATA,IDATA)

C--------------------------      
C- Segment 5
C       Target information
C     + 1 : # of bound neutron
C     + 2 : # of bound proton
C     + 3 : # of free proton
C     + 4 : # of nucleons in an atom
C     + 5 : Interaction on nucleus (1) or free proton (0)

      ISEGM    = 5
      NDATA    = LSEGMTBL(ISEGM)

      IDATA(1)=NUMBNDN
      IDATA(2)=NUMBNDP
      IDATA(3)=NUMFREP
      IDATA(4)=NUMATOM
      IDATA(5)=IBOUND

      call kzrep0(CNAME,ISEGM,'I',NDATA,IDATA)

C--------------------------      
C- Segment 6
C       Cross-section weights ( neutrino )
C     + 1 : 
C      |  : Cross-section weight for neutrino
C     +28 : 

      ISEGM    = 6
      NDATA    = LSEGMTBL(ISEGM)

      DO 30 I=1,NDATA
         RDATA(I)=CRSNEUT(I)
 30   continue

      call kzrep0(CNAME,ISEGM,'F',NDATA,IDATA)

C--------------------------      
C- Segment 7
C       Cross-section weights ( anti neutrino )
C     + 1 : 
C      |  : Cross-section weight for anti neutrino
C     +28 : 
C

      ISEGM    = 7
      NDATA    = LSEGMTBL(ISEGM)

      DO 40 I=1,NDATA
         RDATA(I)=CRSNEUTB(I)
 40   continue

      call kzrep0(CNAME,ISEGM,'F',NDATA,IDATA)
      
C--------------------------      
C- Segment 8
C       Total cross-section
C     + 1 : Energy of neutrino  ( GeV )
C     + 2 : Cross-section       (^10-38cm^2)


      ISEGM    = 8
      NDATA    = LSEGMTBL(ISEGM)

      RDATA(1)=sqrt(PNE(1,1)**2+PNE(2,1)**2+PNE(3,1)**2)
      if (RDATA(1).gt.0.) THEN
         RDATA(2)=fntotpau(ipne(1),rdata(1))
      else
         RDATA(1)=0.
         RDATA(2)=0.
      endif
      
      DO 50 I=3,NDATA
         RDATA(I) = 0.
 50   continue

      call kzrep0(CNAME,ISEGM,'F',NDATA,IDATA)

C--------------------------     
C- Segment 9
C       CCQE related parameters 
C     + 1 : QE model ( MDLQE )
C     + 2 : QE Axial form factor model (MDLQEAF)
C     + 3 : QE Axial Vector mass (XMAQE)
C     + 4 : QE Vector mass (XMVQE)
C     + 5 : Kappa parameter ( KAPP )
C     + 6 : Fermi surface momentum for spectral func ( PFSF )
C     + 7 : Second class current Vector part ( SCCFV )
C     + 8 : Second class current Axial Vector part ( SCCFA )
C     + 9 : Error term for pseudo scalar form factor ( FPQE )
C     +10 : Radiative correction ( IRADCORR )
      
      ISEGM    = 9
      NDATA    = LSEGMTBL(ISEGM)

      IDATA(1)=MDLQE
      IDATA(2)=MDLQEAF
      RDATA(3)=XMAQE
      RDATA(4)=XMVQE
      RDATA(5)=KAPP
      RDATA(6)=PFSF
      RDATA(7)=SCCFV
      RDATA(8)=SCCFA
      RDATA(9)=FPQE
      IDATA(10)=IRADCORR

      call kzrep0(CNAME,ISEGM,'F',NDATA,IDATA)

C--------------------------     
C- Segment 10
C       resonance production related parameters 
C     + 1 : resonance production model ( MDLSPI )
C     + 2 : 1pi production form factor ( NEIFF = IFFSPI )
C     + 3 : Rein-Sehgal form factor type ( NENRTYPE )
C     + 4 : (Graczyk-Sobczyk) 1pi production Axial Vector mass ( XMASPI )
C     + 5 : (Graczyk-Sobczyk) 1pi production Vector mass ( XMVSPI )
C     + 6 : (Graczyk-Sobczyk) 1pi production CA5 parameter ( RCA5ISPI )
C     + 7 : 1pi production background normalization parameter ( RBGSCLSPI )
C     + 8 : Rein-Sehgal resonance production Axial Vector mass ( XMARSRES )
C     + 9 : Rein-Sehgal resonance production Vector mass ( XMVRSRES )
C     +10 : Delta absorption ( pion less delta decay ) flag  ( IPILESSDCY )
C     +11 : Delta absorption ( pion less delta decay ) fraction ( RPILESSDCY )
C
      ISEGM    = 10
      NDATA    = LSEGMTBL(ISEGM)

      IDATA(1)=MDLSPI
      IDATA(2)=NEIFF
      IDATA(3)=NENRTYPE
      RDATA(4)=XMASPI
      RDATA(5)=XMVSPI
      RDATA(6)=RCA5ISPI
      RDATA(7)=RBGSCLSPI
      RDATA(8)=XMARSRES
      RDATA(9)=XMVRSRES
      IDATA(10)=IPILESSDCY
      RDATA(11)=RPILESSDCY

      call kzrep0(CNAME,ISEGM,'F',NDATA,IDATA)

C--------------------------     
C- Segment 11
C        coherent related parameters 
C     + 1 : coherent pi production model ( MDLCOH )
C     + 2 : coherent pi Axial Vector mass ( XMACOH )
C     + 3 : RAD0NU
C     + 4 : fA1COH
C     + 5 : fB1COH
C
      ISEGM    = 11
      NDATA    = LSEGMTBL(ISEGM)

      IDATA(1)=MDLCOH
      RDATA(2)=XMACOH
      RDATA(3)=RAD0NU
      RDATA(4)=fA1COH
      RDATA(5)=fB1COH

      call kzrep0(CNAME,ISEGM,'F',NDATA,IDATA)

C--------------------------     
C- Segment 12
C        DIS related parameters
C     + 1 : DIS production model ( MDLDIS )
C     + 2 : PDF model ( NEPDF )
C     + 3 : Bodek Yang correciton ( NEBODEK )

      ISEGM    = 12
      NDATA    = LSEGMTBL(ISEGM)

      IDATA(1)=MDLDIS
      IDATA(2)=NEPDF
      IDATA(3)=NEBODEK

      call kzrep0(CNAME,ISEGM,'F',NDATA,IDATA)

C--------------------------     
C- Segment 13
C        MEC related parameters
C     + 1 : DUMMY
C
      ISEGM    = 13
      NDATA    = LSEGMTBL(ISEGM)

      IDATA(1)=1

      call kzrep0(CNAME,ISEGM,'F',NDATA,IDATA)

C--------------------------     
C- Segment 14
C        general / pion nuclear effect related parameters
C     + 1 : Fermi surface momentum in GeV ( PFSURF )
C     + 2 : Maximum Fermi momentum in GeV ( PFMAX )
C     + 3 : Initial state potential ( VNUINI )
C     + 4 : Final state potential ( VNUINI )
C     + 5 : formation zone switch ( IFORMLEN )
C     + 6 : formation zone free parameter ( FZMU2 )
C     + 7 : Fermin motion flag ( NEFRMFLG )
C     + 8 : Pauli Blocking flag ( NEPAUFLG )
C     + 9 : Nuclear effect on/off flag ( NENEFO16 )
C     +10 : Low energy pion scattering model ( NENEFMODL )
C     +11 : High energy pion scattering model ( NENFMODH )
C     +12 : Kinematical model for high energy pion scattering ( NENEFKINH )
C     +13 : Nucleon ejection after pion absoption ( NEABSPIEMIT )
C     +14 : Pion QE-scattering ( P<500MeV/c ) correction factor (FEFQE)
C     +15 : Pion QE-scattering ( P>500MeV/c ) correction factor (FEFQEH)
C     +16 : Pion hadron prod.  ( P>500MeV/c ) correction factor (FEFINEL)
C     +17 : Pion absorption    ( P<500MeV/c ) correction factor (FEFABS)
C     +18 : Pion forward scat. ( P>500MeV/c ) correction factor (FEFCOH)
C     +19 : Pion CEX amplitude ( P<500MeV/c ) correction factor (FEFCX)
C     +20 : Pion CEX amplitude ( P>500MeV/c ) correction factor (FEFCXH)
C     +21 : Pion QE scat. Isotoropic decay frac. ( P>500MeV/c ) (FEFQEHF)
C     +22 : Pion forward scat. frac. ( P>500MeV/c ) (FEFFOHF)
C     +23 : Fraction of inelastic scat. incl. c-ex. (p>500MeV/c)(FEFCXHF)
C     +24 : pion scattering coulomb scattering flag on/off      (FEFCOUL)
C     +25 : Correction factor on overall mean free path         (FEFALL)
C
      ISEGM    = 14
      NDATA    = LSEGMTBL(ISEGM)

      RDATA(1)=PFSURF
      RDATA(2)=PFMAX
      RDATA(3)=VNUINI
      RDATA(4)=VNUFIN
      IDATA(5)=IFORMLEN
      RDATA(6)=FZMU2
      IDATA(7)=NEFRMFLG
      IDATA(8)=NEPAUFLG
      IDATA(9)=NENEFO16
      IDATA(10)=NENEFMODL
      IDATA(11)=NENEFMODH
      IDATA(12)=NENEFKINH
      IDATA(13)=NEABSPIEMIT
      RDATA(14)=FEFQE
      RDATA(15)=FEFQEH
      RDATA(16)=FEFINEL
      RDATA(17)=FEFABS
      RDATA(18)=FEFCOH
      RDATA(19)=FEFCX
      RDATA(20)=FEFCXH
      RDATA(21)=FEFQEHF
      RDATA(22)=FEFCOHF
      RDATA(23)=FEFCXHF
      RDATA(24)=FEFCOUL
      RDATA(25)=FEFALL

      call kzrep0(CNAME,ISEGM,'F',NDATA,IDATA)

C--------------------------     
C- Segment 15
C        nucleon nuclear effect related parameters
C     + 1 : Nucleon rescattering on/off ( NUCRESCAT ) 
C     + 2 : Nucleon rescattering factor ( XNUCFACT )
C     + 3 : Nucleon rescattering Pauli blocking - Local/Global ( NUCRESFLG )

      ISEGM    = 15
      NDATA    = LSEGMTBL(ISEGM)

      IDATA(1)=NUCRESCAT
      RDATA(2)=XNUCFACT
      IDATA(3)=NUCRESFLG

      call kzrep0(CNAME,ISEGM,'F',NDATA,IDATA)

 999  continue
      if (ierr.ne.0) then
         write(*,*) 'Failed to create BANK : ',CNAME
      endif

      return

      end
