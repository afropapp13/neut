
      SUBROUTINE NERDNETARG

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

      integer*4 IERR,ISEGM,NDATA,I,NSEGMFILE

      CNAME='NEUTINFO'

      ierr = 0
      CALL KZBLOC(CNAME,IERR)
      if (ierr.eq.0) then 
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
C     +10 : Length of 8th segmetn  : Total cross-section
C

      ISEGM   = 0
      call kzget0(CNAME,ISEGM,NDATA,IDATA)

      DO 10 I=1,8
         IF ((i.ne.6).and.(i.ne.7)) THEN
            if (IDATA(2+I).ne.LSEGMTBL(I)) then
               write (*,*) CNAME,':Bank length for seg #',i,
     $              ' (',IDATA(2+I),') is different from expected:',
     $              ' (',LSEGMTBL(I),').'
               stop
            endif
         endif
 10   continue

C-- new version of the bank --
      IF (NDATA.GT.10) THEN
         DO 20 I=9,NDATA-2
            if (IDATA(2+I).ne.LSEGMTBL(I)) then
               write (*,*) CNAME,':Bank length for seg #',i,
     $              ' (',IDATA(2+I),') is different from expected:',
     $              ' (',LSEGMTBL(I),').'
               stop
            endif
 20      continue
      endif
      NSEGMFILE = NDATA-2

C--------------------------      
C- Segment 1
C    Version of libraries
C     + 1 : NeutCore
C     + 2 : Nuceff
C     + 3 : Nucleon rescattering
C     

      ISEGM = 1
      
      call kzget0(CNAME,ISEGM,NDATA,IDATA)
      if (NDATA.ne.LSEGMTBL(ISEGM)) then
         write (*,*) CNAME,':Bank length for seg #',ISEGM,
     $        ' (',NDATA,') is different from expected',
     $        ' (',LSEGMTBL(ISEGM),').'
         stop
      endif

      corev = IDATA(1)
      nucev = IDATA(2)
      nuccv = IDATA(3) 

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

      call kzget0(CNAME,ISEGM,NDATA,IDATA)
      if (NDATA.ne.LSEGMTBL(ISEGM)) then
         goto 9999
      endif

      SK_GEOMETRY = IDATA(1)

      
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

      call kzget0(CNAME,ISEGM,NDATA,IDATA)
      if (NDATA.ne.LSEGMTBL(ISEGM)) then
         goto 9999
      endif

      NEFRMFLG = IDATA(1) 
      NEPAUFLG = IDATA(2) 
      NENEFO16 = IDATA(3) 
      NEMODFLG = IDATA(4) 
      PFSURF   = RDATA(5) 
      PFMAX    = RDATA(6) 
      VNUINI   = RDATA(7) 
      VNUFIN   = RDATA(8) 
      IFORMLEN = IDATA(9) 
      NUCRESCAT= IDATA(10)
      XNUCFACT = RDATA(11)

C--------------------------     
C- Segment 4
C       Reaction models
C     + 1 : QE model
C     + 2 : Single pi model
C     + 3 : Coherent pi model
C     + 4 : Deep Inelastic scattering model
C     + 5 : QE potential model
C     + 6 : QE Axial mass
C     + 7 : Single pi Axial mass
C

      ISEGM    = 4

      call kzget0(CNAME,ISEGM,NDATA,IDATA)
      if (NDATA.ne.LSEGMTBL(ISEGM)) then
         goto 9999
      endif

      MDLQE   = IDATA(1)
      MDLSPI  = IDATA(2)
      MDLCOH  = IDATA(3)
      MDLDIS  = IDATA(4)
      MDLQEAF = IDATA(5)     
      XMAQE   = RDATA(6)
      XMASPI  = RDATA(7)

      MODELDIS = MDLDIS
      MODELCOH = MDLCOH

C--------------------------      
C- Segment 5
C       Target information
C     + 1 : # of bound neutron
C     + 2 : # of bound proton
C     + 3 : # of free proton
C     + 4 : # of nucleons in an atom

      ISEGM    = 5

      call kzget0(CNAME,ISEGM,NDATA,IDATA)
      if (NDATA.ne.LSEGMTBL(ISEGM)) then
         goto 9999
      endif

      NUMBNDN = IDATA(1)
      NUMBNDP = IDATA(2)
      NUMFREP = IDATA(3)
      NUMATOM = IDATA(4)
      IBOUND  = IDATA(5)

C--------------------------      
C- Segment 6
C       Cross-section weights ( neutrino )
C     + 1 : 
C      |  : Cross-section weight for neutrino
C     +28 : 

      ISEGM    = 6

      DO 25 I=1,500
         IDATA(I)=0
 25   continue

      call kzget0(CNAME,ISEGM,NDATA,IDATA)
C      if (NDATA.ne.LSEGMTBL(ISEGM)) then
C         goto 9999
C      endif

      DO 30 I=1,NDATA
         CRSNEUT(I) = RDATA(I)
 30   continue

C--------------------------      
C- Segment 7
C       Cross-section weights ( anti neutrino )
C     + 1 : 
C      |  : Cross-section weight for anti neutrino
C     +28 : 
C

      ISEGM    = 7

      DO 35 I=1,500
         IDATA(I)=0
 35   continue

      call kzget0(CNAME,ISEGM,NDATA,IDATA)
C      if (NDATA.ne.LSEGMTBL(ISEGM)) then
C         goto 9999
C      endif

      DO 40 I=1,NDATA
         CRSNEUTB(I) = RDATA(I)
 40   continue

C--------------------------      
C- Segment 8
C       Total cross-section
C     + 1 : Energy of neutrino  ( GeV )
C     + 2 : Cross-section       (^10-38cm^2)

      ISEGM    = 8

      call kzget0(CNAME,ISEGM,NDATA,IDATA)
      if (NDATA.ne.LSEGMTBL(ISEGM)) then
         goto 9999
      endif

      CRSENERGY = RDATA(1)
      TOTCRSNE  = RDATA(2)
      
      DO 50 I=3,NDATA
         DIFCRSNE(I-2) = RDATA(I)
 50   continue


      IF (NSEGMFILE.eq.8) then
         goto 99
      endif

C--------------------------     
C- Segment 9
C       CCQE related parameters 
C     + 1 : QE model ( MDLQE )
C     + 2 : QE Vector mass 
C     + 3 : QE Axial form factor model (MDLQEAF)
C     + 4 : QE Axial Vector mass (XMAQE)
C     + 5 : QE Vector mass (XMVQE)
C     + 6 : Kappa parameter ( KAPP )
C     + 7 : Fermi surface momentum for spectral func ( PFSF )
C     + 8 : Second class current Vector part ( SCCFV )
C     + 9 : Second class current Axial Vector part ( SCCFA )
C     +10 : Error term for pseudo scalar form factor ( FPQE )
C     +11 : Radiative correction ( IRADCORR )

      ISEGM    = 9

      call kzget0(CNAME,ISEGM,NDATA,IDATA)
      if (NDATA.ne.LSEGMTBL(ISEGM)) then
         goto 9999
      endif

      MDLQE   =IDATA(1)
      MDLQEAF =IDATA(2)
      XMAQE   =RDATA(3)
      XMVQE   =RDATA(4)
      KAPP    =RDATA(5)
      PFSF    =RDATA(6)
      SCCFV   =RDATA(7)
      SCCFA   =RDATA(8)
      FPQE    =RDATA(9)
      IRADCORR=IDATA(10)

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

      call kzget0(CNAME,ISEGM,NDATA,IDATA)
      if (NDATA.ne.LSEGMTBL(ISEGM)) then
         goto 9999
      endif

      MDLSPI    =IDATA(1)
      NEIFF     =IDATA(2)
      NENRTYPE  =IDATA(3)
      XMASPI    =RDATA(4)
      XMVSPI    =RDATA(5)
      RCA5ISPI  =RDATA(6)
      RBGSCLSPI =RDATA(7)
      XMARSRES  =RDATA(8)
      XMVRSRES  =RDATA(9)
      IPILESSDCY=IDATA(10)
      RPILESSDCY=RDATA(11)

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

      call kzget0(CNAME,ISEGM,NDATA,IDATA)
      if (NDATA.ne.LSEGMTBL(ISEGM)) then
         goto 9999
      endif

      MDLCOH=IDATA(1)
      XMACOH=RDATA(2)
      RAD0NU=RDATA(3)
      fA1COH=RDATA(4)
      fB1COH=RDATA(5)

C--------------------------     
C- Segment 12
C        DIS related parameters
C     + 1 : DIS production model ( MDLDIS )
C     + 2 : PDF model ( NEPDF )
C     + 3 : Bodek Yang correciton ( NEBODEK )

      ISEGM    = 12

      call kzget0(CNAME,ISEGM,NDATA,IDATA)
      if (NDATA.ne.LSEGMTBL(ISEGM)) then
         goto 9999
      endif

      MDLDIS =IDATA(1)
      NEPDF  =IDATA(2)
      NEBODEK=IDATA(3)

C--------------------------     
C- Segment 13
C        MEC related parameters
C     + 1 : DUMMY
C
      ISEGM    = 13

      call kzget0(CNAME,ISEGM,NDATA,IDATA)
      if (NDATA.ne.LSEGMTBL(ISEGM)) then
         goto 9999
      endif

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

      call kzget0(CNAME,ISEGM,NDATA,IDATA)
      if (NDATA.ne.LSEGMTBL(ISEGM)) then
         goto 9999
      endif

      PFSURF     =RDATA(1)
      PFMAX      =RDATA(2)
      VNUINI     =RDATA(3)
      VNUFIN     =RDATA(4)
      IFORMLEN   =IDATA(5)
      FZMU2      =RDATA(6)
      NEFRMFLG   =IDATA(7)
      NEPAUFLG   =IDATA(8)
      NENEFO16   =IDATA(9)
      NENEFMODL  =IDATA(10)
      NENEFMODH  =IDATA(11)
      NENEFKINH  =IDATA(12)
      NEABSPIEMIT=IDATA(13)
      FEFQE      =RDATA(14)
      FEFQEH     =RDATA(15)
      FEFINEL    =RDATA(16)
      FEFABS     =RDATA(17)
      FEFCOH     =RDATA(18)
      FEFCX      =RDATA(19)
      FEFCXH     =RDATA(20)
      FEFQEHF    =RDATA(21)
      FEFCOHF    =RDATA(22)
      FEFCXHF    =RDATA(23)
      FEFCOUL    =RDATA(24)
      FEFALL     =RDATA(25)

C--------------------------     
C- Segment 15
C        nucleon nuclear effect related parameters
C     + 1 : Nucleon rescattering on/off ( NUCRESCAT ) 
C     + 2 : Nucleon rescattering factor ( XNUCFACT )
C     + 3 : Nucleon rescattering Pauli blocking - Local/Global ( NUCRESFLG )

      ISEGM    = 15

      call kzget0(CNAME,ISEGM,NDATA,IDATA)
      if (NDATA.ne.LSEGMTBL(ISEGM)) then
         goto 9999
      endif

      NUCRESCAT = IDATA(1)
      XNUCFACT  = RDATA(2)
      NUCRESFLG = IDATA(3)

 99   return
      
 999  continue
      if (ierr.ne.0) then
         write(*,*) 'Failed to create BANK : ',CNAME
      endif

      return

 9999 continue
      write (*,*) CNAME,':Bank length for seg #',ISEGM,
     $     ' (',NDATA,') is different from expected',
     $     ' (',LSEGMTBL(ISEGM),').'
      stop
      
      end
