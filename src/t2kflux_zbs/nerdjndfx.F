      SUBROUTINE nerdjndfx

      implicit none

#ifdef FLUX_10A
#include "uhdef.fh"
#else
#ifdef FLUX_10C
#include "uhdef_10c.fh"
#else
#ifdef FLUX_11A
#include "uhdef_11a.fh"
#else
#ifdef FLUX_11B
#include "uhdef_11b.fh"
#else
#ifdef FLUX_13
#include "uhdef_13_uwfunc.fh"
#else
#include "beamntpl.h"
#endif
#endif
#endif
#endif      
#endif

#include "beamvectbl.h"

#define ngFDMAX 12

      INTEGER IDATA(500)
      REAL    RDATA(500)
      EQUIVALENCE(IDATA(1),RDATA(1))
      CHARACTER*32 CNAME

      integer*4 IPBANK,iseg,ndata,i,iofst
      data ipbank/0/
C
C     -- MAKE EMB
C
      CNAME='JNDFXVEC'

      CALL KZBLOC(CNAME,IPBANK)
      IF (IPBANK.EQ.0) GOTO 1000

C     
C     -- SET DETECTOR INFORMATION 
C     
C        IND+ 1 : I : Energy of neutrino
C           + 2 : I : pid of the parent particle
C           + 3 : I : decay mode of the parent particle
C           + 4 : R : momentum of parent particle
C           + 5 : R : decay position in global coord. (x)
C           + 6 : R : decay position in global coord. (y)
C           + 7 : R : decay position in global coord. (z)
C           + 8 : R : direction of parent particle in global coord. (x)
C           + 9 : R : direction of parent particle in global coord. (y)
C           +10 : R : direction of parent particle in global coord. (z)
C           +11 : R : cosine to the beam direction of the parent particle
C           +12 : R : normalization factor
C           +13 : R : momentum of the parent particle (at the generation point)
C           +14 : R : generated position in global coord. (x)
C           +15 : R : generated position in global coord. (y)
C           +16 : R : generated position in global coord. (z)
C           +17 : R : direction of parent particle at generation point
C                                             in global coord. (x)
C           +18 : R : direction of parent particle at generation point
C                                             in global coord. (y)
C           +19 : R : direction of parent particle at generation point
C                                             in global coord. (z)
C           +20 : R : cosine to the beam direction of the parent particle
C                       at generation point
C           +21 : R : position of the neutrino at the detector 
C                                             in detector coord. (r)
C           +22 : R : position of the neutrino at the detector 
C                                             in detector coord. (x)
C           +23 : R : position of the neutrino at the detector 
C                                             in detector coord. (y)
C           +24 : R : position of the neutrino at the detector 
C                                             in global coord. (x)
C           +25 : R : position of the neutrino at the detector 
C                                             in global coord. (y)
C           +26 : R : position of the neutrino at the detector 
C                                             in global coord. (z)
C           +27 : I : ID of the front(near) detector

      ISEG=1
#if defined(FLUX_10C)||defined(FLUX_11A)||defined(FLUX_11B)||defined(FLUX_13)
      CALL KZGET0(CNAME,ISEG,NDATA,IDATA)
#else
      CALL KZGET1(CNAME,ISEG,NDATA,IDATA)
#endif

#ifdef FLUX_10A
      if (NDATA.NE.71) THEN
         write(*,*) 'Length of the bank is different from 71(10A-)'
         STOP
      endif
#else
#if defined(FLUX_10C)||defined(FLUX_11A)||defined(FLUX_11B)||defined(FLUX_13)
#else
      if (NDATA.NE.27) THEN
         write(*,*) 'Length of the bank is different from 27(-07A)'
         STOP
      endif
#endif
#endif      

      ENUFD     = RDATA( 1)
      PPIDFD    = IDATA( 2)
      MODEFD    = IDATA( 3)
      ppiFD     = RDATA( 4)
      xpiFD(1)  = RDATA( 5)
      xpiFD(2)  = RDATA( 6)
      xpiFD(3)  = RDATA( 7)
      npiFD(1)  = RDATA( 8)
      npiFD(2)  = RDATA( 9)
      npiFD(3)  = RDATA(10)
      cospibmFD = RDATA(11)
      normFD    = RDATA(12)
      ppi0FD    = RDATA(13)
      xpi0FD(1) = RDATA(14)
      xpi0FD(2) = RDATA(15)
      xpi0FD(3) = RDATA(16)
      npi0FD(1) = RDATA(17)
      npi0FD(2) = RDATA(18)
      npi0FD(3) = RDATA(19)
      cospi0bmFD= RDATA(20)
      rFD       = RDATA(21)
      xnuFD     = RDATA(22)
      ynuFD     = RDATA(23)
      nnuFD(1)  = RDATA(24)
      nnuFD(2)  = RDATA(25)
      nnuFD(3)  = RDATA(26)
      idFD      = IDATA(27)

#ifdef FLUX_10A
      nvtx0FD   = IDATA(28)
      gipartFD  = IDATA(29)
      gpos0FD(1)= RDATA(30)
      gpos0FD(2)= RDATA(31)
      gpos0FD(3)= RDATA(32)
      gvec0FD(1)= RDATA(33)
      gvec0FD(2)= RDATA(34)
      gvec0FD(3)= RDATA(35)
      gamom0FD  = RDATA(36)
      spidFD    = IDATA(37)
      pgenFD    = IDATA(38)
      psi0FD    = RDATA(39)
      xsi0FD(1) = RDATA(40)
      xsi0FD(2) = RDATA(41)
      xsi0FD(3) = RDATA(42)
      nsi0FD(1) = RDATA(43)
      nsi0FD(2) = RDATA(44)
      nsi0FD(3) = RDATA(45)
      cossi0bmFD= RDATA(46)
      xsiFD(1)  = RDATA(47)
      xsiFD(2)  = RDATA(48)
      xsiFD(3)  = RDATA(49)
      smechFD   = IDATA(50)
      intgtFD   = IDATA(51)
      prvtxFD(1)= RDATA(52)
      prvtxFD(2)= RDATA(53)
      prvtxFD(3)= RDATA(54)
      smedFD    = IDATA(55)
      gppidFD   = IDATA(56)
      xgpi0FD(1)= RDATA(57)
      xgpi0FD(2)= RDATA(58)
      xgpi0FD(3)= RDATA(59)
      xgpiFD(1) = RDATA(60)
      xgpiFD(2) = RDATA(61)
      xgpiFD(3) = RDATA(62)
      pgpi0FD   = RDATA(63)
      gpmechFD  = IDATA(64)
      gpmedFD   = IDATA(65)
      prmechFD  = IDATA(66)
      prmedFD   = IDATA(67)
      prdghtFD  = IDATA(68)
      sdghtFD   = IDATA(69)
      gpdghtFD  = IDATA(70)
      chainFD   = IDATA(71)
#endif
#if defined(FLUX_10C)||defined(FLUX_11A)||defined(FLUX_11B)||defined(FLUX_13)
      nvtx0FD   = IDATA(28)
      gipartFD  = IDATA(29)
      gpos0FD(1)= RDATA(30)
      gpos0FD(2)= RDATA(31)
      gpos0FD(3)= RDATA(32)
      gvec0FD(1)= RDATA(33)
      gvec0FD(2)= RDATA(34)
      gvec0FD(3)= RDATA(35)
      gamom0FD  = RDATA(36)
      EnuSKFD   = RDATA(37)
      normSKFD  = RDATA(38)
      anormFD   = RDATA(39)
      ngFD      = IDATA(40)
#endif
      if (ngFD.gt.ngFDMAX) then
         ngFD = ngFDMAX
      endif

      DO 100 I=1,ngFD
         gpxFD(I)   = RDATA(41+(I-1)*9)
         gpyFD(I)   = RDATA(42+(I-1)*9)
         gpzFD(I)   = RDATA(43+(I-1)*9)
         gcosbmFD(I)= RDATA(44+(I-1)*9)
         gvxFD(I)   = RDATA(45+(I-1)*9)
         gvyFD(I)   = RDATA(46+(I-1)*9)
         gvzFD(I)   = RDATA(47+(I-1)*9)
         gpidFD(I)  = IDATA(48+(I-1)*9)
         gmecFD(I)  = IDATA(49+(I-1)*9)
 100  continue

#if defined(FLUX_11A)||defined(FLUX_11B)||defined(FLUX_13)
      IOFST = 50 + (I-1)*9
      DO 110 I=1,ngFD
         gmatFD(I)   = IDATA(IOFST  +(I-1)*5)
         gdistcFD(I) = RDATA(IOFST+1+(I-1)*5)
         gdistalFD(I)= RDATA(IOFST+2+(I-1)*5)
         gdisttiFD(I)= RDATA(IOFST+3+(I-1)*5)
         gdistfeFD(I)= RDATA(IOFST+4+(I-1)*5)
 110  continue
#endif         

#if defined(FLUX_13)
      idExitFD    =IDATA(55+(ngFD-1)*9)
      ngExitFD    =IDATA(56+(ngFD-1)*9)
                   
      posExitFD(1)=IDATA(57+(ngFD-1)*9)
      posExitFD(2)=IDATA(58+(ngFD-1)*9)
      posExitFD(3)=IDATA(59+(ngFD-1)*9)
                   
      momExitFD(1)=IDATA(60+(ngFD-1)*9)
      momExitFD(2)=IDATA(61+(ngFD-1)*9)
      momExitFD(3)=IDATA(62+(ngFD-1)*9)

#endif
      
      RETURN
C     
C     ++ ERROR MESSAGE
C     
 1000 WRITE(6,900)
  900 FORMAT(' *** ERROR *** ( BANK READ IN NERDJNUFX )')
      RETURN
      END
      
      
