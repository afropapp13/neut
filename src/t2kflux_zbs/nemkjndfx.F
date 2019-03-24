      SUBROUTINE nemkjndfx

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

      INTEGER IDATA(500)
      REAL    RDATA(500)
      EQUIVALENCE(IDATA(1),RDATA(1))
      CHARACTER*32 CNAME

      integer*4 NWMAST,LSEGM,NSEGM,ierr,iseg
C
C     -- MAKE EMB
C
      CNAME='JNDFXVEC'

C-- Version 

      
      NWMAST=1
      NSEGM=1
#ifdef FLUX_10A
      LSEGM=71
      CALL KZBCR1(CNAME,NWMAST,'1I',NSEGM,LSEGM,
     $     '1F 2I 23F 3I 7F 2I 11F 2I 3F 2I 7F 8I',IERR)
      IDATA(1)=2
#else
#if defined(FLUX_10C)||defined(FLUX_11A)||defined(FLUX_11B)||defined(FLUX_13)
      CALL KZBCR0(CNAME,IERR)
      IDATA(1)=3
#else
      LSEGM=27
      CALL KZBCR1(CNAME,NWMAST,'1I',NSEGM,LSEGM,
     $     '1F 2I 23F 1I',IERR)      
      IDATA(1)=1
#endif
#endif
      IF (IERR.NE.0) GOTO 1000

#if defined(FLUX_10C)||defined(FLUX_11A)||defined(FLUX_11B)||defined(FLUX_13)
      CALL KZREP0(CNAME,0,'I',1,IDATA)
#else
      CALL KZREP1(CNAME,0,NWMAST,IDATA)
#endif

C     
C     -- SET DETECTOR INFORMATION 
C     
C        IND+ 1 : R : Energy of neutrino
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
C                                             at generation point
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

      RDATA( 1)=ENUFD
      IDATA( 2)=PPIDFD
      IDATA( 3)=MODEFD
      RDATA( 4)=ppiFD
      RDATA( 5)=xpiFD(1)
      RDATA( 6)=xpiFD(2)
      RDATA( 7)=xpiFD(3)
      RDATA( 8)=npiFD(1)
      RDATA( 9)=npiFD(2)
      RDATA(10)=npiFD(3)
      RDATA(11)=cospibmFD
      RDATA(12)=normFD
      RDATA(13)=ppi0FD
      RDATA(14)=xpi0FD(1)
      RDATA(15)=xpi0FD(2)
      RDATA(16)=xpi0FD(3)
      RDATA(17)=npi0FD(1)
      RDATA(18)=npi0FD(2)
      RDATA(19)=npi0FD(3)
      RDATA(20)=cospi0bmFD
      RDATA(21)=rFD
      RDATA(22)=xnuFD
      RDATA(23)=ynuFD
      RDATA(24)=nnuFD(1)
      RDATA(25)=nnuFD(2)
      RDATA(26)=nnuFD(3)
      IDATA(27)=idfd
      
#if defined(FLUX_10A)
      IDATA(28)=nvtx0FD
      IDATA(29)=gipartFD
      RDATA(30)=gpos0FD(1)
      RDATA(31)=gpos0FD(2)
      RDATA(32)=gpos0FD(3)
      RDATA(33)=gvec0FD(1)
      RDATA(34)=gvec0FD(2)
      RDATA(35)=gvec0FD(3)
      RDATA(36)=gamom0FD
      IDATA(37)=spidFD
      IDATA(38)=pgenFD
      RDATA(39)=psi0FD
      RDATA(40)=xsi0FD(1)
      RDATA(41)=xsi0FD(2)
      RDATA(42)=xsi0FD(3)
      RDATA(43)=nsi0FD(1)
      RDATA(44)=nsi0FD(2)
      RDATA(45)=nsi0FD(3)
      RDATA(46)=cossi0bmFD
      RDATA(47)=xsiFD(1)
      RDATA(48)=xsiFD(2)
      RDATA(49)=xsiFD(3)
      IDATA(50)=smechFD
      IDATA(51)=intgtFD
      RDATA(52)=prvtxFD(1)
      RDATA(53)=prvtxFD(2)
      RDATA(54)=prvtxFD(3)
      IDATA(55)=smedFD
      IDATA(56)=gppidFD
      RDATA(57)=xgpi0FD(1)
      RDATA(58)=xgpi0FD(2)
      RDATA(59)=xgpi0FD(3)
      RDATA(60)=xgpiFD(1)
      RDATA(61)=xgpiFD(2)
      RDATA(62)=xgpiFD(3)
      RDATA(63)=pgpi0FD 
      IDATA(64)=gpmechFD
      IDATA(65)=gpmedFD
      IDATA(66)=prmechFD
      IDATA(67)=prmedFD
      IDATA(68)=prdghtFD
      IDATA(69)=sdghtFD
      IDATA(70)=gpdghtFD
      IDATA(71)=chainFD
#endif
#if defined(FLUX_10C)||defined(FLUX_11A)||defined(FLUX_11B)
      IDATA(28)=nvtx0FD
      IDATA(29)=gipartFD
      RDATA(30)=gpos0FD(1)
      RDATA(31)=gpos0FD(2)
      RDATA(32)=gpos0FD(3)
      RDATA(33)=gvec0FD(1)
      RDATA(34)=gvec0FD(2)
      RDATA(35)=gvec0FD(3)
      RDATA(36)=gamom0FD
      RDATA(37)=EnuSKFD
      RDATA(38)=normSKFD
      RDATA(39)=anormFD
      IDATA(40)=ngFD

#ifdef FLUX_11B
      if (ngFD.gt.ngFDMAX) then
         ngFD = ngFDMAX
      endif
#endif

      DO 100 I=1,ngFD
         RDATA(41+(I-1)*9)=gpxFD(I)
         RDATA(42+(I-1)*9)=gpyFD(I)
         RDATA(43+(I-1)*9)=gpzFD(I)
         RDATA(44+(I-1)*9)=gcosbmFD(I)
         RDATA(45+(I-1)*9)=gvxFD(I)
         RDATA(46+(I-1)*9)=gvyFD(I)
         RDATA(47+(I-1)*9)=gvzFD(I)
         IDATA(48+(I-1)*9)=gpidFD(I)
         IDATA(49+(I-1)*9)=gmecFD(I)
#if defined(FLUX_11A)||defined(FLUX_11B)
         IDATA(50+(I-1)*9)=gmatFD(I)   
         RDATA(51+(I-1)*9)=gdistcFD(I) 
         RDATA(52+(I-1)*9)=gdistalFD(I)
         RDATA(53+(I-1)*9)=gdisttiFD(I)
         RDATA(54+(I-1)*9)=gdistfeFD(I)
#endif
#endif
 100  continue

#if defined(FLUX_13)
      IDATA(55+(ngFD-1)*9)=idExitFD
      IDATA(56+(ngFD-1)*9)=ngExitFD

      RDATA(57+(ngFD-1)*9)=posExitFD(1)
      RDATA(58+(ngFD-1)*9)=posExitFD(2)
      RDATA(59+(ngFD-1)*9)=posExitFD(3)

      RDATA(60+(ngFD-1)*9)=momExitFD(1)
      RDATA(61+(ngFD-1)*9)=momExitFD(2)
      RDATA(62+(ngFD-1)*9)=momExitFD(3)
#endif         
      ISEG=1
      
#if defined(FLUX_10C)
      LSEGM = 49 + ngFD*9
      CALL KZREP0(CNAME,ISEG,'F',LSEGM,IDATA)
#else
#if defined(FLUX_11A)||defined(FLUX_11B)
      LSEGM = 54 + ngFD*9
      CALL KZREP0(CNAME,ISEG,'F',LSEGM,IDATA)
#else      
#if defined(FLUX_13)
      LSEGM = 49 + ngFD*9 +  8
      CALL KZREP0(CNAME,ISEG,'F',LSEGM,IDATA)
#else
      CALL KZREP1(CNAME,ISEG,LSEGM,IDATA(1))
#endif
#endif
#endif

      RETURN
C     
C     ++ ERROR MESSAGE
C     
 1000 WRITE(6,900)
  900 FORMAT(' *** ERROR *** ( BANK CREATION IN nemkjnufx )')
      RETURN
      END
      
      
