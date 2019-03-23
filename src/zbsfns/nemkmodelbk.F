***********************************************************************
*     ------------------------
      SUBROUTINE NEMKMODELBK
*     ------------------------
*     
*     ( purpose )
*       Create NEUT bank for reaction models etc...
*
*     ( input ) 
*       COMMON block
*
*     ( output )
*       ZBS bank
*
*     ( creation date and author )
*       2007.08.22 G.Mitsuka 
*                  makes new bank for flag of reaction models
*
***********************************************************************
      Implicit none

#include "necard.h"
#include "skhead.h"
#include "nrcard.h"
#include "neutparams.h"
#include "neutmodel.h"
#include "nefillver.h"
#include "rscons.h"

      INTEGER*4 IERR
      INTEGER*4 LSEGMC,LSEGMF,LSEGMM,LSEGMM2,LSEGMN,LSEGMC2
      INTEGER*4 NSEGMC,NSEGMF,NSEGMM,NSEGMM2,NSEGMN,NSEGMC2
      INTEGER*4 NWMASTC,NWMASTF,NWMASTM,NWMASTM2,NWMASTN,NWMASTC2


      ! neut version
      INTEGER IDATAN(3)
      REAL    RDATAN(3)
      CHARACTER*32 CNAMEN
      EQUIVALENCE (IDATAN(1),RDATAN(1))

      ! flux and geometry
      INTEGER IDATAF(6)
      REAL    RDATAF(6)
      CHARACTER*32 CNAMEF
      EQUIVALENCE (IDATAF(1),RDATAF(1))

      ! nuclear effects
      INTEGER IDATAC(11)
      REAL    RDATAC(11)
      CHARACTER*32 CNAMEC
      EQUIVALENCE (IDATAC(1),RDATAC(1))

      ! reaction models
      INTEGER IDATAM(10)
      REAL    RDATAM(10)
      CHARACTER*32 CNAMEM
      EQUIVALENCE (IDATAM(1),RDATAM(1))

      ! reaction models
      INTEGER IDATAM2(17)
      REAL    RDATAM2(17)
      CHARACTER*32 CNAMEM2
      EQUIVALENCE (IDATAM(1),RDATAM(1))

      ! nuclear effects #2 
      INTEGER IDATAC2(12)
      REAL    RDATAC2(12)
      CHARACTER*32 CNAMEC2
      EQUIVALENCE (IDATAC2(1),RDATAC2(1))


C--- NEUT version
      IDATAN(1)=corev
      IDATAN(2)=nucev
      IDATAN(3)=nuccv

      CNAMEN='NEUTVER'

      NWMASTN=3
C      NSEGMN=corev
      NSEGMN=1
      LSEGMN=3

      CALL KZBCR1(CNAMEN,NWMASTN,'3I',NSEGMN,LSEGMN,'3I',IERR)      
      IF (IERR.NE.0) GOTO 1000
      CALL KZREP1(CNAMEN,0,NWMASTN,IDATAN(1))


C--- Flux and geometry
      IDATAF(1)=SK_GEOMETRY
      IDATAF(2)=0
      RDATAF(3)=0
      IDATAF(4)=0
      RDATAF(5)=0
      IDATAF(6)=ITAUFLGCORE

      CNAMEF='NEUTFLUX'

      NWMASTF=6
C      NSEGMF=GEOMVER
      NSEGMF=1
      LSEGMF=6

      CALL KZBCR1(CNAMEF,NWMASTF,'6I',NSEGMF,LSEGMF,'6I',IERR)
      IF (IERR.NE.0) GOTO 1000
      CALL KZREP1(CNAMEF,0,NWMASTF,IDATAF(1))
      CALL KZREP1(CNAMEF,1,NWMASTF,IDATAF(1))

C-- Obsolete : 2016/02/28 : Use NEUTINFO defined in nemknetarg.F 
C      
C
CC--- Configuration of nuclear efffects, etc...
C      IDATAC(1)=NEFRMFLG
C      IDATAC(2)=NEPAUFLG
C      IDATAC(3)=NENEFO16
C      IDATAC(4)=NEMODFLG
C      RDATAC(5)=PFSURF
C      RDATAC(6)=PFMAX
C      RDATAC(7)=VNUINI
C      RDATAC(8)=VNUFIN
C      IDATAC(9)=IFORMLEN
C      IDATAC(10)=NUCRESCAT
C      RDATAC(11)=XNUCFACT
C
C      CNAMEC='NEUTCONFIG'
C
C      NWMASTC=11
CC      NSEGMC=NEFRMFLG
C      NSEGMC=1
C      LSEGMC=11
C
C      CALL KZBCR1(CNAMEC,NWMASTC,'4I4F2I1F',NSEGMC,LSEGMC,
C     $            '4I4F2I1F',IERR)      
C      IF (IERR.NE.0) GOTO 1000
C      CALL KZREP1(CNAMEC,0,NWMASTC,IDATAC(1))
C      CALL KZREP1(CNAMEC,1,NWMASTC,IDATAC(1))
C
C
CC--- Reaction model
C      IDATAM(1)=MODELDIS
C      IDATAM(2)=MODELCOH
C
C      CNAMEM='NEUTMODEL'
C
C      NWMASTM=2
CC      NSEGMM=MODELDIS
C      NSEGMM=1
C      LSEGMM=7
C
C      CALL KZBCR1(CNAMEM,NWMASTM,'2I',NSEGMM,LSEGMM,'5I2F',IERR)      
C      IF (IERR.NE.0) GOTO 1000
C      CALL KZREP1(CNAMEM,0,NWMASTM,IDATAM(1))
C
C      IDATAM(1)=MDLQE
C      IDATAM(2)=MDLSPI
C      IDATAM(3)=MDLCOH
C      IDATAM(4)=MDLDIS
C      IDATAM(5)=MDLQEAF
C      RDATAM(6)=XMAQE
C      RDATAM(7)=XMASPI
C      CALL KZREP1(CNAMEM,1,LSEGMM,IDATAM(1))      
C
CC--- Reaction model #2
C      RDATAM2(1)=XMARES
C      RDATAM2(2)=XMVQE
C      RDATAM2(3)=XMVSPI
C      RDATAM2(4)=XMVRES
C      RDATAM2(5)=KAPP
C      RDATAM2(6)=XMACOH
C      RDATAM2(7)=RAD0NU
C      RDATAM2(8)=fA1COH
C      RDATAM2(9)=fb1COH
C      IDATAM2(10)=IFFSPI
C      IDATAM2(11)=NRTYPESPI
C      RDATAM2(12)=RCA5ISPI
C      RDATAM2(13)=RBGSCLSPI
C      RDATAM2(14)=SCCFV
C      RDATAM2(15)=SCCFA
C      RDATAM2(16)=FPQE
C      RDATAM2(17)=PFSF
C
C      CNAMEM2='NEUTMODEL2'
C
C      NWMASTM2=17
CC      NSEGMM=MODELDIS
C      NSEGMM2=1
C      LSEGMM2=17
C
C      CALL KZBCR1(CNAMEM2,NWMASTM2,'9F2I6F',
C     $            NSEGMM2,LSEGMM2,'9I2I6F',IERR)      
C      IF (IERR.NE.0) GOTO 1000
C      CALL KZREP1(CNAMEM2,0,NWMASTM2,IDATAM2(1))
C
CC--- Configuration of nuclear efffects weights
C      RDATAC2(1)=FEFQE
C      RDATAC2(2)=FEFQEH
C      RDATAC2(3)=FEFINEL
C      RDATAC2(4)=FEFABS
C      RDATAC2(5)=FEFCOH
C      RDATAC2(6)=FEFQEHF
C      RDATAC2(7)=FEFCOHF
C      RDATAC2(8)=FEFCX
C      RDATAC2(9)=FEFCXHF
C      RDATAC2(10)=FEFCXH
C      RDATAC2(11)=FEFCOUL
C      RDATAC2(12)=FEFALL
C
C      CNAMEC2='NUCEFPARAM'
C
C      NWMASTC2=12
CC      NSEGMC=NEFRMFLG
C      NSEGMC2=1
C      LSEGMC2=12
C
C      CALL KZBCR1(CNAMEC2,NWMASTC2,'12F',NSEGMC2,LSEGMC2,
C     $            '12F',IERR)      
C      IF (IERR.NE.0) GOTO 1000
C      CALL KZREP1(CNAMEC2,0,NWMASTC2,IDATAC2(1))
C
C

#ifdef VERBOSE
      print*, 'NEUTVER-> ',IDATAN
      print*, 'NEUTFLUX-> ',
     &     IDATAF(1),IDATAF(2),RDATAF(3),IDATAF(4),
     &     RDATAF(5),IDATAF(6)
      print*, 'NEUTCONFIG-> ',
     &     IDATAC(1),IDATAC(2),IDATAC(3),IDATAC(4),
     &     RDATAC(5),RDATAC(6),RDATAC(7),RDATAC(8),IDATAC(9),
     &     IDATAC(10),RDATAC(11)
      print*, 'NEUTMODEL-> ',IDATAM
#endif

      RETURN

C     
C     ++ ERROR MESSAGE
C     
 1000 WRITE(6,900)
  900 FORMAT(' *** ERROR *** ( BANK CREATION IN NEMKMODELBK )')
      RETURN
      END
