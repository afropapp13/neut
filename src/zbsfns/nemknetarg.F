
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
      PARAMETER (NSEGM=8)

      integer*4 LSEGMTBL(10)
      DATA LSEGMTBL/ 3, 6,11, 7, 5,27,27,10, 0, 0/

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
C     +10 : Length of 8th segmetn  : Total cross-section
C

      ISEGM   = 0
      NDATA   = 10
      IDATA(1)= NDATA
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
C     + 5 : QE potential model
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
C     +27 : 

      ISEGM    = 6
      NDATA    = LSEGMTBL(ISEGM)

      DO 30 I=1,27
         RDATA(I)=CRSNEUT(I)
 30   continue

      call kzrep0(CNAME,ISEGM,'F',NDATA,IDATA)

C--------------------------      
C- Segment 7
C       Cross-section weights ( anti neutrino )
C     + 1 : 
C      |  : Cross-section weight for anti neutrino
C     +27 : 
C

      ISEGM    = 7
      NDATA    = LSEGMTBL(ISEGM)

      DO 40 I=1,27
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

 999  continue
      if (ierr.ne.0) then
         write(*,*) 'Failed to create BANK : ',CNAME
      endif

      return

      end
