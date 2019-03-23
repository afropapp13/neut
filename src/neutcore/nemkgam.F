************************************************************************
*     ------------------------------------------------
      SUBROUTINE NEMKGAM
*     ------------------------------------------------
C     
*     (Purpose)
C     MAKE VECTORS OF NUCLEAR EXCITATION GAMMA
C     
C     (Input)
*     COMMON /VCWORK/
*     COMMON /VCVRTX/
C     
C     (Output)
*     COMMON /VCWORK/
*     COMMON /VCVRTX/
C     
*     (Creation Date and Author)
C     23-JUN-96    Y.Hayato (S.K. version)
C     Modified:
C     7-OCT-99     D.Casper  Use all de-excitation modes from Ejiri paper
C     2013.03.12   P.Sinclair Skip if event is 2p2h
C     
      implicit none

#include "vcvrtx.h"
#include "vcwork.h"
      
#include "posinnuc.h"

      INTEGER  NEMODGAM
      EXTERNAL NEMODGAM

      INTEGER  MODEGAM
      REAL     DUMMY

      integer*4 i,j

      real*4    vmod
      external  vmod
      

      integer max_emit, n_emit
      parameter (max_emit = 5)	
      logical pion_absorbed   
      integer id_emit(max_emit)
      real p_emit(4,max_emit)
      integer proton_code, neutron_code, gamma_code
      parameter (proton_code = 2212)
      parameter (neutron_code = 2112)
      parameter (gamma_code = 22)
C     
C     + Ignore if target isn't neutron or proton
C     
      if (ipvc(2).ne.proton_code .and. ipvc(2).ne.neutron_code) return
  
C     Skip if mode is MEC
      if (((ipvc(2).eq.proton_code).or.(ipvc(2).eq.neutron_code)).and.
     $    ((ipvc(3).eq.proton_code).or.(ipvc(3).eq.neutron_code)).and.
     $    ((ipvc(5).eq.proton_code).or.(ipvc(5).eq.neutron_code)).and.
     $    ((ipvc(6).eq.proton_code).or.(ipvc(6).eq.neutron_code))) then
        return
      endif

      if (ibound.eq.0) then
         return
      endif

      do 100 i=1,NVC
         if (icrnvc(i).eq.1) goto 110
 100  continue
C      -- there is no active particle : necessary to return
      return
 110  continue

C     
C     + Determine if target is free proton, or whether any pions were absorbed
C     
      MODEGAM = NEMODGAM(DUMMY)
      if (modegam.le.0 .or. modegam.gt.2) then
         return
      else if (modegam.eq.1) then
         pion_absorbed = .false.
      else if (modegam.eq.2) then
         pion_absorbed = .true.
      endif
C     
C     + Generate de-excitation products
C     
      call nudeex(pion_absorbed,ipvc(2),0,n_emit, id_emit, p_emit)
C     
C     + Avoid particle stack overflow
C     
      if (n_emit + nvc .gt. maxvc) then
         write (*,*) 'nemkgam: de-excitation products truncated'
         n_emit = max(0,maxvc-nvc)
      endif
C     
C     + Store de-excitation products
C     
      do i = 1, n_emit
C     type *,'nemkgam: ',id_emit(i),p_emit(4,i)
         ipvc(nvc+i) = id_emit(i)
         if (id_emit(i).eq.gamma_code)then
            amasvc(nvc+i) = 0.
         else
            amasvc(nvc+i) 
     $           = sqrt(p_emit(4,i)**2 - vmod(p_emit(1,i),3)**2)
         endif
         call ucopy(p_emit(1,i),pvc(1,nvc+i),3)
         call ucopy(posvc,posivc(1,nvc+i),3)
         call ucopy(posvc,posfvc(1,nvc+i),3)
         iorgvc(nvc+i) = 2      ! parent = initial nucleon  TO DO: FIX THIS 
         iflgvc(nvc+i) = 0      ! flag = "consider later procedure"
         icrnvc(nvc+i) = 1      ! flag = "chase"
         timvc(nvc+i) = 0.      ! starting time
         ivtivc(nvc+i) = 1      ! initial vertex=1
         ivtfvc(nvc+i) = 1      ! final vertex=1
         CALL MCMASS(IPVC(NVC+i),AMASVC(NVC+i))
         do 200 j=1,3
            POSNUC(j,nvc+i) = POSNUC(j,1)
 200     continue

      enddo	
C     
C     + update particle counter
C     
      nvc = nvc + n_emit
      
      RETURN
      END
