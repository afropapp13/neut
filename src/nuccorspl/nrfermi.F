      subroutine nrfermi(r)
*
*     ( Creation Date and Author )
*       2007.11.05 ; G.Mitsuka - add support other target than 16O
*                                other target is treated by scaling
*                                to 16O
*
      implicit none
      real*4 r

#include "nrdis.h"
c     ommon /nrdis/probmi(31),probdi(31),rhotab(30),pnorm
#include "nrint.h"
C      common /nrint/uin,chint,ein,pin(3),upin(3)
#include "nrcms.h"
C      common /nrcms/ rhon,pfermi,unucl,ecms2,up4cms(3),beta(3),prot
C      logical prot
#include "nrprat.h"
C      common /nrpart/prat
#include "nrnuclparam.h"
#include "neutparams.h" 
#include "nework.h"
C      dimension rec(3),p4cms(3)
      real*4 rec(3),p4cms(3)

      real*4    toten2
      parameter (toten2=931.**2)

      real*4    rbin,rmin,rmax
      integer*4 nbin

      real*4   ranf
      external ranf

      integer*4 irlo,irhi
      real*4    rrem,rhonf,arec,e4cms,p4cmsm
      integer*4 j

      real*4 RMSRADOXY
      parameter (RMSRADOXY = 2.730)

       REAL effrmgas
       


C     SET PARAMETER FOR EACH TARGET NUCLEI
      prat = AP/APN
      if(prat.gt.1.) then
         print*, 'NUMBNDP is larger than NUMATOM, see necard.h'
         stop
      endif
      rbin = 0.2
      rmin = 0.
      rmax = 6.*(NRRMSRAD/RMSRADOXY)
      nbin=nint((rmax-rmin)/rbin)

C     Asmita R.:For MEC events fermi surface momentum is function of radial distance r (22-07-2013)
C     while for other modes its read from card
      pfermi=NRPFSURF*1000. ! NRPFSURF is GeV unit
        if(abs(MODENE).EQ.2) then 
        	pfermi=EFFRMGAS(0,0,r)
        end if
      
     
c     decide if interaction is on proton (prot=.true.) or neutron (.false.)
      if(ranf().lt.prat)then
         prot=.true.
      else
         prot=.false.
      end if
C     NUCLEAR DENSITY AT THIS VALUE OF R
C     INCLUDE 0.1   FACTOR is TO CHANGE FROM MB TO FMSQ (cms is in fm!)
      call nrhis(r-rmin,rbin,nbin,irlo,irhi,rrem)
      RHON=(RHOTAB(IRHI)*RREM+RHOTAB(IRLO)*(1.0-RREM))*0.1
      RHONF=RHON/PNORM

C     FERMI MOMENTUM AT THIS VALUE OF R
C     PFERMI=250.0*((RHONF*10.)**0.3333)
C     IF(PFERMI.LE.0.0) PFERMI=1.0

      AREC=pfermi*RANF()**0.3333
      call nrranve(arec,rec,2.,-1.)
      unucl=sqrt(toten2-pfermi**2)

C     ENERGY AND 3 MOMENTUM PARTS   OF CMS 4 MOMENTUM
      E4CMS=EIN+sqrt(unucl**2+arec**2)
      DO J=1,3
         P4CMS(J)=PIN(J)+REC(J)
      END DO
      P4CMSM=SQRT(P4CMS(1)**2+P4CMS(2)**2+P4CMS(3)**2)
      DO J=1,3
         UP4CMS(J)=P4CMS(J)/P4CMSM
         BETA(J)= -P4CMS(J)/E4CMS
      END DO
C     SQUARE OF CMS ENERGY
      ECMS2=E4CMS**2-p4cmsm**2

      return
      end
