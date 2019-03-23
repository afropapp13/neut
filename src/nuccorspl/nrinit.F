      subroutine nrinit(anuc,icont)
*                                                                                                                                         
*     ( Creation Date and Author )
*       2007.11.05 ; G.Mitsuka - add support other target than 16O
*                                other target is treated by scaling
*                                to 16O
*
	  implicit none
	  real*4 anuc

#include "nrdis.h"
C     common /nucdis/ probmi(31),probdi(31),rhotab(30),pnorm
#include "nrcms.h"
C	  common /nuccms/rhon,pfermi,unucl,ecms2,up4cms(3),beta(3),prot
C	  logical prot

#include "nrnuclparam.h"

	  real*4    probd(30)
C	  dimension probm(30),probd(30)

C     DOVER ET AL	DISTRIBUTION
	  DATA PROBD/
     $         0.0,0.2,4.7,20.2,67.6,186.2,384.2,802.4,1409.0,2285.8,
     $         3264.0,3913.1,4665.6,4901.0,4170.9,4050.0,2973.6,1976.8,
     $         1296.0,985.5,704.0,370.4,232.3,177.7,80.6,35.7,20.0,10.0,
     $         5.0,2.5/


	  real*4    probs,rbin,rmin,rmax
	  integer*4 nbin,i,icont
      real*4 probm(70), RMSRADOXY
      parameter (RMSRADOXY = 2.730)


          

C SET UP MATTER DISTRIBUTION  IN NUCLEUS
C     UNITS ARE IN FERMIS

	  RBIN = 0.2
	  RMIN = 0.0
          RMAX = 6.*(NRRMSRAD/RMSRADOXY)

C	  Asmita R.: pfermi is now calculated in nrprton
C         pfermi=NRPFSURF*1000. ! NRPFSURF is GeV uniti

	  PROBS=0.0
C     e.g. NBIN=30(O),27(C),38(AR)
	  NBIN=NINT((RMAX-RMIN)/RBIN)
          print*, '#of bins in nucleus for nuccorspl = ',NBIN
	  call nrrodis(rmin,rbin,nbin,probm,rhotab)
          if(icont.eq.0) then
             print*, 'THIS MODE IS NOT GUARANTEED NOW (see nrinit.F)'
             stop
C             CALL nrCFD(PROBMI,NBIN,PROBM)
C             CALL nrCFD(PROBDI,NBIN,PROBD)
          endif
          
	  DO I=1,NBIN
C             print*, 'I:PROBM:RHOTAB',I,PROBM(I),RHOTAB(I)
             PROBS=PROBS+PROBM(I)
	  END DO
	  PNORM=ANUC/PROBS
	  DO I=1,NBIN
             RHOTAB(I)=PNORM*RHOTAB(I)
	  END DO

	  return
	  end
