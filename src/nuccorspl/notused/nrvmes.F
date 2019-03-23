	  subroutine nrvmes(jpar,ides)
	  implicit none
	  integer*4 jpar,ides

#include "nrint.h"
C	  common /nucint/uin,chint,ein,pin(3),upin(3)
#include "nrcms.h"
C	  common /nuccms/ rhon,pfermi,unucl,ecms2,up4cms(3),beta(3),prot
#include "nrhad.h"
C	  common /nuchad/ihap,nel(3)
#include "nrhadr.h"
C	  common /hadron/ inot,itrans,idec,ine1,ine2,iel,iabs,icx
C	  logical prot

	  real*4    xeta(20),xome(20)
C	  dimension xeta(20),xome(20)
C     PARTIAL AND	TOTAL CROSS SECTIONS
C     PI-P TO ETA+N CROSS	SECTION	(MB) AT     20*20MEV/C VALUES
C     OF THE CMS MOM. OF THE ETA
      DATA XETA/
	1	   0.,.1,.25,.55,1.,1.25,1.4,1.6,1.7,1.75,
	2	   1.8,1.8,1.75,1.7,1.6,1.4,1.3,1.15,1.0,.9/
C     PI-P TO OME+N CROSS	SECTION	(MB) AT     20*20MEV/C VALUES
C     OF THE CMS MOM. OF THE OME
      DATA XOME/
	1	   0.,.1,.15,.28,.45,.65,.90,1.05,1.20,1.35,
	2	   1.5,1.65,1.8,1.95,2.1,2.25,2.4,2.55,2.7,2.85/


	  real*4 pout1l(3),poutel(3),rv(3)
C	  dimension pout1l(3),poutel(3),rv(3)

	  real*4 upi,step
	  parameter (upi=139.6,  step=0.2)

	  real*4 pcmsm2,pcmsm1,pinm,prem,pint,pdec,pela,ptot
	  real*4 finch,xr,yr,dist,fsum,poutnlm,poutelm
	  real*4 pout1lm,pout2lm,eoutel,eout1l
	  integer*4 ilo,ihi,inel,j

	  real*4   ranf
	  external ranf

	  ides=1
	  IF(ECMS2.LE.(UPI+UNUCL)**2) GO TO 200
	  IF(ECMS2.LE.(uin+unucl)**2) GO TO 200

	  PCMSM2=SQRT((ECMS2-(UPI-UNUCL)**2)*(ECMS2-(UPI+UNUCL)**2)/(4.0*ECMS2))
	  PCMSM1=SQRT((ECMS2-(UIn-UNUCL)**2)*(ECMS2-(UIn+UNUCL)**2)/(4.0*ECMS2))
	  pinm=sqrt(pin(1)**2+pin(2)**2+pin(3)**2)

	  call nrhis(pcmsm1,20.,20,ilo,ihi,prem)
	  if(jpar.eq.1)then
		 PINT=1.5/.71*
	1		  (XETA(IHI)*PREM+XETA(ILO)*(1.-PREM))*RHON*(PCMSM2/PCMSM1)**2
		 PDEC=UIn/(PINM*230000.0)
		 if(pcmsm1.lt.100.)then
			pela=.43*pcmsm1
		 else if(pcmsm1.gt.270.)then
			pela=156.-.42*pcmsm1
			if(pela.lt.0.)pela=0.
		 else
			pela=43.
		 end if
		 pela=pela*rhon
		 PTOT=PDEC+PINT+pela
	  else
		 pint=0.5*
	1		  (xome(ihi)*prem+xome(ilo)*(1.-prem))*RHON*(PCMSM2/PCMSM1)**2
		 pdec=uin/(20.*pinm)
		 pela=0
		 ptot=pdec+pint
	  end if

	  FINCH=0.6666*pint
C     RHON FACTOR	ALREADY	IN

	  XR=ranf()+1.e-6
	  YR=ranf()*ptot
	  DIST=-LOG(XR)/ptot

C**********IF NOTHING HAPPENS
	  IF(STEP.LT.DIST)GO TO 200
C**********IF SOMETHING HAPPENS
c     
C**********DOES IT DECAY?
	  FSUM=pDEC
	  IF(YR.LT.FSUM)THEN
		 ihap=1
		 IDEC=IDEC+1
		 GO TO 99
	  END IF


C     VECTOR MESON MUST HAVE INTERACTED

	  call nrranve(1.,rv,2.,-1.)

C****************DOES VECTOR	MESON ELASTICALLY SCATTER?
	  FSUM=FSUM+pELA
	  IF(YR.LT.FSUM)THEN

		 call nrkin(rv,pcmsm1,uin,poutnlm,poutel,poutelm,eoutel)
		 IF(POUTNLM.LT.PFERMI) go to 200
         DO J=1,3
			PIN(J)=POUTEL(J)
			UPIN(J)=PIN(J)/POUTELM
         END DO
         EIN=EOUTEL
		 ihap=1
		 nel(1)=1
         GO TO 200
	  END IF

C****************VECTOR MESON MAKES A PION*********

	  call nrkin(rv,pcmsm2,upi,pout2lm,pout1l,pout1lm,eout1l)
	  IF(POUT2LM.LT.PFERMI)GO TO 200
	  ihap=1
	  INEL=INEL+1
C*****IS SCATTERING OFF P OR N?
	  FSUM=FSUM+FINCH
	  if(prot)then
C*****PROTON
		 IF(YR.LT.FSUM)THEN
			CHINT=1.0
		 ELSE
			CHINT=0.0
		 END IF
	  ELSE
C*****NEUTRON
		 IF(YR.LT.FSUM) THEN
			CHINT=-1.0
		 ELSE
			CHINT=0.0
		 END IF
	  END IF
	  DO J=1,3
		 PIN(J)=POUT1L(J)
		 UPIN(J)=PIN(J)/POUT1LM
	  END DO
	  EIN=EOUT1L
	  UIn=upi
 99	  ides=0
 200  return
	  end
