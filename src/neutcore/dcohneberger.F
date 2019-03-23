************************************************************************
*     ----------------------------
      FUNCTION DCOHNEBERGER(E,IPAR,IMODE,X,Y,Z,PHI)
*     ----------------------------
*
*     (Purpose)
*       RETURN DIFFERENTIAL CROSSSECTION WHICH DEPENDS X,Y,Z,PHI.
*
*
*     (Input)
*       E     : INCIDENT NEUTRINO ENERGY ( GEV )
*       IPAR  : Neutrino particle type
*       IMODE : Interaction mode ( +-16 or +-36 )
*       COMMON /nework/
*
*     (Output)
*       DCOHNE : DIFFERENTIAL CROSSECTION
*     
*     (Creation Date and Author)
*       1994.08.30 ; G.DION & K.Okumura
*       1996.02.05 ; K.Kaneyuki   DCOH -> DCOHNE
*       2007.01.08 ; G.Mitsuka,  Add lepton mass correction,
*                                rename from DCOHNE to DCOHNEREIN
*       2007.11.05 ; G.Mitsuka,  support target nucleus besides 16O
*                                Atomic number is given by "NUMATOM"
*
*     (Comment)
*       THIS FUNCTION NEED FOR SUBROUTINE 'COHVCT' WHICH IS A MONTE-CALRO
*       SUMILATION PROGRAM OF COHERENT PI0 PRODUCTION.
*
************************************************************************
      IMPLICIT NONE

#include "necard.h"
#include "nework.h"
#include "neutmodel.h"
#include "pdf804/w50516.inc"
C     COMMON /NVECT/ IMOD,LVECT,IP(20),PIN(3,20)
C

      REAL*4    DCOHNEBERGER
      INTEGER*4 IPAR,IMODE

      REAL*4    E,X,Y,Z,PHI

      real XLEPMASS

      REAL*4 TINY
      PARAMETER (     TINY    = 1.E-4               )

      REAL*4 PI
      PARAMETER (     PI      = 3.141592653589793D0 )

      REAL*4 HBARC
      PARAMETER (     HBARC   = .197327053          )

      REAL*4 GFERMI
      PARAMETER (     GFERMI  = 0.89719167E-7       )

      REAL*4 FPIFAC
      PARAMETER (     FPIFAC  = 0.93                )
 
C     Now included in neutmodel.h
C      REAL*4 AM
C      PARAMETER (     AM      = 1.                  )
C
C      REAL*4 RAD0NU
C      PARAMETER (     RAD0NU  = 1.0                 )

      REAL*4 BRNTC2
      PARAMETER (     BRNTC2  = 1.E-24              )

      REAL*4 BRNTF2
      PARAMETER (     BRNTF2  = 1.E+02              )

      REAL*4 FMTCTM
      PARAMETER (     FMTCTM  = 1.E-13              )
      
      REAL*4 DUNIT
      PARAMETER (     DUNIT   = 1.E+12              )
C
      REAL*4 KIN, ARG1, ARG2, u,v,q,cosa, Jacobian, t
      REAL*4 FRONT, THETACAB

      INTEGER*4 ITYPE
      common/BERGER/itype
      
      INTEGER*4 IPAROUT

      REAL*4 PRMASS,PIMASS

      REAL*4 ATOMN,CORR

      REAL*4 QSMX,QSMN,Q2min,Q2max,q2
      REAL*4 XMAX,XMIN,YMAX,YMIN,Y2MPE2
      REAL*4 EPITOT,EPIKIN

      REAL*4 DSIGEL

C
C     PION AND PROTON MASS
C
      CALL MCMASSGV(IPNE(5),PIMASS)
      CALL MCMASSGV(2112,PRMASS)

      if (imode.ne.modene) then
         write(*,*) 
     $        'dcohneberger: Interaction mode flags are inconsistent'
         write(*,*) 'Imode = ',imode, '/ modene=', modene
         stop
      endif

      if ((abs(imode).ne.16).and.(abs(imode).ne.36)) then
         write(*,*) 
     $        'dcohneberger: Incorrect interaction mode?.'
         write(*,*) 'Imode = ',imode
         stop
      endif
C
C     SET LEPTON MASS ACCORDING TO THE INTERACTION MODE
C
      if (abs(imode).eq.36) then
         XLEPMASS = 0
      else
         IPAROUT = abs(IPAR) - 1
         CALL MCMASSGV(IPAROUT,XLEPMASS)
      endif

      ATOMN = float(NUMATOM)
      DCOHNEBERGER = 0.
      CORR = 1.

      if(Y.lt.PIMASS/E) goto 100
      if(Y.gt.1./(1.+PRMASS*X/(2.*E))) goto 100

      QSMX=4.*E*E*(1.-Y)
      QSMN=0.
      if(QSMX.lt.0.) goto 100

      XMAX=QSMX/(2.*PRMASS*Y*E)
      XMIN=QSMN/(2.*PRMASS*Y*E)
      if(X.gt.XMAX) goto 100
      if(X.lt.XMIN) goto 100

      if(abs(modene).eq.36) goto 101

C     **************************
C      ! for CC cross section(G.Mitsuka, Jan 2006)
C     PCAC CORRECTION FACTOR 

      ymin = PIMASS/E
      ymax = 1.-XLEPMASS/E

      Q2min = XLEPMASS**2.*Y/(1-Y)
      Q2max = 2.*PRMASS*E*(1.-XLEPMASS/E)
      Q2 = 2.*PRMASS*E*X*Y
      THETACAB = 13.02*PI/180

      corr = 1.-Q2min/(2.*(Q2+PIMASS**2.))
      corr = corr*corr + Y/4.*Q2min*(Q2-Q2min)/(Q2+PIMASS**2.)**2.   
C      corr = corr*2*COS(THETACAB)**2
C
C--- cross-section became half for CC but to keep consistency 
C    with the previously generated events
C    also, dcohmxne uses 1/2 values for CC.
C      
      corr = corr*COS(THETACAB)**2
      if(Y.lt.ymin.or.Y.gt.ymax.or.
     &     Q2.lt.Q2min.or.Q2.gt.Q2max) corr=0.
C     **************************
 101  continue 
C
C  ((OUTGOING PION MOMENTUM) / (INCIDENT NEUTRINO ENERGY))**2
      Y2MPE2=Y**2-(PIMASS/E)**2

C  CONSTANT FACTOR
      FRONT = (GFERMI*PIMASS)**2
      FRONT = FRONT*FPIFAC**2/(4*PI**2)
      FRONT = FRONT/HBARC**6

C  Convert in unit of 10^-38 cm^2  
      FRONT = FRONT*1.E-24*1.E+38

C  COS ALPHA WITH LEPTON MASS
      cosa = (2*E**2*(1-Y)-2*PRMASS*E*X*Y-XLEPMASS**2)/
     & (2*E*SQRT(ABS(E**2*(1-Y)**2-XLEPMASS**2))) 
      IF(ABS(cosa).GT.1) goto 100

C  Q WITH LEPTON MASS
      q = SQRT(ABS(2*E**2*(1-Y)*(1-SQRT(ABS(1-XLEPMASS**2/
     & (E**2*(1-Y)**2)))*cosa)+(E*Y)**2)-XLEPMASS**2)  
      IF(ISNAN(q).OR.q.EQ.0) goto 100    

C  u AND v AS DESCRIBED IN B&S PAPER
      u = (E*(2-Y)+q)/(2*E)
      v = (E*(2-Y)-q)/(2*E)
      Jacobian = 2*E**2*Y**2*SQRT(1-(PIMASS/(E*Y))**2)*
     & SQRT(1+2*PRMASS*X/(E*Y)) * 2*PRMASS*E*Y

C  FULL KINEMATIC FACTOR
      KIN = E*u*v*Jacobian/q
     & *(XMACOH**2/(2.*PRMASS*E*X*Y+XMACOH**2))**2

C  PION KINEMATIC NEERGY
      EPIKIN = E*Y-PIMASS
      IF(EPIKIN.LE.0) goto 100

C  CALCULATE |t| WITH A SCALING
      ARG1=Y+PRMASS*X/E-PIMASS**2/(2*E**2*Y)
      ARG1=Y*(ARG1-Z*(1+PRMASS*X/E)*SQRT(ABS(Y2MPE2)))

      ARG2=(1.-Z**2)*Y2MPE2*(2.*PRMASS/E)*X*Y
      ARG2=ARG2*(1.-Y*(1.+PRMASS*X/(2.*E)))
      ARG2=SQRT(ABS(ARG2))*COS(PHI)
      t = 2*E**2*(ARG1+ARG2)*(ATOMN/12.)**(2/3)
      
      IF(ISNAN(t)) goto 100  

C  PION-NEUCLEON INELASTIC CROSS SECTION SUBROUTINE (UNIT FM)
      CALL NECROSEL(EPIKIN,t,DSIGEL)     
      IF (DSIGEL.LT.0.00) goto 100

C  DIFFERENTIAL XS WITH A SCALING
      DCOHNEBERGER = FRONT*KIN*DSIGEL*corr*(ATOMN/12.)**2
    
C      write(*,*) 'DCOHNER = ',DCOHNEBERGER
C      write(*,*) '-------------------------------------------'

 100  RETURN
      END
