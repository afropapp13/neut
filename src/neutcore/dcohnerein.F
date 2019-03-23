************************************************************************
*     ----------------------------
      FUNCTION DCOHNEREIN(E,IPAR,IMODE,X,Y,Z,PHI)
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

      REAL*4    DCOHNEREIN
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
      REAL*8 ARG, SOLN, ARG1, ARG2, ARG3
      REAL*8 FRONT

      INTEGER*4 ITYPE
      common/REIN/itype

      INTEGER*4 IPAROUT

      REAL*4 PRMASS,PIMASS

      REAL*4 ATOMN,CORR

      REAL*4 QSMX,QSMN,Q2min,Q2max,q2
      REAL*4 XMAX,XMIN,YMAX,YMIN,Y2MPE2
      REAL*4 EPITOT,EPIKIN

      REAL*4 BTE,P2OS

      REAL*4 FEMTO,R
      REAL*4 DSIGDT

C
C     PION MASS
C
      CALL MCMASSGV(IPNE(5),PIMASS)
      CALL MCMASSGV(2112,PRMASS)

      if (imode.ne.modene) then
         write(*,*) 
     $        'dcohnerein: Interaction mode flags are inconsistent'
         write(*,*) 'Imode = ',imode, '/ modene=', modene
         stop
      endif

      if ((abs(imode).ne.16).and.(abs(imode).ne.36)) then
         write(*,*) 
     $        'dcohnerein: Incorrect interaction mode?.'
         write(*,*) 'Imode = ',imode
         stop
      endif

      if (abs(imode).eq.36) then
         XLEPMASS = 0
      else
         IPAROUT = abs(IPAR) - 1
         CALL MCMASSGV(IPAROUT,XLEPMASS)
      endif

      ATOMN = float(NUMATOM)
      DCOHNEREIN = 0.
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
      ! for CC cross section(G.Mitsuka, Jan 2006)
      ymin = PIMASS/E
      ymax = 1.-XLEPMASS/E

      Q2min = XLEPMASS**2.*Y/(1-Y)
      Q2max = 2.*PRMASS*E*(1.-XLEPMASS/E)
      Q2 = 2.*PRMASS*E*X*Y

      corr = 1.-Q2min/(2.*(Q2+PIMASS**2.))
      corr = corr*corr + Y/4.*Q2min*(Q2-Q2min)/(Q2+PIMASS**2.)**2.      
      if(Y.lt.ymin.or.Y.gt.ymax.or.
     &     Q2.lt.Q2min.or.Q2.gt.Q2max) corr=0.
C     **************************

 101  continue

C     
C
C  ((OUTGOING PION MOMENTUM) / (INCIDENT NEUTRINO ENERGY))**2
      Y2MPE2=Y**2-(PIMASS/E)**2
C
C  BTE : (RADIUS OF NEUCLEI)*E**2
      BTE   = RAD0NU**2*ATOMN**(2./3.)*E**2/(3.*HBARC**2)
C
      FRONT = PRMASS*(ATOMN*PIMASS*GFERMI)**2
      FRONT = FRONT*E**3
      FRONT = FRONT*FPIFAC**2/PI**2
      FRONT = FRONT/HBARC**8
      FRONT = FRONT*EXP(BTE*(PIMASS/E)**2)
C      FRONT = FRONT*DFLOAT(FMTCTM)**2
C
C     UNIT : 10**-38 cm**2 = 10**-12 fm **2 
C      
      FRONT = FRONT*DUNIT
C      write(6,*) ' dcohne : front = ',front,fmtctm
C     
      SOLN = SQRT(ABS(Y*(Y+2.*PRMASS*X/E)*ABS(Y2MPE2)))
      SOLN = SOLN*(1.-Y)/(1.+2.*PRMASS*E*X*Y/XMACOH**2)**2
C
C  PION TOTAL ENERGY
      EPITOT = E*Y
C
C  PION KINEMATIC NEERGY
      EPIKIN = E*Y-PIMASS
C
      ARG1=Y+PRMASS*X/E
      ARG1=Y*(ARG1-Z*(1+PRMASS*X/E)*SQRT(ABS(Y2MPE2)))
      ARG1=-2.*BTE*ARG1
C
      ARG2=(1.-Z**2)*Y2MPE2*(2.*PRMASS/E)*X*Y
      ARG2=ARG2*(1.-Y*(1.+PRMASS*X/(2.*E)))
      ARG2=-2.*BTE*SQRT(ABS(ARG2))*COS(PHI)
C
C  PION-NEUCLEON INELASTIC CROSS SECTION SUBROUTINE (UNIT FM)
      CALL NECROSIN(EPIKIN, FEMTO)
C      write(6,*) ' dcohne : FEMTO = ',EPIKIN,FEMTO
C
C  EXP(ARG3) = FABS : PION ABSOLUTION COEFFICIENT
      ARG3=-9.*ATOMN**(1./3.)*FEMTO/(16.*PI*RAD0NU**2)
C
      ARG = DBLE(ARG1+ARG2+ARG3)
C      write(6,*) ' dcohne : arg = ',arg,arg1,arg2,arg3
C     
      IF (ARG.LT.-20.D0) THEN
         DCOHNEREIN = 0.0
         RETURN
      ELSE
         SOLN = SOLN*DEXP(ARG)
      ENDIF
C
C      P2OS = ABS(EPITOT**2 - PIMASS**2)
C      P2OS = P2OS/(EPITOT + ((PRMASS**2+PIMASS**2)/(2.*PRMASS)))**2
C      P2OS = P2OS/4.
C
      P2OS = 1./4.
C
C  PION-NEUCLEON TOTAL CROSS SECTION
      CALL NECROSTO(EPIKIN,FEMTO)
C      write(6,*) ' dcohne : femto = ',femto,epikin
C     
C  CASE REGGE FORM OF TOTAL CROSS SECTION
C      FEMTO = 2.4 + 1.2/SQRT(EPIKIN+TINY)
C
C  R : REAL(SCATTERING AMPLITUDE) / IM(SCATTERING AMPLITUDE)
      CALL NEGETR(EPIKIN,EPITOT,FEMTO,R)
C      write(6,*) ' dcohne : r = ',r,epikin,epitot,femto
C
      DSIGDT = P2OS*FEMTO**2*(1.+R**2)/(4.*PI)
C
C      write(6,*) ' dcohne : dcohne = ',front,soln,dsigdt,
C     &            front*soln*dsigdt
      DCOHNEREIN = FRONT*SOLN*DSIGDT*corr
C
 100  RETURN
      END
