************************************************************************
*     ---------------------------------------------------
      DOUBLE PRECISION FUNCTION DNDXDY(IP,it,itype,EE,XX,YY)
*     ---------------------------------------------------
*
*     (Purpose)
*     ++ FUNCTION TO GIVE D(SIGMA)/DX/DY   -- BJORKEN SCALING
*                           ( 10**-38 cm2/GeV )
*
*     (Input)
*       IP     : NEUTRINO SPECIES
*       IT     : Nucleon SPECIES
*            =2212 proton
*            =2112 neutron
*       EE     : NEUTRINO ENERGY ( GEV )
*       ITYPE  : INTERACTION TYPE
*            =1 CHARGED CURRENT
*            =0 NEUTRAL CURRENT
*     
*     (Output)
*       XX     : Q**2/(2Mv) BJORKEN x
*       YY     : v/E        BJORKEN y
*
*     (Creation Date and Author)
*     
************************************************************************
      implicit none
#include "necard.h"
      integer ip,it,itype
      REAL*8 EE,XX,YY
      real*8 hc,gf,pi
      real*8 xmp,xmn
      real*8 consp,consn,u,d
      PARAMETER (HC    = 0.19732858D-13)
      PARAMETER (GF    = 1.16637D-5*HC)
      PARAMETER (PI    = 3.141592653589793238D0)
      PARAMETER (xmp    = 0.93827231D0)
      PARAMETER (xmn    = 0.93956563D0)
      parameter (consp =  1.D+38*GF*GF*xmp/PI)
      parameter (consn =  1.D+38*GF*GF*xmn/PI)

      real*8 const
      real*8 f2,xf3
      real*8 amw2
      parameter( amw2=80.2D0*80.2D0 )

      real*8 am,amlep
      real*8 cor1,cor2,corw
      real*8 COEFF2, COEFF2_1, COEFF2_2, R, fnrworld

      real*8 q2, Qtmp
      real*8 e,x,y
      real*8 a,b
      real*8 ymax,xmax
      real*8 rtt,qsmx
      real cost
      real plep,elep

      E=EE
      X=XX
      Y=YY

      DNDXDY=0.D0

      if(itype.eq.0) then
        amlep=0.D0
      elseif(iabs(ip).eq.12) then
        amlep=0.51099906D-3
      elseif(iabs(ip).eq.14) then
        amlep=0.105658389D0
      elseif(iabs(ip).eq.16) then
        amlep=1.7771D0
      else
        return
      endif

      if(it.eq.2212) then
        am = xmp
        const = consp
      elseif(it.eq.2112) then
        am = xmn
        const = consn
      else
        return
      endif

      YMAX=1.-AMLEP/E
      IF(YMAX.LT.0.) THEN
*      PRINT *,' ***** 1-ML/E < 0.******',X,Y
      RETURN
      END IF

      IF(Y.GT.YMAX) THEN
*      PRINT *,' ***** Y > 1-ML/E ****** X,Y,YMAX=',X,Y,YMAX
      RETURN
      END IF

*     momentum of outgoing lepton
      RTT=E*E*(1.-Y)**2-AMLEP**2
      IF(RTT.LT.0.) THEN
*      PRINT *,'*** RTT<0 X,Y,YMAX,RTT=',X,Y,YMAX,RTT
      RTT=0.
      END IF

*     q2 maximum
      QSMX=2.*E*E*(1.-Y)+2.*E*SQRT(RTT)-AMLEP**2
      IF(QSMX.LT.0.) THEN
*      PRINT *,' ***** Q2(X,Y,E,ML) < 0 ********',X,Y
      RETURN
      END IF

*     X maximum
      XMAX=QSMX/(AM**2+2.*AM*Y*E)
      IF(X.GE.XMAX) THEN
*      PRINT *,' ***** X > Q2MAX/(M**2+2MYE) *****',X,Y
      RETURN
      END IF

      XMAX=2.*Y*E/(2.*Y*E+AM)
      IF(X.GE.XMAX) THEN
*      PRINT *,' ****** W < 0.938GEV ******',X,Y
      RETURN
      END IF

      Q2=X*(2.*AM*Y*E+AM**2)
      ELEP=E-Y*E
      PLEP=ELEP**2-AMLEP**2

      IF(PLEP.LT.0.) THEN
*      PRINT *,' ***** P(E,Y,ML) < 0 ********',X,Y
      RETURN
      END IF

      PLEP=SQRT(PLEP)
      COST=(2.*E*ELEP-AMLEP*AMLEP-Q2)/2./E/PLEP

      IF(ABS(COST).GT.1.) THEN
*      PRINT *,' ****** COS(NEU,LEP) > 1 ******',COST,X,Y,AMLEP
      RETURN
      END IF

      Q2=2.D0*am*e*y*x
      if(nebodek.eq.0) then
!     --> Original
         COR1=(AMLEP**2)*Y/(4.*am*E*(X+1.E-30))-am*X*Y/(2.*E)
         COR1=COR1-(AMLEP**2)/(4.*E*E)-(AMLEP**2)/(2.*am*E*(X+1.E-30))
         COR2=-(AMLEP**2)/(4.*am*E*(X+1.E-30))
         
         corw=(AMW2/(AMW2+Q2))**2

         call structg(ip,it,e,x,y,f2,xf3,u,d)

         A=(1.-Y+Y**2/2.+COR1)*f2
         B=Y*(1-Y/2.+COR2)*xf3
      elseif(nebodek.eq.1) then
!     --> Bodek correction(Aug'06 G.Mitsuka)
         if(Q2.ge.0.35) then
            R = fnrworld(x,Q2)
         else if(Q2.lt.0.35) then
            Qtmp = 0.35
            R = 3.207*Q2*fnrworld(x,Qtmp)/(Q2**2 + 1.)
         endif
         COEFF2_1 = 1. - Y - am*X*Y/(2.*E) - (AMLEP**2)/(4.*E**2)
     &        - (AMLEP**2)/(2.*am*E*(X+1.E-30))
         COEFF2_2 = ( Y**2/2. + (Y*AMLEP**2)/(4.*am*E*(X+1.E-30)) )*
     &        (1. + (4.*am*x**2)/Q2)/(1. + R)
         COEFF2 = COEFF2_1 + COEFF2_2
         
         COR2=-(AMLEP**2)/(4.*am*E*(X+1.E-30))
         
         corw=(AMW2/(AMW2+Q2))**2

         call structg(ip,it,e,x,y,f2,xf3,u,d)

         A=COEFF2*f2
         B=Y*(1-Y/2.+COR2)*xf3
      endif

C      COR1=(AMLEP**2)*Y/(4.*am*E*(X+1.E-30))-am*X*Y/(2.*E)
C      COR1=COR1-(AMLEP**2)/(4.*E*E)-(AMLEP**2)/(2.*am*E*(X+1.E-30))
C      COR2=-(AMLEP**2)/(4.*am*E*(X+1.E-30))

C      Q2=2.D0*am*e*y*x
C      corw=(AMW2/(AMW2+Q2))**2

C      call structg(ip,it,e,x,y,f2,xf3,u,d)

C      A=(1.-Y+Y**2/2.+COR1)*f2
C      B=Y*(1-Y/2.+COR2)*xf3

      IF(IP.EQ. 12.OR.IP.EQ. 14.OR.IP.EQ. 16)
     & DNDXDY=CONST*E*(A+B)*corw
      IF(IP.EQ.-12.OR.IP.EQ.-14.OR.IP.EQ.-16)
     & DNDXDY=CONST*E*(A-B)*corw
*---
      IF(DNDXDY.LE.0.) DNDXDY=0.
C     
      RETURN
      END
