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
*     ????.??.?? ; Previous versions
*     2006.08.?? ; G.Mitsuka
*                  Add Bodek-Yang corrections 
*     2016.02.29 ; C. Bronner
*                  Fix definition of x when Bodek-Yang corrections are used.
*                  If Bodek-Yang is not used, keep the target mass corrections
*     2016.07.20 ; C. Bronner
*                  Target mass corrections moved to PDF calculation (grv98_lo.F)
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

C     Assign lepton mass based on neutrino and interaction type
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

C     Assign target nucleon mass
      if(it.eq.2212) then
        am = xmp
        const = consp
      elseif(it.eq.2112) then
        am = xmn
        const = consn
      else
        return
      endif

C     Check that the value of y leaves enough energy to produce the lepton
      YMAX=1.-AMLEP/E
      IF(YMAX.LT.0.) THEN
*      PRINT *,' ***** 1-ML/E < 0.******',X,Y
      RETURN
      END IF

      IF(Y.GT.YMAX) THEN
*      PRINT *,' ***** Y > 1-ML/E ****** X,Y,YMAX=',X,Y,YMAX
      RETURN
      END IF

C     Check value of y is not too large
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
      XMAX=QSMX/(2.*AM*Y*E)
      IF(X.GE.XMAX) THEN
*      PRINT *,' ***** X > Q2MAX/(2MYE) *****',X,Y
      RETURN
      END IF

      
      XMAX=1      
      IF(X.GE.XMAX) THEN
*      PRINT *,' ****** X>1 ******',X,Y
      RETURN
      END IF

C     Check those values of x and y give meaningful lepton kinematics
      Q2=X*(2.*AM*Y*E)
      
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

C     Compute double differential cross section in x and y.
C     Formula taken from C.H. Albright and C. Jarlskog, Nucl. Phys B84, 467 (1975), modified as this paper looks at the case of heavy lepton.
      Q2=2.D0*am*e*y*x
      if(nebodek.eq.0) then
!     --> Original
*        Factors in front of the structure functions
         COR1=(AMLEP**2)*Y/(4.*am*E*(X+1.E-30))-am*X*Y/(2.*E)
         COR1=COR1-(AMLEP**2)/(4.*E*E)-(AMLEP**2)/(2.*am*E*(X+1.E-30))
         COR2=-(AMLEP**2)/(4.*am*E*(X+1.E-30))
         
         corw=(AMW2/(AMW2+Q2))**2
*        Get structure functions
         call structg(ip,it,itype,e,x,y,f2,xf3)

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
     &        (1. + (4.*(am**2)*x**2)/Q2)/(1. + R)

         COEFF2 = COEFF2_1 + COEFF2_2
         
         COR2=-(AMLEP**2)/(4.*am*E*(X+1.E-30))
         
         corw=(AMW2/(AMW2+Q2))**2
*        Get structure functions
         call structg(ip,it,itype,e,x,y,f2,xf3)

         A=COEFF2*f2
         B=Y*(1-Y/2.+COR2)*xf3
      endif


      IF(IP.EQ. 12.OR.IP.EQ. 14.OR.IP.EQ. 16)
     & DNDXDY=CONST*E*(A+B)*corw
      IF(IP.EQ.-12.OR.IP.EQ.-14.OR.IP.EQ.-16)
     & DNDXDY=CONST*E*(A-B)*corw
*---
      IF(DNDXDY.LE.0.) DNDXDY=0.
C     
      RETURN
      END
