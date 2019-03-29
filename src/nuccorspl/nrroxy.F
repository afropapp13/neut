          real FUNCTION nrROXY(R)
*
*     ( Creation Date and Author )
*       2007.11.05 ; G.Mitsuka - parameters are given by nrsettarget
*
          implicit none

#include "nrnuclparam.h"
          real*4 R
          real*4 w,c,z,d,g,e
          real*4 gr,RO,dro

          parameter (d=.102, g=2.76, e=.35)

          w=NRWPARM
          c=NRC
          z=NRAF
          RO=(1+W*(R/C)**2)/(1+EXP((R-C)/Z))
          IF(R.EQ.0)THEN
                 DRO=1
          ELSE
                 GR=G*R
                 DRO=SIN(GR)/GR
          END IF
          DRO=DRO*D*EXP(-(E*R)**2)
          nrROXY=RO-DRO
          if(nrROXY.ge.1.) nrROXY=0.9999999
          if(nrROXY.le.0.) nrROXY=1.E-7
          RETURN
          END
