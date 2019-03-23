*************************************************
*     -------------------------------------------
      DOUBLE PRECISION FUNCTION FNRWORLD(x,q2)
*     -------------------------------------------
*
*     (Purpose)
*     Calcurate Rworld function used in the relation
*     F1 and F2
*     (See Bodek paper:hep-ex/0203008)
*
*     (Input)
*     x  : Q**2/(2Mv) BJORKEN x
*     q2 : Q**2
*
*     (Output)
*     Rworld : Longitudinal R
*     
*     (Creation and Date)
*     G.Mitsuka, Aug'06
*************************************************
      implicit none

      real*8 x, q2, theta, r1, r2

      theta = 1. + 12.*Q2/(Q2 + 1.)*(0.125*0.125)/(0.125*0.125 + x**2)
      r1 = 0.0635/log(Q2/0.04)*theta
      r2 = 0.5747/Q2 - 0.3534/(Q2**2 + 0.09)
      fnrworld = r1 + r2

      return
      end
