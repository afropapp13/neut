*************************************************
*     -------------------------------------------
      DOUBLE PRECISION FUNCTION FNR1998(x,q2)
*     -------------------------------------------
*
*     (Purpose)
*     Calcurate R1998 function used in the relation
*     F1 and F2
*     (See E143 paper:arXiv:hep-ex/9808028v1 26 Aug 1998)
*
*     (Input)
*     x  : Q**2/(2Mv) BJORKEN x
*     q2 : Q**2
*
*     (Output)
*     R1998 : Longitudinal R
*     
*     (Creation and Date)
*     J.Xia, Jan'20
*************************************************
      implicit none

      real*8 x, q2, theta, ra, rb, rc

      theta = 1. + 12.*Q2/(Q2 + 1.)*(0.125*0.125)/(0.125*0.125 + x**2)
      ra = 0.0485/LOG(Q2/0.04)*theta+
     &     0.5470/((Q2**4.+2.0621**4)**0.25)*(1-0.3804*x+0.5090*x**2)
     &     *x**(-0.0285)
      rb = 0.0481/LOG(Q2/0.04)*theta+
     &     (0.6114/Q2-0.3509/(Q2**2.+0.3**2))*(1-0.4611*x+0.7172*x**2)
     &     *x**(-0.0317)
      rc = 0.0577/LOG(Q2/0.04)*theta+
     &     0.4644*((Q2-12.3708*x+43.1043*x**2-41.7415*x**3)**2
     &     +1.8288**2)**(-0.5)


      fnr1998 = (ra+rb+rc)/3.D0

      return
      end
