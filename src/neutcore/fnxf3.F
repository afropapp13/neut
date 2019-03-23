**********************************************************************
*     ------------------
      FUNCTION FNXF3(X,Y,E)
*     ------------------
*     
*     ( purpose )
*       FUNCTION TO GIVE XF3(X,Y,E)
*
*     ( input ) 
*     X : BJORKEN SCALING PARAMETER
*     E : ENERGY OF NEUTRINO
*
*     ( output ) 
*       FNXF3 : XF3(X,E)
*
*     ( creation date and author )
*       1983.11.13 ; First version by M.Nakahata
*       1995.03.03 ; K.Kaneyuki for S.K.
*                    XF3 -> FNXF3
*
**********************************************************************
      implicit none

      real xdat(11),axf3(11),bxf3(11)
      integer i,j,k,l
      real xx,x,y,E,Q
      real alpha,beta,x1,x2
      real f1,f2,pmass,fnxf3
      data  pmass/0.938/
      data  xdat/0.015,0.045,0.080,0.125,0.175,0.225,0.275,
     &   0.350,0.450,0.550,0.650/
      data  axf3/0.02,0.06,0.04,0.02,-0.04,0.02,-0.1,-0.1,
     &  -0.12,-0.26,-0.46/
      data  bxf3/-1.32,-0.58,-0.44,-0.22,-0.08,-0.26,0.04,
     &  -0.22,-0.56,-0.66,-0.66/
 
      xx = x
      Q = 2*PMASS*E*x*y
      
      if (xx.le.xdat(1)) then 
         x1 = xdat(1)
         x2 = xdat(2)
         alpha = xx - x1
         beta  = x2 - xx
         k = 1
       end if
 
       if (xx.gt.xdat(11)) then
         x1 = xdat(10)
         x2 = xdat(11)
         alpha = xx - x1
         beta  = x2 -xx
         k = 10
       endif
   
       do i = 1,10
        if (((xx-xdat(i))*(xx-xdat(i+1))).le.0)  then
          x1 = xdat(i)
          x2 = xdat(i+1)
          alpha = xx - x1
          beta =  x2 - xx
          k = i
        endif
       enddo

       f1 = exp(bxf3(k))*(q**axf3(k))
       f2 = exp(bxf3(k+1))*(q**axf3(k+1))
      
       fnxf3 = (beta*f1+alpha*f2)/(alpha+beta)

       return 
       end  









