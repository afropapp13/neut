**********************************************************************
*     ------------------
      FUNCTION FNF2(X,Y,E)
*     ------------------
*     
*     ( purpose )
*       FUNCTION TO GIVE F2(X,Y,E)
*
*     ( input ) 
*     X : BJORKEN SCALING PARAMETER
*     E : ENERGY OF NEUTRINO
*
*     ( output ) 
*       FNF2 : F2(X,E)
*
*     ( creation date and author )
*       1983.11.13 ; First version by M.Nakahata
*       1995.03.03 ; K.Kaneyuki for S.K.
*                    F2 -> FNF2
*
**********************************************************************
      implicit none

      real xdat(11),af2(11),bf2(11)
      integer i,j,k,l
      real x,y,E,Q
      real alpha,beta,x1,x2
      real f1,f2,pmass,fnf2
      data  pmass/0.938/
      data  xdat/0.015,0.045,0.080,0.125,0.175,0.225,0.275,
     &   0.350,0.450,0.550,0.650/
      data  af2/0.160,0.120,0.060 ,0.04,0.0,-0.04,-0.04,-0.06,
     &  -0.16,-0.22,-0.28/
      data  bf2/0.160,0.180,0.280,0.200,0.140,0.140,-0.02,-0.26,
     &  -0.40,-0.78,-1.240/ 

      Q = 2*PMASS*E*x*y
      
      if (x.le.xdat(1)) then 
         x1 = xdat(1)
         x2 = xdat(2)
         alpha = x - x1
         beta  = x2 - x
         k = 1
       end if
 
       if (x.gt.xdat(11)) then
         x1 = xdat(10)
         x2 = xdat(11)
         alpha = x - x1
         beta  = x2 -x
         k = 10
       endif
   
       do i = 1,10
        if (((x-xdat(i))*(x-xdat(i+1))).le.0)  then
          x1 = xdat(i)
          x2 = xdat(i+1)
          alpha = x - x1
          beta =  x2 - x
          k = i
        endif
       enddo

       f1 = exp(bf2(k))*(q**af2(k))
       f2 = exp(bf2(k+1))*(q**af2(k+1))
      
       fnf2 = (beta*f1+alpha*f2)/(alpha+beta)

       return 
      end  


