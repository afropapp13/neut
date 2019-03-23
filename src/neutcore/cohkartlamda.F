******************************************
      real function cohkartlamda(r1,r2,r3)
*
*     (Purpose)
*       Return lamda
*
*     (Input)
*       r1,r2,r3
*
*     (Output)
*       cohkartlamda = r1**2 + r2**2 + r3**2 -2*(r1*r2 + r1*r3 + r2*r3)
*
*     (Creation Date and Author)
*       2007.01.08 ; G.Mitsuka -- first creation
*       2007.11.05 ; G.Mitsuka -- split this function to one program
*
      implicit none
      real r1,r2,r3

      cohkartlamda = r1**2+r2**2+r3**2-2*(r1*r2+r2*r3+r3*r1)

      return
      end
