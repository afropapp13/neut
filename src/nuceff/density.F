        SUBROUTINE DENSITY(R,dens)
*    PURPOSE ;  CALCULATE NUCLEON DENSITY IN O16
*
*    BASED ON WOOD-SAXON MODEL 
*
*    INPUT  ; R : POSITION (RADIUS )
*    OUTPUT ;   DENS : NUCLEON DENSITY 


       implicit NONE
       real rho,dens,c,aa,R
       INTEGER A,Z

       parameter(aa = 1.80)
       parameter(c = 2.69)
       parameter(rho = 1.80)
       parameter(A = 16)
       parameter(Z = 8)
      
       dens = Z*rho/(1+exp((R-c)/aa))/A

       RETURN
       END
