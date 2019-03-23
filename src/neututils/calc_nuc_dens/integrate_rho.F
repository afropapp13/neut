CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Purpose: Integrand for nuclear density integral in
C              spherical coordinates
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      real*4 function efabrho_rsq(r)
      
      implicit none
      
      real*4   efabrho
      external efabrho

      real*4 r

      efabrho_rsq = efabrho(r)*r*r

      end
