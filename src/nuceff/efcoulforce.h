      
C     Parameters defined in src/neututils/calc_nuc_dens

      integer nNuclei, nRadPoints, nNucleiPar, nRadPointsPar
      parameter (nNucleiPar=26)
      parameter (nRadPointsPar=150)

      integer nucZ(nNucleiPar)
      real*4 radialPosition(nRadPointsPar)
      real*4 integratedDensity(nRadPointsPar, nNucleiPar)

      common /nucdens/ nNuclei, nRadPoints, nucZ, radialPosition, integratedDensity
