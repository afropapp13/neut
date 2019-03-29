*     -----------------------------------------------------
      FUNCTION EFFRMMOM_L(RADIUS,PMOM)
*     -----------------------------------------------------
*
*     ( purpose )
*        PROBABILITY TO NUCLEON TO TAKE FERMI MOMENTUM PF for Local Fermigas
*
*     ( input )
*        RADIUS  : radius   (in fermi)
*        PMOM    : momentum (in MeV)
*
*     ( output )
*       effrmmom_L: probability to take PF
*
      REAL*4 EFFRMMOM_L

      REAL*4 DUM1,DUM2,R

      REAL*4   EFFRMGAS
      external EFFRMGAS

      REAL*4 PFMAX

      R = RADIUS
      PFMAX = EFFRMGAS(DUM1,DUM2,R)

      if ((PMOM.GT.PFMAX).OR.(PFMAX.LE.0)) THEN
         EFFRMMOM_L = 0.
      else
         EFFRMMOM_L = (3*(PMOM**2))/(PFMAX**3)
      endif

      RETURN
      END
