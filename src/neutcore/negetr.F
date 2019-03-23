************************************************************************
*     ----------------------------------------
      SUBROUTINE NEGETR(EPIKIN,EPITOT,FEMTO,R)
*     ----------------------------------------
*
*     (Purpose)
*      THIS SUBROUTINE IS FOR 'COHVCT' THAT IS COHERENT PI0 PRODUCTION
*      PROGRAM.
*
*     (Input)
*       EPIKIN :
*       EPITOT :
*       FEMTO  :      
*     
*     (Output)
*       R :
*     
*     (Creation Date and Author)
*       1994.08.30 ; G.DION & K.Okumura
*       1996.02.05 ; K.Kaneyuki   GETR -> NEGETR
*
*     (Comment)
*
************************************************************************

      PARAMETER (     HBARC   = .197327053          )
      PARAMETER (     PI      = 3.141592653589793D0 )

      CALL MCMASSGV(111,PIMASS)
      CALL MCMASSGV(2112,PRMASS)
      
       FIMPOL = 0.005/2.
       FREPOL = 0.073/2.
C      FIMPOL = 0.005/2.   Forward Scattering Ampl. Real Pole term.
C      FREPOL = 0.073/2.   Forward Scattering Ampl. Imaginary Pole term.

      IF (EPIKIN.LT.0.0) THEN
         FRE = 0.0
      ELSE IF (EPIKIN.LT.0.13) THEN
         FRE = FREPOL + (1.1/0.13)*EPIKIN
      ELSE IF (EPIKIN.LT.0.280) THEN
         FRE = 1.100 + (-2.4/0.15)*(EPIKIN-0.13)
      ELSE IF (EPIKIN.LT.0.800) THEN
         FRE = -1.300 + (2./0.73)*(EPIKIN-0.280)
      ELSE IF (EPIKIN.LT.2.4) THEN
         FRE = 0.0 + (-1.5/1.6)*(EPIKIN-0.800)
      ELSE
         FRE = -1.6
      ENDIF

      FRE = FRE/2. + FREPOL

      RFACTO = EPITOT**2 - PIMASS**2
      RFACTO = RFACTO/(EPITOT+((PRMASS**2+PIMASS**2)/(2.*PRMASS)))
      RFACTO = RFACTO*(PRMASS/2.)
      RFACTO = SQRT(ABS(RFACTO))/(4.*PI*HBARC)

      FIM = RFACTO*FEMTO + FIMPOL

      R = FRE/FIM

      IF (EPIKIN.GT.5.) R = 0.1

      RETURN

      END
