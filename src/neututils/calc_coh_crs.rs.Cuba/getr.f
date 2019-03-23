      double precision function GETR(EPIKIN,EPITOT,FEMTO)

      implicit none

      REAL*8 EPIKIN, EPITOT, FEMTO
      REAL*8 PIMASS, PRMASS, HBARC, PI
      REAL*8 FREPOL, FIM, FRE, FIMPOL, RFACTO

      PARAMETER (     PIMASS  = 0.135D0             )
      PARAMETER (     PRMASS  = 0.939D0             )
      PARAMETER (     HBARC   = .197327053D0        )
      PARAMETER (     PI      = 3.141592653589793D0 )

c      FIMPOL = 0.005/2.   Forward Scattering Ampl. Real Pole term.
c      FREPOL = 0.073/2.   Forward Scattering Ampl. Imaginary Pole term. 
      FIMPOL = 0.
      FREPOL = 0.

      IF (EPIKIN.LE.0.0) THEN
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

      GETR = FRE/FIM

      IF (EPIKIN.GT.5.) GETR = 0.1
      IF (FEMTO.eq.0.) GETR = 0.

      RETURN
      END
