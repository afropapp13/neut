************************************************************************
*     ---------------------------------
      COMPLEX FUNCTION RHOCALC( RES1, J1,
     &                          RES2, J2,
     &                          ALPSH, BETSH, DELSH,
     &                          AM, ETABW)
*     ---------------------------------
C
C     (Purpose)
C       Calculate the rho^{r,s}_{m,m'} for two interfering resonances
C       r,s  = two resonance (r = s is allowed), e.g. Delta and S11(1450)
C       m,m' = two isospins to be summed over
C
C       See Rein-Sehgal 1981 p.150 and p.151
C
C     (Input)
C      RES1  :  Number of the first resonance 
C               (see numbering in rsclfm.F and add 1!)
C      AMRES1:  Amplitudes of first resonance
C          J1:  Isospin of first resonance
C      ETABW1:  Eta Breit-Wigner of first resonance
C
C      RES2  :  Number of the second resonance 
C               (see numbering in rsclfm.F and add 1!)
C      AMRES2:  Amplitudes of second resonance
C          J2:  Isospin of second resonance
C      ETABW2:  Eta Breit-Wigner of second resonance
C
C      ALPSH: Alpha as defined in equation A.23 Rein-Sehgal
C      BETSH: Beta as defined in equation A.23 Rein-Sehgal
C      DELSH: Delta as defined in equation A.23 Rein-Sehgal
C
C     (Output)    
C      The rho^{RES1, RES2}_{J1,J2} of two resonances
C
C     (Creation Date and Author)
C       2016.Oct C. Wret: Add in to allow using multiple interfering 
C                         resonances in ejecting the pion in the
C                         resonance rest-frame
C                         See Rein-Sehgal 1981 Eq A20 and onwards
C                         Here we are calculating the "rho" above Eq A20
C
************************************************************************
      IMPLICIT NONE
C     Resonance number (found in rsclfm.F)
      INTEGER RES1,RES2
C     Isospin (J = 3 means 3/2)
      INTEGER J1,J2

      INTEGER NCTR

C     Array that holds the Breit-Wigner etas
      COMPLEX ETABW(31)

      REAL ALPSH,BETSH, DELSH

C     The alpha, beta and delta as defined in equation A.23 Rein-Sehgal

C     Actually only need to pass 2*6 helicity amplitudes
C     Holds the amplitudes
      REAL AM(6,31)
C     AMRES1(1) = AM3C for resonance 1
C     AMRES1(2) = AM1C for resonance 1
C     AMRES1(3) = AP3C for resonance 1
C     AMRES1(4) = AP1C for resonance 1
C     AMRES1(5) = A0PC for resonance 1
C     AMRES1(6) = A0MC for resonance 1
      REAL AM3C(31), AM1C(31), AP3C(31), AP1C(31), A0PC(31), A0MC(31)

C     The idea behind this function is that the way we perform 
C     the summation only depends on the isospin of the resonances, not
C     on their amplitudes or Breit-Wigners

C     Some error checking
      IF (J1.NE.3 .AND. J1.NE.1 .AND. J1.NE.-1 .AND. J1.NE.-3) THEN
        WRITE(*,*) "ERROR, J1 IS NOT 3, 1, -1, -3; J1 = ", J1
        WRITE(*,*) "EXITING"
        STOP
      END IF

      IF (J2.NE.3 .AND. J2.NE.1 .AND. J2.NE.-1 .AND. J2.NE.-3) THEN
        WRITE(*,*) "ERROR, J2 IS NOT 3, 1, -1, -3; J2 = ", J2
        WRITE(*,*) "EXITING"
        STOP
      END IF

      IF (RES1.GT.31 .OR. RES1.LT.1) THEN
        WRITE(*,*) "INVALID RES1 GIVEN; RES1 = ", RES1
        WRITE(*,*) "EXITING"
        STOP
      END IF

      IF (RES2.GT.31 .OR. RES2.LT.1) THEN
        WRITE(*,*) "INVALID RES2 GIVEN; RES2 = ", RES2
        WRITE(*,*) "EXITING"
        STOP
      END IF

C     However, rho^{r,s}_{m,m'} does depend on the amplitudes and the 
C     Breit-Wigner's, so we need to pass those to the function
C     But we don't know which ones to pass until we've got the isospins,
C     So pass all of the resonance amplitudes

      DO 99 NCTR=1,31
      AM3C(NCTR) = AM(1, NCTR)
      AM1C(NCTR) = AM(2, NCTR)
      AP3C(NCTR) = AM(3, NCTR)
      AP1C(NCTR) = AM(4, NCTR)
      A0PC(NCTR) = AM(5, NCTR)
      A0MC(NCTR) = AM(6, NCTR)
 99   CONTINUE     

C     Will always have the Breit-Wigner factor
      RHOCALC = 0.5*CONJG(ETABW(RES1))


C     3/2 xx
      IF (J1 .EQ. 3) THEN

C     3/2 3/2 
        IF (J2 .EQ. 3) THEN
          RHOCALC = RHOCALC*BETSH*BETSH*AM3C(RES1)*AM3C(RES2)
C     3/2 1/2
        ELSE IF (J2 .EQ. 1) THEN
          RHOCALC = RHOCALC*BETSH*DELSH*AM3C(RES1)*A0PC(RES2)
C     3/2 -1/2
        ELSE IF (J2 .EQ. -1) THEN
          RHOCALC = RHOCALC*(-1.)*BETSH*ALPSH*AM3C(RES1)*AP1C(RES2)
C     3/2 -3/2 -- THIS DOESN'T SEEM TO HAPPEN?
C       ELSE IF (J2 .EQ. -3) THEN

        END IF ! end the if on J2

C     1/2 xx
      ELSE IF (J1 .EQ. 1) THEN

C     1/2 3/2
        IF (J2 .EQ. 3) THEN
          RHOCALC = RHOCALC*BETSH*DELSH*A0PC(RES1)*AM3C(RES2)
C     1/2 1/2
        ELSE IF (J2 .EQ. 1) THEN
          RHOCALC = RHOCALC*(BETSH*BETSH*AM1C(RES1)*AM1C(RES2)
     &            + DELSH*DELSH*A0PC(RES1)*A0PC(RES2))
C     1/2 -1/2
        ELSE IF (J2 .EQ. -1) THEN
          RHOCALC = RHOCALC*(BETSH*DELSH*AM1C(RES1)*A0MC(RES2)
     &            - DELSH*ALPSH*A0PC(RES1)*AP1C(RES2))
C     1/2 -3/2
        ELSE IF (J2 .EQ. -3) THEN
          RHOCALC = RHOCALC*(-1.)*BETSH*ALPSH*AM1C(RES1)*AP3C(RES2)

        END IF ! end the if on J2


C     -1/2 xx
      ELSE IF (J1 .EQ. -1) THEN

C     -1/2 3/2
        IF (J2 .EQ. 3) THEN
          RHOCALC = RHOCALC*BETSH*BETSH*AM1C(RES1)*AM3C(RES2)
C     -1/2 1/2
        ELSE IF (J2 .EQ. 1) THEN
          RHOCALC = RHOCALC*(DELSH*BETSH*A0MC(RES1)*AM1C(RES2)
     &            - ALPSH*DELSH*AP1C(RES1)*A0PC(RES2))
C     -1/2 -1/2
        ELSE IF (J2 .EQ. -1) THEN
          RHOCALC = RHOCALC*(ALPSH*ALPSH*AP1C(RES1)*AP1C(RES2)
     &            + DELSH*DELSH*A0MC(RES1)*A0MC(RES2))
C     -1/2 -3/2
        ELSE IF (J2 .EQ. -3) THEN
          RHOCALC = RHOCALC*(-1.)*ALPSH*DELSH*A0MC(RES1)*AP3C(RES2)

        END IF ! end the if on J2


C     -3/2 xx
      ELSE IF (J1 .EQ. -3) THEN

C     -3/2 3/2 THIS DOESN'T SEEM TO HAPPEN
C       IF (J2 .EQ. 3) THEN
C     -3/2 1/2
        IF (J2 .EQ. 1) THEN
          RHOCALC = RHOCALC*(-1.)*ALPSH*BETSH*AP3C(RES1)*AM1C(RES2)
C     -3/2 -1/2
        ELSE IF (J2 .EQ. -1) THEN
          RHOCALC = RHOCALC*(-1.)*ALPSH*DELSH*AP3C(RES1)*A0MC(RES2)
C     -3/2 -3/2
        ELSE IF (J2 .EQ. -3) THEN
          RHOCALC = RHOCALC*ALPSH*ALPSH*AP3C(RES1)*AP3C(RES2)

        END IF ! end the if on J2

      END IF ! end the if on J1

      RHOCALC = RHOCALC*ETABW(RES2)

      RETURN

      END
