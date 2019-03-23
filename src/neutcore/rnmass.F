************************************************************************
*     ---------------------------------
      FUNCTION RNMASS(AMP,AM1,AM2,WID1)
*     ---------------------------------
*
*     (Purpose)
*       DETERMINE MASS ,CONSIDERING BREIT-WIGNER FORMULA
*       THE HIGHEST BOUND OF THE MASS IS CONSIDERED BUT THE LOWEST 
*       MASS IS NOT
*     
*     (Input)
*       AMP  ; MASS OF PARENT PARTICLE     ( NON RESONANCE)
*       AM1  ;    "    DECAYING PARTICLE 1 (*****RESONANCE*****)
*       AM2  ;    "    THE PARTNER OF PARTICLE 1
*       WID1 ; FULL WIDTH OF THE MASS OF THE PARTICLE 1
*
*     (Output)
*       RNMASS 
*
*     (Creation Date and Author)
*       1985.10.19 ; T.KAJITA
*       1995.03.03 ; K. KANEYUKI  FOR S.K.
*                    RSMASS -> RNMASS
*     
************************************************************************

    1 RNT=(RLU(DUM)-0.5)*3.14159
      RNMASS=TAN(RNT)*WID1/2. + AM1
C
C -- COSIDER KINEMATICS
C
      IF((RNMASS+AM2).GE.AMP)GO TO 1

      RETURN
      END
