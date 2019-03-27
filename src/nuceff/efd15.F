*     -----------------------------------------------------
      COMPLEX FUNCTION EFD15(Q)
*     -----------------------------------------------------
*
*     ( purpose )
*       OBTAIN PHASE SHIFT OF D15 WAVE IN PI-N SCATTERING
*
*     ( input )
*        Q      : relative momentum in C-M system
*
*     ( output )
*        EFD15  : (EXP(2*I*DELTA)-1)/2I/Q

      COMPLEX I
      I=(0.0,1.0)
C  PUT VALUE
      AMPI=139.0
      AMN=939.0
      BB=0.112E-2*(Q/AMPI)**5
      CC=-0.270E-3*(Q/AMPI)**7
      DD=0.19E-4*(Q/AMPI)**9
      X=0.43
      GAMMA=155.0
      OMG0=1670.0
      Q0=560.0
C
      OMG=SQRT(Q**2+AMPI**2)+SQRT(Q**2+AMN**2)
      IF(OMG.NE.OMG0) GO TO 1
      EFD15=I/Q
      GO TO 2
C
    1 TAND=BB+CC+DD+X*GAMMA*OMG0/(OMG0*OMG0-OMG*OMG)*(Q/Q0)**5
      EFD15=(1.0+I*TAND)/(1.0-I*TAND)
      EFD15=-(EFD15-1.0)/2.0*I/Q
    2 RETURN
      END
