*     -----------------------------------------------------
      COMPLEX FUNCTION EFP31(Q)
*     -----------------------------------------------------
*
*     ( purpose )
*       OBTAIN PHASE SHIFT OF P31 WAVE IN PI-N SCATTERING
*
*     ( input )
*           Q   : Relative momentum in C-M system(MeV)
*
*     ( output )
*        EFP31  : (EXP(2*I*DELTA)-1)/2I/Q**3
C
      COMPLEX I
      I=(0.0,1.0)
C  PUT VALUE
      AMPI=139.0
      AMN=939.0
      B=-0.0291/AMPI**3
      C=0.00345/AMPI**5
      D=-0.00015/AMPI**7
      X=0.22
      GAMMA=200.0
      OMG0=1850.0
      Q0=678.0
C
      OMG=SQRT(Q**2+AMPI**2)+SQRT(Q**2+AMN**2)
      IF(OMG.NE.OMG0) GO TO 1
      EFP31=I/Q**3
      GO TO 2
    1 TAND=Q**3*(B+C*Q*Q+D*Q**4+X*GAMMA*OMG0/Q0**3/(OMG0*OMG0-OMG*OMG))
      EFP31=(1.0+I*TAND)/(1.0-I*TAND)
      EFP31=-(EFP31-1.0)/2.0*I/Q**3
    2 RETURN
      END
