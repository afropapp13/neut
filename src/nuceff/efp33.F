*     -----------------------------------------------------
      COMPLEX FUNCTION EFP33(Q)
*     -----------------------------------------------------
*
*     ( purpose )
*       OBTAIN PHASE SHIFT OF P33 WAVE IN PI-N SCATTERING
*
*     ( input )
*           Q   : Relative momentum in C-M system(MeV)
*
*     ( output )
*        EFP33  : (EXP(2*I*DELTA)-1)/2I/Q**3

      COMPLEX I
      I=(0.0,1.0)
C  PUT VALUE
      AMPI=139.0
      AMN=939.0
      B=0.114/AMPI**3
      C=-0.0154/AMPI**5
      D=0.00072/AMPI**7
      X=0.99
      GAMMA=116.0
      OMG0=1233.0
      Q0=228.0
C
      OMG=SQRT(Q**2+AMPI**2)+SQRT(Q**2+AMN**2)
      IF(OMG.NE.OMG0) GO TO 1
      EFP33=I/Q**3
      GO TO 2
    1 TAND=Q**3*(B+C*Q*Q+D*Q**4+X*GAMMA*OMG0/Q0**3/(OMG0*OMG0-OMG*OMG))
      EFP33=(1.0+I*TAND)/(1.0-I*TAND)
      EFP33=-(EFP33-1.0)/2.0*I/Q**3
    2 RETURN
      END
