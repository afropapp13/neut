*     -----------------------------------------------------
      COMPLEX FUNCTION EFP11(Q)
*     -----------------------------------------------------
*
*     ( purpose )
*       OBTAIN PHASE SHIFT OF P11 WAVE IN PI-N SCATTERING
*
*     ( input )
*           Q   : Relative momentum in C-M system(MeV)
*
*     ( output )
*        EFP11  : (EXP(2*I*DELTA)-1)/2I/Q**3

      COMPLEX I
      I=(0.0,1.0)
C  PUT VALUE
      AMPI=139.0
      AMN=939.0
      B=-0.0571/AMPI**3
      C=0.0254/AMPI**5
      D=-0.0029/AMPI**7
      X=0.61
      GAMMA=230.0
      OMG0=1435.0
      Q0=393.0
C
      OMG=SQRT(Q**2+AMPI**2)+SQRT(Q**2+AMN**2)
      IF(OMG.NE.OMG0) GO TO 1
      EFP11=I/Q**3
      GO TO 2
    1 TAND=Q**3*(B+C*Q*Q+D*Q**4+X*GAMMA*OMG0/Q0**3/(OMG0*OMG0-OMG*OMG))
      EFP11=(1.0+I*TAND)/(1.0-I*TAND)
      EFP11=-(EFP11-1.0)/2.0*I/Q**3
    2 RETURN
      END
