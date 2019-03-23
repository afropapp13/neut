*     -----------------------------------------------------
      COMPLEX FUNCTION EFS11(Q)
*     -----------------------------------------------------
*     
*     ( purpose )
*       OBTAIN PHASE SHIFT OF S11 WAVE IN PI-N SCATTERING
*
*     ( input ) 
*           Q   : Relative momentum in C-M system(MeV)
*
*     ( output ) 
*        EFS11  : (EXP(2*I*DELTA)-1)/2I/Q**3

      COMPLEX I
      I=(0.0,1.0)
C  PUT VALUE
      AMPI=139.0
      AMN=939.0
      B=0.168/AMPI
      C=-0.0354/AMPI**3
      D=0.0027/AMPI**5
      X=0.44
      GAMMA=105.0
      OMG0=1550.0
      Q0=477.0
C
      OMG=SQRT(Q**2+AMPI**2)+SQRT(Q**2+AMN**2)
      IF(OMG.NE.OMG0) GO TO 1
      EFS11=I/Q
      GO TO 2
    1 TAND=Q*(B+C*Q*Q+D*Q**4+X*GAMMA*OMG0/Q0/(OMG0*OMG0-OMG*OMG))
      EFS11=(1.0+I*TAND)/(1.0-I*TAND)
      EFS11=-(EFS11-1.0)/2.0*I/Q
    2 RETURN
      END
