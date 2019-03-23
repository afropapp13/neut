*     -----------------------------------------------------
      COMPLEX FUNCTION EFS31(Q)
*     -----------------------------------------------------
*     
*     ( purpose )
*       OBTAIN PHASE SHIFT OF S31 WAVE IN PI-N SCATTERING
*
*     ( input ) 
*           Q   : Relative momentum in C-M system(MeV)
*
*     ( output ) 
*        EFS31  : (EXP(2*I*DELTA)-1)/2I/Q**3
C
      COMPLEX I
      I=(0.0,1.0)
C  PUT VALUE
      AMPI=139.0
      AMN=939.0
      B=-0.112/AMPI
      C=-0.0307/AMPI**3
      D=0.0021/AMPI**5
      X=0.31
      GAMMA=170.0
      OMG0=1655.0
      Q0=550.0
C
      OMG=SQRT(Q**2+AMPI**2)+SQRT(Q**2+AMN**2)
      IF(OMG.NE.OMG0) GO TO 1
      EFS31=I/Q
      GO TO 2
    1 TAND=Q*(B+C*Q*Q+D*Q**4+X*GAMMA*OMG0/Q0/(OMG0*OMG0-OMG*OMG))
      EFS31=(1.0+I*TAND)/(1.0-I*TAND)
      EFS31=-(EFS31-1.0)/2.0*I/Q
    2 RETURN
      END
