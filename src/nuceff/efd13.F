*     -----------------------------------------------------
      COMPLEX FUNCTION EFD13(Q)      
*     -----------------------------------------------------
*     
*     ( purpose )
*       OBTAIN PHASE SHIFT OF D13 WAVE IN PI-N SCATTERING
*
*     ( input ) 
*          Q   : RELATIVE MOMENTUM IN C-M SYSTEM
*
*     ( output ) 
*         EFD13:  (EXP(2*I*DELTA)-1)/2I/Q


      COMPLEX I
      I=(0.0,1.0)
C     PUT VALUE
      AMPI=139.0
      AMN=939.0
      BB=0.109E-2*(Q/AMPI)**5
      CC=-0.031E-3*(Q/AMPI)**7
      DD=0.003E-4*(Q/AMPI)**9
      X=0.54
      GAMMA=125.0
      OMG0=1525.0
      Q0=459.0
C     
      OMG=SQRT(Q**2+AMPI**2)+SQRT(Q**2+AMN**2)
      IF(OMG.NE.OMG0) GO TO 1
      EFD13=I/Q
      GO TO 2
C     
    1 TAND=BB+CC+DD+X*GAMMA*OMG0/(OMG0*OMG0-OMG*OMG)*(Q/Q0)**5
      EFD13=(1.0+I*TAND)/(1.0-I*TAND)
      EFD13=-(EFD13-1.0)/2.0*I/Q
    2 RETURN
      END
