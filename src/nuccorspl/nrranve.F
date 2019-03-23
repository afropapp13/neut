      SUBROUTINE nrrANVE(VR,V,Q1,Q2)
      implicit none
      real*4 VR,V(3),Q1,Q2
      real*4 ct,phi,st

      real*4   ranf
      external ranf
C
C  ASSIGNS TO V(3) A RANDOM VECTOR OF LENGTH VR
C      REAL*4 V(3)
      CT=Q1
      PHI=Q2
      IF(CT.GT.1.)CT=2.*RANF()-1.
      IF(PHI.LT.0.)PHI=6.28318*RANF()
      if (abs(ct).eq.1.) then
         ST = 0.
      else
         ST=SQRT(1.-CT**2)
      endif
      V(1)=VR*ST*COS(PHI)
      V(2)=VR*ST*SIN(PHI)
      V(3)=VR*CT
      RETURN
      END
