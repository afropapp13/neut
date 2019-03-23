C     Calculate Q2 limits for given s=2*mi*Enu+mi^2; Enu,W = const      
      SUBROUTINE RES_Q2_limits(s,W,ml,Q2_min,Q2_max)

      implicit none

      Real*4 s, W, ml, Q2_min, Q2_max
      Real*4 ml2, mi, mi2
      Real*4 RES_Lambda, C, D
      EXTERNAL RES_Lambda      

      ml2 = ml**2
      mi  = 0.94
      mi2 = mi**2

      C = 0.5*(W**2+ml2-mi2*(W**2-ml2)/s)
      D = RES_Lambda(1.,ml2/s,(W**2)/s)
      Q2_min = 0.5*(s-mi2)*(1.-SQRT(D))-C
      Q2_max = 0.5*(s-mi2)*(1.+SQRT(D))-C

      END
