C     Calculate Q2 limits for given s=2*mN*Enu+mN^2; Enu,W = const      
C     Not the most trivial thing in the whole world...

      SUBROUTINE RES_Q2_LIMITS(s_r, W_r, ml_r, Q2_min_r, Q2_max_r)

      IMPLICIT NONE

#include "rscons.h"

C     The reals that come from NEUT
      REAL s_r, W_r, ml_r, Q2_min_r, Q2_max_r, xmn_r
      DOUBLE PRECISION s_d, W_d, ml_d, Q2_min_d, Q2_max_d, xmn_d

C     Auxillary variables
C     REAL c, d, res_lambda
      DOUBLE PRECISION c_d, d_d, res_lambda_d
      external res_lambda_d

      xmn_r = xmn

      s_d=DBLE(s_r)
      W_d=DBLE(W_r)
      ml_d=DBLE(ml_r)
      xmn_d=DBLE(xmn_r)

C     write(*,*) "**************"
C     write(*,*) s_d, W_d, ml_d, xmn_d

c     C = 0.5*(W_r**2 + ml_r**2 - (xmn_r**2)*(W_r**2 - ml_r**2)/s_r)
c     D = RES_Lambda(1., (ml_r**2)/s_r, (W_r**2)/s_r)
      C_d = 0.5*(W_d**2 + ml_d**2 - (xmn_d**2)*(W_d**2 - ml_d**2)/s_d)
      D_d = RES_Lambda_d(DBLE(1.), (ml_d**2)/s_d, (W_d**2)/s_d)

C     write(*,*) C_d, D_d

      IF (D_d .LT. 0) THEN
C       write(*,*) "neutcore/res_q2_limits.F:"
C       write(*,*) "D < 0, D = ", D_d
        D_d = 0
      END IF

C     CWRET: changed this to NEGATIVE SIGN
C     to agree with other NEUT QMIN and QMAX calc
C     Q2_min_r = -1.*(0.5*(s_r-XMN_r**2)*(1.-SQRT(D))-C)
C     Q2_max_r = -1.*(0.5*(s_r-XMN_r**2)*(1.+SQRT(D))-C)
      Q2_min_d = -1.*(0.5*(s_d-XMN_d**2)*(1.-SQRT(D_d))-C_d)
      Q2_max_d = -1.*(0.5*(s_d-XMN_d**2)*(1.+SQRT(D_d))-C_d)

      if (q2_min_d .ge. 0) then
        q2_min_d = -1.E-7
      endif
      if (q2_max_d .ge. 0) then
        q2_max_d = -1.E-7
      endif

      q2_min_r = SNGL(q2_min_d)
      q2_max_r = SNGL(q2_max_d)

C     write(*,*) q2_min_d, q2_max_d
C     write(*,*) q2_min_r, q2_max_r

      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Have the RES_Lambda function here
C     a = 1, b = ml^2/s, c = w^2/s
      DOUBLE PRECISION FUNCTION RES_Lambda_d(a,b,c)

      IMPLICIT NONE

      DOUBLE PRECISION a, b, c

C     write(*,*) a,b,c
      RES_Lambda_d = a**2 + b**2 + c**2 - 
     &             2.*(a*b + b*c + a*c)

      RETURN
      END
