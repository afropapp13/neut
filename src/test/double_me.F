      DOUBLE PRECISION FUNCTION double_me(argin)
        DOUBLE PRECISION argin
        double_me = argin * 2
        RETURN
      END

      REAL*4 FUNCTION double_me_f(argin)
        REAL*4 argin
        double_me_f = argin * 2
        RETURN
      END

      REAL FUNCTION double_me_auto(argin)
        REAL argin
        double_me_auto = argin * 2
        RETURN
      END


      SUBROUTINE logical_test(argin)
        LOGICAL argin
        WRITE(*,*) "FORT: ", argin
      END
