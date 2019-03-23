      double precision function CROSIN(EPIKIN)

      implicit none

      REAL*8 EPIKIN, BRNTFM, XFACT, YFACT
      REAL*8 EPI0XC, EPI0YC, ITEMP, BARNS, EPI0
      INTEGER J, N

      PARAMETER(BRNTFM=1.D2)
      DIMENSION EPI0XC(11),EPI0YC(11),ITEMP(11)
      PARAMETER (XFACT = 17.8/5.0, YFACT = 60./8.42)

      DATA EPI0XC / 
     &  0.00,
     &  0.72,   
     &  1.08,   
     &  1.48,   
     &  2.15,   
     &  3.58,   
     &  4.62,   
     &  5.30,   
     &  6.43,   
     &  7.10,   
     &  17.8 /

      DATA EPI0YC / 
     &  0.00,
     &  3.52,
     &  8.00,
     &  2.35,
     &  1.12,
     &  3.58,   
     &  2.87,   
     &  3.62,   
     &  3.20,   
     &  3.35,
     &  3.10 /
      
      N = -1

      EPI0 = EPIKIN*XFACT

      IF (EPIKIN.GE.4.99) THEN
         CROSIN = 2. +.2*SQRT(5./EPIKIN)
         RETURN
      ENDIF

      ITEMP(1)  = SIGN(1.D0, EPI0 - EPI0XC(2) )
      ITEMP(2)  = SIGN(1.D0, EPI0 - EPI0XC(3) )
      ITEMP(3)  = SIGN(1.D0, EPI0 - EPI0XC(4) )
      ITEMP(4)  = SIGN(1.D0, EPI0 - EPI0XC(5) )
      ITEMP(5)  = SIGN(1.D0, EPI0 - EPI0XC(6) )
      ITEMP(6)  = SIGN(1.D0, EPI0 - EPI0XC(7) )
      ITEMP(7)  = SIGN(1.D0, EPI0 - EPI0XC(8) )
      ITEMP(8)  = SIGN(1.D0, EPI0 - EPI0XC(9) )
      ITEMP(9)  = SIGN(1.D0, EPI0 - EPI0XC(10))
      ITEMP(10) = SIGN(1.D0, EPI0 - EPI0XC(11))

      DO J =1,10
         IF (ITEMP(J).EQ.-1) THEN
            N = J
            GOTO 10
         ENDIF
      ENDDO

 10   CONTINUE

      BARNS = (EPI0YC(N+1)-EPI0YC(N))/(EPI0XC(N+1)-EPI0XC(N))
      BARNS = BARNS*(EPI0-EPI0XC(N+1))+EPI0YC(N+1)

      BARNS = BARNS*YFACT*0.001   !CONVERT FROM M-BARNS TO BARNS
      CROSIN = BARNS*BRNTFM

      RETURN      
      END
