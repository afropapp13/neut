      double precision function CROSTO(EPIKIN)

      implicit none
      
      REAL*8 EPIKIN, BRNTFM, XFACT, YFACT
      REAL*8 EPI0XC, EPI0YC, ITEMP, BARNS, EPI0
      INTEGER J, N

      PARAMETER(BRNTFM=1.E2)
      DIMENSION EPI0XC(10),EPI0YC(10),ITEMP(10)
      PARAMETER (XFACT = 17.7/5.0, YFACT = 100./7.05)

      DATA EPI0XC / 
     &  0.00,
     &  0.90,   
     &  2.05,   
     &  2.42,   
     &  2.85,   
     &  3.50,   
     &  4.22,
     &  5.30,   
     &  6.98,   
     &  17.7 /

      DATA EPI0YC /
     &  0.00,
     &  7.72,
     &  1.85,
     &  2.10,
     &  1.95,
     &  2.70,   
     &  2.20,   
     &  2.55,   
     &  2.15,   
     &  1.85 /
      
      EPI0 = EPIKIN*XFACT

      IF (EPIKIN.GT.4.9) THEN
         CROSTO = 2.4 + 1.2/SQRT(EPIKIN)
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

      DO J =1,9
         IF (ITEMP(J).EQ.-1) THEN
            N = J
            GOTO 10
         ENDIF
      ENDDO

 10   CONTINUE

      BARNS = (EPI0YC(N+1)-EPI0YC(N))/(EPI0XC(N+1)-EPI0XC(N))
      BARNS = BARNS*(EPI0-EPI0XC(N+1))+EPI0YC(N+1)
      
      BARNS = BARNS*YFACT*0.001
      CROSTO = BRNTFM*BARNS
      
      RETURN      
      END
