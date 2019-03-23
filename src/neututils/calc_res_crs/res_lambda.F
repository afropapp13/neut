      REAL FUNCTION RES_Lambda(a,b,c)

      implicit none

      Real*4 a, b, c
      RES_Lambda = a**2+b**2+c**2-2.*(a*b+b*c+a*c)

      RETURN
      END
