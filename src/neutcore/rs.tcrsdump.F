      PROGRAM CRSSECT
      
      DO 20 IMOD = 1,7
         DO 10 I=1,1000
            E = I * 0.1
            CALL RSTCRST(IMOD,E,SI1,SI2)
            WRITE(120+IMOD,*) (E),SI1,SI2
            WRITE(20+IMOD,*) LOG10(E),SI1,SI2
 10      CONTINUE
 20   CONTINUE
      
      DO 40 IMOD = 11,17
         DO 30 I=1,1000
            E = I * 0.1
            CALL RSTCRST(IMOD,E,SI1,SI2)
            WRITE(120+IMOD,*) (E),SI1,SI2
            WRITE(20+IMOD,*) LOG10(E),SI1,SI2
 30      CONTINUE
 40   CONTINUE
