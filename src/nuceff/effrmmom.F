*     -----------------------------------------------------
      FUNCTION EFFRMMOM(PF)
*     -----------------------------------------------------
*
*     ( purpose )
*        PROBABILITY TO NUCLEON TO TAKE FERMI MOMENTUM PF
*
*     ( input )
*        PF     : fermi momentum of nucleon(MeV)
*
*     ( output )
*       effrmmom: probability to take PF
*

      DIMENSION DISS(10),DISP(12)
      DATA DISS/4. , 3.8 , 3.4 , 2.6 , 1.9 , 1.2 , 0.75 , 0.4 , 0.2,
     &          0./
      DATA DISP/0.8 , 1.6 , 3.4 , 5.9 , 5.5 , 4.9 , 3.8 , 2.5 , 1.5 ,
     &          0.8 , 0.4 , 0./
      DATA SNORM/2.552E6/,PNORM/1.162E7/
C
      EFFRMMOM=0.
      I=IFIX(PF/25.)+1
      FRAC=(PF-25.*FLOAT(I-1))/25.
C -- S WAVE
      IF(PF.GE.225.)GO TO 100
      PSWAVE=DISS(I)+(DISS(I+1)-DISS(I))*FRAC
      EFFRMMOM=PSWAVE*PF**2/SNORM*0.25
C -- P WAVE
  100 IF(PF.GE.275.)RETURN
      PPWAVE=DISP(I)+(DISP(I+1)-DISP(I))*FRAC
      EFFRMMOM=EFFRMMOM+PPWAVE*PF**2/PNORM*0.75
      RETURN
      END
