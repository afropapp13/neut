C************************************************************************
C*     ----------------------------------
        SUBROUTINE NECROSEL( EPIKIN, TVALUE, DSIGEL)
C*     ----------------------------------
C*
C*     (Purpose)
C*       THIS SUBROUTINE IS FOR 'COHVCT' THAT IS COHERENT PI0 PRODUCTION
C*       PROGRAM.
C*
C*     (Input)
C*       EPIKIN : pion kinetic energy 
C*       TVALUE : value of |t| 
C*       TMIN   : minimum allowed value for |t|
C*
C*     (Output)
C*       DSIGEL : elastic differential cross section as described in the
C*                Berger and Sehgal paper
C*     
C*     (Creation Date and Author)
C*       2015.03.15 ; Paul Martins
C*
C*     (Comment)
C* ELASTIC PI-0 / NUCLEON SCATTERING, OBTAINED BY LINE FITTING THE DATA GIVEN
C* IN THE PAPER BY BERGER AND SEGHAL.
C*
C************************************************************************
      IMPLICIT NONE

#include "necard.h"
#include "neutmodel.h"
      
      REAL*4 EPIKIN
      REAL*4 TVALUE,DSIGEL
      INTEGER*4 i

      REAL*4  EPI0(12),EPI0A1(12),EPI0b1(12),A1,b1

C These are the values in the table of the Berger and Sehgal paper
      DATA EPI0 /0.00, 0.076, 0.080, 0.100, 0.148, 0.162, 0.226, 0.486,
     &  0.584, 0.662, 0.776, 0.870/

      DATA EPI0A1 /0.00,11600,14700,18300,21300,22400,16400,5730,4610,
     &  4570, 4930, 5140/

      DATA EPI0b1 /116, 116, 109.0, 89.8, 91.0, 89.2, 80.8, 54.6, 55.2,
     &  58.4, 60.5, 62.2/

C Set the value of elastic cross section according to Berger and Sehgal
C paper, as function of pion kinetic energy (linear interpolation)
 50   IF (EPIKIN.LT.0.0) goto 100
      IF (EPIKIN.LT.EPI0(2)) THEN
        A1 = (EPI0A1(2)-EPI0A1(1))/(EPI0(2)-EPI0(1))
     &  *(EPIKIN-EPI0(2))+EPI0A1(2)
        b1 = (EPI0b1(2)-EPI0b1(1))/(EPI0(2)-EPI0(1))
     &  *(EPIKIN-EPI0(2))+EPI0b1(2)
        DSIGEL = A1*(1+fA1COH)*EXP(-1*b1*(1+fb1COH)*TVALUE)
        GO TO 101
      ELSE IF (EPIKIN.LT.EPI0(3)) THEN
        A1 = (EPI0A1(3)-EPI0A1(2))/(EPI0(3)-EPI0(2))
     &  *(EPIKIN-EPI0(3))+EPI0A1(3)
        b1 = (EPI0b1(3)-EPI0b1(2))/(EPI0(3)-EPI0(2))
     &  *(EPIKIN-EPI0(3))+EPI0b1(3)
        DSIGEL = A1*(1+fA1COH)*EXP(-1*b1*(1+fb1COH)*TVALUE)
        GO TO 101     
      ELSE IF (EPIKIN.LT.EPI0(4)) THEN
        A1 = (EPI0A1(4)-EPI0A1(3))/(EPI0(4)-EPI0(3))
     &  *(EPIKIN-EPI0(4))+EPI0A1(4)
        b1 = (EPI0b1(4)-EPI0b1(3))/(EPI0(4)-EPI0(3))
     &  *(EPIKIN-EPI0(4))+EPI0b1(4)
        DSIGEL = A1*(1+fA1COH)*EXP(-1*b1*(1+fb1COH)*TVALUE)
        GO TO 101  
      ELSE IF (EPIKIN.LT.EPI0(5)) THEN
        A1 = (EPI0A1(5)-EPI0A1(4))/(EPI0(5)-EPI0(4))
     &  *(EPIKIN-EPI0(5))+EPI0A1(5)
        b1 = (EPI0b1(5)-EPI0b1(4))/(EPI0(5)-EPI0(4))
     &  *(EPIKIN-EPI0(5))+EPI0b1(5)
        DSIGEL = A1*(1+fA1COH)*EXP(-1*b1*(1+fb1COH)*TVALUE)
        GO TO 101
      ELSE IF (EPIKIN.LT.EPI0(6)) THEN
        A1 = (EPI0A1(6)-EPI0A1(5))/(EPI0(6)-EPI0(5))
     &  *(EPIKIN-EPI0(6))+EPI0A1(6)
        b1 = (EPI0b1(6)-EPI0b1(5))/(EPI0(6)-EPI0(5))
     &  *(EPIKIN-EPI0(6))+EPI0b1(6)
        DSIGEL = A1*(1+fA1COH)*EXP(-1*b1*(1+fb1COH)*TVALUE)
        GO TO 101
      ELSE IF (EPIKIN.LT.EPI0(7)) THEN
        A1 = (EPI0A1(7)-EPI0A1(6))/(EPI0(7)-EPI0(6))
     &  *(EPIKIN-EPI0(7))+EPI0A1(7)
        b1 = (EPI0b1(7)-EPI0b1(6))/(EPI0(7)-EPI0(6))
     &  *(EPIKIN-EPI0(7))+EPI0b1(7)
        DSIGEL = A1*(1+fA1COH)*EXP(-1*b1*(1+fb1COH)*TVALUE)
        GO TO 101
      ELSE IF (EPIKIN.LT.EPI0(8)) THEN
        A1 = (EPI0A1(8)-EPI0A1(7))/(EPI0(8)-EPI0(7))
     &  *(EPIKIN-EPI0(8))+EPI0A1(8)
        b1 = (EPI0b1(8)-EPI0b1(7))/(EPI0(8)-EPI0(7))
     &  *(EPIKIN-EPI0(8))+EPI0b1(8)
        DSIGEL = A1*(1+fA1COH)*EXP(-1*b1*(1+fb1COH)*TVALUE)
        GO TO 101
      ELSE IF (EPIKIN.LT.EPI0(9)) THEN
        A1 = (EPI0A1(9)-EPI0A1(8))/(EPI0(9)-EPI0(8))
     &  *(EPIKIN-EPI0(9))+EPI0A1(9)
        b1 = (EPI0b1(9)-EPI0b1(8))/(EPI0(9)-EPI0(8))
     &  *(EPIKIN-EPI0(9))+EPI0b1(9)
        DSIGEL = A1*(1+fA1COH)*EXP(-1*b1*(1+fb1COH)*TVALUE)
        GO TO 101
      ELSE IF (EPIKIN.LT.EPI0(10)) THEN
        A1 = (EPI0A1(10)-EPI0A1(9))/(EPI0(10)-EPI0(9))
     &  *(EPIKIN-EPI0(10))+EPI0A1(10)
        b1 = (EPI0b1(10)-EPI0b1(9))/(EPI0(10)-EPI0(9))
     &  *(EPIKIN-EPI0(10))+EPI0b1(10)
        DSIGEL = A1*(1+fA1COH)*EXP(-1*b1*(1+fb1COH)*TVALUE)
        GO TO 101
      ELSE IF (EPIKIN.LT.EPI0(11)) THEN
        A1 = (EPI0A1(11)-EPI0A1(10))/(EPI0(11)-EPI0(10))
     &  *(EPIKIN-EPI0(11))+EPI0A1(11)
        b1 = (EPI0b1(11)-EPI0b1(10))/(EPI0(11)-EPI0(10))
     &  *(EPIKIN-EPI0(11))+EPI0b1(11)
        DSIGEL = A1*(1+fA1COH)*EXP(-1*b1*(1+fb1COH)*TVALUE)
        GO TO 101
      ELSE IF (EPIKIN.LT.EPI0(12)) THEN
        A1 = (EPI0A1(12)-EPI0A1(11))/(EPI0(12)-EPI0(11))
     &  *(EPIKIN-EPI0(12))+EPI0A1(12)
        b1 = (EPI0b1(12)-EPI0b1(11))/(EPI0(12)-EPI0(11))
     &  *(EPIKIN-EPI0(12))+EPI0b1(12)
        DSIGEL = A1*(1+fA1COH)*EXP(-1*b1*(1+fb1COH)*TVALUE)
        GO TO 101
       ELSE IF (EPIKIN.GE.EPI0(12)) THEN
        A1 = (EPI0A1(12)-EPI0A1(11))/(EPI0(12)-EPI0(11))
     &  *(EPIKIN-EPI0(12))+EPI0A1(12)
        b1 = (EPI0b1(12)-EPI0b1(11))/(EPI0(12)-EPI0(11))
     &  *(EPIKIN-EPI0(12))+EPI0b1(12)
        DSIGEL = A1*(1+fA1COH)*EXP(-1*b1*(1+fb1COH)*TVALUE)
      ENDIF

C      write(*,*) 'fA1 param in necrosel.F = ', fA1COH
C Convert DSIGEL from mbarns to barns
C 101    write(*,*) 'A1 = ', A1,' b1 = ',b1, ' DSIGEL = ',DSIGEL*0.001
 101  DSIGEL = DSIGEL*0.001

 100  RETURN
      END
