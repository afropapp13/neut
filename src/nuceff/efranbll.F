*     -----------------------------------------------------
      SUBROUTINE EFRANBLL(AA,X)
*     -----------------------------------------------------
*     
*     ( purpose )
*       GENERATE RANDUM NUCLEON POINT IN NUCLEI
*
*     ( input ) 
*       AA     : NUMBER OF NUCLEONS (DUMMY)
*
*     ( output ) 
*       X(3)   : VECTOR (FERMI)
C     

      real*4 AA
      real*4 X(3)

      write(*,*) 'nuceff/efrabll.F was called.'
      write(*,*) 'Obsolete function. Stop.'
      stop

CCC---   
CCC---   
CCC---         DIMENSION X(3)
CCC---   #include<efpion.h>
CCC---         DATA ICALL/0/
CCC---   C     
CCC---         ICALL=ICALL+1
CCC---         IF(ICALL.NE.1)GO TO 200
CCC---   C     
CCC---         RMAX=2.5*C
CCC---         DISMAX=0.
CCC---         DO 100 I=1,50
CCC---            R=RMAX/50.*FLOAT(I)
CCC---            DEN=EFRHO(R)
CCC---            IF(DEN.LT.DISMAX)GO TO 100
CCC---            DISMAX=DEN
CCC---     100 CONTINUE
CCC---   C     
CCC---   C  200 R=RNDM(DUM)*RMAX
CCC---     200 R=RLU(IDUM)*RMAX
CCC---   C      HIT=RNDM(DUM)*DISMAX
CCC---         HIT=RLU(IDUM)*DISMAX
CCC---         DEN=EFRHO(R)
CCC---         IF(HIT.GT.DEN)GO TO 200
CCC---   C     
CCC---         CALL RNDIR(X)
CCC---         X(1)=X(1)*R
CCC---         X(2)=X(2)*R
CCC---         X(3)=X(3)*R
CCC---         RETURN
CCC---         END
CCC---
      
      end
