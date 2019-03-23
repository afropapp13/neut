C   19/10/85 511071035  MEMBER NAME  DCY2BD   (FORT)        FORTRAN77
C
      SUBROUTINE DCY2BD (AM0,AM1,AM2,P)
C
C   GET MOMENTUM OF DECAYED PARTICLE IN C.M.S. ONLY 2 BODY DECAY IS
C  CONSIDERED.
C      INPUT  : AM0 - MASS OF DECAYED PARTICLE
C               AM1 -    "    DECAY PRODUCT 1
C               AM2 -    "          "       2
C      OUTPUT : P   - MOMENTUM
C
C    T.KAJITA   19-OCT-85
C
      IF (AM0 .GT. (AM1+AM2) ) GO TO 10
      WRITE(6,600)AM0,AM1,AM2
 600  FORMAT(' ','** ERROR ** ',F7.1,'<',F7.1,'+',F7.1)
      RETURN
C
 10   P= SQRT( ( AM0**2-(AM1+AM2)**2 )*( AM0**2-(AM1-AM2)**2 )  /
     &                                         (4.*AM0**2)  )
C
      RETURN
      END
